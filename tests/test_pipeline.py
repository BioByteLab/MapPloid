from genmap_fixture import mock_genmap
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import csv
import json
import shutil
import unittest

import pysam
import mapploid.pipeline as p


@unittest.skipUnless(shutil.which('samtools'), 'samtools required')
class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.ref=self.root/'ref.fa'
        self.ref.write_text('>chr1\n'+'A'*2000+'\n')
        self.bg=self.root/'scores.bedgraph'
        self.bg.write_text('chr1\t0\t1951\t1\n')
        self.enterContext(mock_genmap(self.bg))
        self.assets=self.root/'assets'
        self.bam=self.root/'sample.bam'
        with pysam.AlignmentFile(str(self.bam),'wb',header={
            'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':2000}]}):
            pass
        pysam.index(str(self.bam))

    def tearDown(self):
        self.temp.cleanup()

    def prepare(self):
        return p.prepare(p.PrepareConfig(self.ref,self.assets,50,1,min_region_length=300))

    def result(self, output):
        with output.open() as h:
            return next(csv.DictReader(h,delimiter='\t'))

    def test_small_reference_reuses_assets_without_touching_source_fasta_index(self):
        self.prepare()
        assets=p.load_assets(self.assets)
        self.assertEqual(len(assets['stages']),1)
        self.assertLess(assets['stages'][0]['bases'],20_000_000)
        self.assertFalse(Path(str(self.ref)+'.fai').exists())
        before=(self.assets/'assets.json').stat().st_mtime_ns
        self.prepare()
        self.assertEqual(before,(self.assets/'assets.json').stat().st_mtime_ns)

    def test_zero_coverage_is_uncertain_and_completed_run_is_reused(self):
        self.prepare()
        config=p.CallConfig(self.bam,self.assets,self.root/'out')
        result=p.call(config)
        row=self.result(result)
        self.assertEqual(row['status'],'UNCERTAIN')
        self.assertEqual(p.read_json(result.parent/'details/run.json')['result']['reason'],'NO_COVERAGE')
        modified=result.stat().st_mtime_ns
        self.assertEqual(p.call(config),result)
        self.assertEqual(modified,result.stat().st_mtime_ns)
        result.write_text('damaged')
        with self.assertRaisesRegex(ValueError,'modified'):
            p.call(config)

    def test_no_mappable_regions_is_valid_uncertain_result(self):
        self.bg.write_text('chr1\t0\t1951\t0.5\n')
        self.prepare()
        result=p.call(p.CallConfig(self.bam,self.assets,self.root/'out'))
        self.assertEqual(p.read_json(result.parent/'details/run.json')['result']['reason'],'NO_ELIGIBLE_REGIONS')
        self.assertEqual((self.root/'out/regions.bed').read_text(),'')

    def test_asset_content_change_is_rejected(self):
        self.prepare()
        manifest=p.read_json(self.assets/'assets.json')
        (self.assets/manifest['stages'][0]['bed']).write_text('chr1\t0\t10\n')
        with self.assertRaisesRegex(ValueError,'modified'):
            p.call(p.CallConfig(self.bam,self.assets,self.root/'out'))

    def test_public_call_preserves_v20_short_read_validation(self):
        from mapploid.cli import main
        self.prepare()
        asset_stamp=(self.assets/'assets.json').stat().st_mtime_ns
        header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':2000}]}
        with pysam.AlignmentFile(str(self.bam),'wb',header=header) as handle:
            for length in (20,75):
                read=pysam.AlignedSegment()
                read.query_name=f'read_{length}';read.query_sequence='A'*length
                read.flag=0;read.reference_id=0;read.reference_start=100
                read.mapping_quality=60;read.cigar=[(0,length)]
                read.query_qualities=pysam.qualitystring_to_array('I'*length)
                handle.write(read)
        pysam.index(str(self.bam))
        output=self.root/'out'
        status=main(['call','--bam',str(self.bam),'--assets',str(self.assets),
                     '--out',str(output),'--sample','short_mixed'])
        # The manuscript v20 source retains its existing k-mer/read-length check.
        self.assertEqual(status,1)
        self.assertEqual(self.result(output/'result.tsv')['status'],'FAILED')
        self.assertIn('shorter',p.read_json(output/'details/run.json')['result']['reason'])
        self.assertEqual((self.assets/'assets.json').stat().st_mtime_ns,asset_stamp)

    def test_tool_failure_is_failed_not_biological_uncertainty(self):
        self.prepare()
        with self.assertRaisesRegex(RuntimeError,'samtools'):
            p.call(p.CallConfig(self.bam,self.assets,self.root/'out',samtools_executable='no-such-samtools'))
        self.assertEqual(p.read_json(self.root/'out/details/run.json')['state'],'FAILED')
        self.assertEqual(self.result(self.root/'out/result.tsv')['status'], 'FAILED')
        self.assertTrue((self.root/'out/report.html').is_file())

    def test_changed_input_cannot_reuse_previous_call(self):
        self.prepare()
        config=p.CallConfig(self.bam,self.assets,self.root/'out')
        p.call(config)
        import os
        info=self.bam.stat()
        os.utime(self.bam,ns=(info.st_atime_ns,info.st_mtime_ns+1_000_000))
        with self.assertRaisesRegex(ValueError,'inputs or settings changed'):
            p.call(config)

    def test_expansion_early_stop_and_resume_skip_completed_stage(self):
        self.prepare()
        manifest=p.read_json(self.assets/'assets.json')
        manifest['stages']=[]
        for i,stop in enumerate((500,1000,1500),1):
            path=self.assets/f'stage{i}.bed'
            path.write_text(f'chr1\t50\t{stop}\n')
            manifest['sha256'][path.name]=p.sha256(path)
            manifest['stages'].append(dict(name=str(i),bed=path.name,bases=stop-50))
        p.atomic_json(self.assets/'assets.json',manifest)
        counts_calls=[]
        model_calls=[]
        def counts(config):
            counts_calls.append(config.output_dir.parent.name)
            config.output_dir.mkdir(parents=True,exist_ok=True)
            path=config.output_dir/'run_manifest.json'
            p.atomic_json(path,dict(filtering=dict(retained_sites=30000),extraction=dict(acgt_observations=100000)))
            data=config.output_dir/'counts.tsv'
            data.write_text('fixture')
            return dict(manifest=path,biallelic_sites=data)
        def model(config):
            model_calls.append(config.output_dir.parent.name)
            if len(model_calls)==2:
                raise RuntimeError('interrupted fixture')
            config.output_dir.mkdir(parents=True,exist_ok=True)
            row=p.empty_evidence('LOW_DELTA_BIC')
            row['warnings']='LOW_BOOTSTRAP_SUPPORT;INSUFFICIENT_CHROMOSOME_EVIDENCE;CHROMOSOME_HETEROGENEITY'
            if config.output_dir.parent.name=='stage_2':
                row.update(predicted_ploidy='3x',uncertainty_reasons='NONE',selected_sites=30000)
            path=config.output_dir/'run_manifest.json'
            p.atomic_json(path,dict(result=row))
            return dict(manifest=path)
        config=p.CallConfig(self.bam,self.assets,self.root/'out')
        with patch.object(p,'build_bam_counts',side_effect=counts),patch.object(p,'run_ploidy_model',side_effect=model):
            with self.assertRaisesRegex(RuntimeError,'interrupted'):
                p.call(config)
            self.assertEqual(counts_calls,['stage_1','stage_2'])
            result=p.call(config)
        self.assertEqual(counts_calls,['stage_1','stage_2','stage_2'])
        self.assertEqual(self.result(result)['predicted_ploidy'],'3x')
        self.assertEqual(self.result(result)['status'],'CALLED')
        self.assertIn('LOW_BOOTSTRAP_SUPPORT',p.read_json(result.parent/'details/run.json')['result']['warnings'])
        self.assertEqual(p.read_json(result.parent/'details/run.json')['result']['reason'],'NONE')
        self.assertEqual((self.root/'out/regions.bed').read_text(),'chr1\t50\t1000\n')
        self.assertFalse((self.root/'out/details/stage_3').exists())


if __name__=='__main__':
    unittest.main()
