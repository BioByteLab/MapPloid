from genmap_fixture import mock_genmap
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import csv
import json
import shutil
import unittest

import pysam
from mapploid import __version__
from mapploid.cli import main
import mapploid.batch as batch
from mapploid.pipeline import PrepareConfig, prepare, load_assets


class BamListTests(unittest.TestCase):
    def test_relative_paths_spaces_comments_and_bom(self):
        with TemporaryDirectory() as temporary:
            root=Path(temporary);listing=root/'bams.txt'
            listing.write_text('\ufeff# samples\n\nbams space/alpha.bam\n'+str(root/'beta.bam')+'\n')
            rows=batch.read_bam_list(listing)
            self.assertEqual([x['sample'] for x in rows],['alpha','beta'])
            self.assertEqual(rows[0]['bam'],root/'bams space/alpha.bam')

    def test_ambiguous_lists_fail_before_output_creation(self):
        with TemporaryDirectory() as temporary:
            root=Path(temporary);listing=root/'bams.txt';output=root/'out'
            for text,pattern in [('# empty\n','contains no paths'),
                                 ('a.bam\n'+str(root/'a.bam')+'\n','duplicate BAM path'),
                                 ('one/a.bam\ntwo/a.bam\n','duplicate sample name'),
                                 ('a.fastq\n','BAM file path'),
                                 ('bad name.bam\n','invalid sample name')]:
                listing.write_text(text)
                with self.subTest(text=text),self.assertRaisesRegex(ValueError,pattern):
                    batch.call_batch(listing,root/'missing-assets',output)
                self.assertFalse(output.exists())

    def test_cli_rejects_two_input_types_and_sample_name_for_batch(self):
        base=['call','--assets','assets','--out','out']
        for flags in [[],['--bam','a.bam','--bam-list','bams.txt'],['--bam-list','bams.txt','--sample','a']]:
            with self.subTest(flags=flags),redirect_stderr(StringIO()),self.assertRaises(SystemExit) as caught:
                main(base+flags)
            self.assertEqual(caught.exception.code,2)
        self.assertEqual(__version__,'0.10.5')

    def test_invalid_batch_filters_fail_before_opening_list(self):
        with TemporaryDirectory() as temporary:
            root=Path(temporary);output=root/'out'
            with self.assertRaisesRegex(ValueError,'allele'):
                batch.call_batch(root/'missing.txt',root/'missing-assets',output,allele_min=.5)
            self.assertFalse(output.exists())


@unittest.skipUnless(shutil.which('samtools'),'samtools required')
class BatchIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary=TemporaryDirectory();self.root=Path(self.temporary.name)
        ref=self.root/'ref.fa';ref.write_text('>chr1\n'+'A'*2000+'\n')
        bg=self.root/'scores.bg';bg.write_text('chr1\t0\t1951\t1\n')
        self.enterContext(mock_genmap(bg))
        self.assets=self.root/'assets'
        prepare(PrepareConfig(ref,self.assets,50,1,min_region_length=300))
        self.listing=self.root/'bams.txt';self.output=self.root/'batch'

    def tearDown(self):self.temporary.cleanup()

    def bam(self,name):
        path=self.root/f'{name}.bam'
        with pysam.AlignmentFile(str(path),'wb',header={
            'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':2000}]}):pass
        pysam.index(str(path));return path

    def rows(self):
        with (self.output/'result.tsv').open() as handle:return list(csv.DictReader(handle,delimiter='\t'))

    def command(self):
        return ['call','--bam-list',str(self.listing),'--assets',str(self.assets),'--out',str(self.output)]

    def test_batch_cli_continues_after_failure_and_retries_without_recomputing_success(self):
        self.bam('alpha');self.bam('beta')
        self.listing.write_text('alpha.bam\nmissing.bam\nbeta.bam\n')
        with redirect_stdout(StringIO()),redirect_stderr(StringIO()):self.assertEqual(main(self.command()),1)
        rows=self.rows();self.assertEqual([x['sample'] for x in rows],['alpha','missing','beta'])
        self.assertEqual([x['status'] for x in rows],['UNCERTAIN','FAILED','UNCERTAIN'])
        self.assertEqual(json.loads((self.output/'samples/alpha/details/run.json').read_text())['result']['reason'],'NO_COVERAGE')
        self.assertTrue(json.loads((self.output/'samples/missing/details/run.json').read_text())['result']['reason'])
        for row in (rows[0],rows[2]):self.assertTrue((self.output/row['result_file']).is_file())
        alpha=self.output/'samples/alpha/result.tsv';stamp=alpha.stat().st_mtime_ns
        self.bam('missing')
        with redirect_stdout(StringIO()),redirect_stderr(StringIO()):self.assertEqual(main(self.command()),0)
        self.assertEqual(alpha.stat().st_mtime_ns,stamp)
        self.assertEqual([x['status'] for x in self.rows()],['UNCERTAIN']*3)
        self.assertEqual(json.loads((self.output/'details/run.json').read_text())['state'],'COMPLETE')
        self.assertEqual(json.loads((self.output/'samples/alpha/details/run.json').read_text())['request']['model_site_limits'],dict(min_sites=25000,max_sites=100000))
        with (self.output/'samples/alpha/result.tsv').open() as handle:single=next(csv.DictReader(handle,delimiter='\t'))
        for key,value in single.items():self.assertEqual(self.rows()[0][key],value)

    def test_interruption_preserves_completed_sample_and_resumes(self):
        self.bam('alpha');self.bam('beta');self.listing.write_text('alpha.bam\nbeta.bam\n')
        real_call=batch.call
        def interrupted(config):
            if config.sample=='beta':raise KeyboardInterrupt()
            return real_call(config)
        with patch.object(batch,'call',side_effect=interrupted),redirect_stdout(StringIO()),redirect_stderr(StringIO()):
            self.assertEqual(main(self.command()),130)
        self.assertEqual(len(self.rows()),1)
        self.assertEqual(json.loads((self.output/'details/run.json').read_text())['state'],'INTERRUPTED')
        alpha=self.output/'samples/alpha/result.tsv';stamp=alpha.stat().st_mtime_ns
        with redirect_stdout(StringIO()),redirect_stderr(StringIO()):self.assertEqual(main(self.command()),0)
        self.assertEqual(alpha.stat().st_mtime_ns,stamp)
        self.assertEqual(len(self.rows()),2)

    def test_old_reference_assets_work_but_changed_batch_settings_are_rejected(self):
        manifest=self.assets/'assets.json';data=json.loads(manifest.read_text());data['version']='0.10.1';manifest.write_text(json.dumps(data))
        self.assertEqual(load_assets(self.assets)['version'],'0.10.1')
        self.bam('alpha');self.listing.write_text('alpha.bam\n')
        with redirect_stdout(StringIO()),redirect_stderr(StringIO()):self.assertEqual(main(self.command()),0)
        before=(self.output/'result.tsv').read_bytes()
        with redirect_stdout(StringIO()),redirect_stderr(StringIO()):
            self.assertEqual(main(self.command()+['--min-mapq','30']),1)
        self.assertEqual((self.output/'result.tsv').read_bytes(),before)
        self.listing.write_text('alpha.bam\n# edited\n')
        with redirect_stdout(StringIO()),redirect_stderr(StringIO()):self.assertEqual(main(self.command()),1)
        self.assertEqual((self.output/'result.tsv').read_bytes(),before)

    def test_real_covered_bam_and_empty_bam_produce_called_and_uncertain_rows(self):
        import numpy as np
        ref=self.root/'three.fa';ref.write_text(''.join(f'>chr{i}\n'+ 'A'*12000+'\n' for i in range(1,4)))
        scores=self.root/'three.bg';scores.write_text(''.join(f'chr{i}\t0\t11951\t1\n' for i in range(1,4)))
        self.enterContext(mock_genmap(scores))
        assets=self.root/'three-assets'
        prepare(PrepareConfig(ref,assets,50,1))
        header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':f'chr{i}','LN':12000} for i in range(1,4)]}
        rng=np.random.default_rng(20261006)
        covered=self.root/'triploid.bam'
        with pysam.AlignmentFile(str(covered),'wb',header=header) as handle:
            for chrom in range(3):
                for start in range(0,12000,1000):
                    for j in range(30):
                        read=pysam.AlignedSegment();read.query_name=f'{chrom}_{start}_{j}'
                        read.query_sequence=''.join(np.where(rng.random(1000)<1/3,'C','A'))
                        read.flag=0;read.reference_id=chrom;read.reference_start=start
                        read.mapping_quality=60;read.cigar=[(0,1000)]
                        read.query_qualities=pysam.qualitystring_to_array('I'*1000);handle.write(read)
        pysam.index(str(covered))
        empty=self.root/'empty.bam'
        with pysam.AlignmentFile(str(empty),'wb',header=header):pass
        pysam.index(str(empty))
        self.listing.write_text('triploid.bam\nempty.bam\n')
        command=['call','--bam-list',str(self.listing),'--assets',str(assets),'--out',str(self.output)]
        with redirect_stdout(StringIO()),redirect_stderr(StringIO()):self.assertEqual(main(command),0)
        rows=self.rows()
        self.assertEqual([x['status'] for x in rows],['CALLED','UNCERTAIN'])
        self.assertEqual(rows[0]['predicted_ploidy'],'3x')
        self.assertGreaterEqual(int(rows[0]['selected_sites']),25000)
        self.assertEqual(json.loads((self.output/'samples/empty/details/run.json').read_text())['result']['reason'],'NO_COVERAGE')
