"""Reference preparation and single/BAM-list ploidy calls."""
import argparse
import os
from pathlib import Path
import sys

from . import __version__


def add_help(parser):
    parser.add_argument('-h','--help',action='help',help='show this help message and exit')


def build_parser():
    parser = argparse.ArgumentParser(
        prog='mapploid', description='Infer sample ploidy: prepare reference regions, then analyze one or more BAM samples.',
        formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False,
        epilog='Prepare the reference first:\n  mapploid prepare --reference genome.fa --out reference_assets\n'
               'Analyze one sample:\n  mapploid call --bam sample.bam --assets reference_assets --out result\n'
               'Analyze a batch:\n  mapploid call --bam-list bams.txt --assets reference_assets --out batch_results')
    add_help(parser)
    parser.add_argument('--version',action='version',version=__version__,help='show the software version and exit')
    sub = parser.add_subparsers(dest='command',required=True,title='commands')
    prepare = sub.add_parser('prepare',help='prepare reusable analysis regions from a reference FASTA',add_help=False,
                             description='Prepare once per reference and region settings. Use the resulting --out directory as call --assets.')
    add_help(prepare)
    prepare.add_argument('--reference',type=Path,required=True,metavar='FASTA',help='reference genome as an uncompressed FASTA file')
    prepare.add_argument('--kmer',type=int,default=150,metavar='BP',help='GenMap k-mer length in bp (default: 150)')
    prepare.add_argument('--min-region-length',type=int,default=10_000,metavar='BP',help='minimum continuous high-mappability region length in bp (default: 10000)')
    prepare.add_argument('--out',type=Path,required=True,metavar='DIR',help='directory to create, e.g. reference_assets; pass this directory to call --assets')
    prepare.add_argument('--threads',type=int,default=8,metavar='N',help='threads for reference preparation (default: 8)')
    prepare.add_argument('--exclude-bed',type=Path,action='append',default=[],metavar='BED',help='BED file of regions to exclude, e.g. organellar sequences; may be repeated')
    call = sub.add_parser('call',help='infer ploidy for one BAM or a list of BAMs',add_help=False,
                          formatter_class=argparse.RawDescriptionHelpFormatter,
                          description='Infer sample ploidy using the reference directory created by prepare.',
                          epilog='Examples:\n  mapploid call --bam sample.bam --assets reference_assets --out result\n'
                                 '  mapploid call --bam-list bams.txt --assets reference_assets --out batch_results\n\n'
                                 'bams.txt: one BAM path per line; relative paths are resolved from the list directory.\n'
                                 'Batch outputs: batch_results/result.tsv summarizes all samples;\n'
                                 '               batch_results/samples/SAMPLE/ contains each sample\'s results and logs.')
    add_help(call)
    inputs = call.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--bam',type=Path,metavar='FILE',help='one sample BAM file, coordinate-sorted with a BAI/CSI index')
    inputs.add_argument('--bam-list',type=Path,metavar='FILE',help='text file with one BAM path per line; each BAM must be sorted and indexed')
    call.add_argument('--assets',type=Path,required=True,metavar='DIR',help='reference directory created by prepare --out, e.g. reference_assets; shared by samples using the same reference')
    call.add_argument('--out',type=Path,required=True,metavar='DIR',help='output directory: result.tsv for one sample, or a summary and sample subdirectories for a batch')
    call.add_argument('--sample',metavar='NAME',help='name for --bam only (default: BAM filename without .bam); batches use each BAM filename')
    call.add_argument('--keep-details', action='store_true', help='retain stage tables and intermediate site data after completion')
    call.add_argument('--export-figures', action='store_true', help='also export the report evidence figure as PDF and 300-dpi PNG')
    filters = call.add_argument_group('site filters (optional)')
    filters.add_argument('--allele-min',type=float,default=.15,metavar='FLOAT',help='minimum frequency of each of the two leading alleles (default: 0.15)')
    filters.add_argument('--allele-max',type=float,default=.85,metavar='FLOAT',help='maximum frequency of each of the two leading alleles (default: 0.85)')
    filters.add_argument('--min-depth',type=int,metavar='N',help='minimum site depth after quality filtering (default: automatic)')
    filters.add_argument('--max-depth',type=int,metavar='N',help='maximum site depth after quality filtering (default: automatic)')
    filters.add_argument('--min-mapq',type=int,default=20,metavar='N',help='minimum read mapping quality, 0-255 (default: 20)')
    filters.add_argument('--min-baseq',type=int,default=20,metavar='N',help='minimum base quality, 0-255 (default: 20)')
    filters.add_argument('--min-allele-reads',type=int,default=2,metavar='N',help='minimum supporting reads for each allele (default: 2)')
    filters.add_argument('--pair-purity',type=float,default=.95,metavar='FLOAT',help='minimum fraction of valid base counts explained by the two leading alleles (default: 0.95)')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command=='call' and args.bam_list is not None and args.sample is not None:
        parser.error('--sample is only available with --bam; batches use each BAM filename as the sample name')
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
        os.environ.setdefault(key,'1')
    from .pipeline import PrepareConfig, CallConfig, prepare, call
    try:
        if args.command=='prepare':
            path = prepare(PrepareConfig(
                reference=args.reference,output_dir=args.out,kmer_length=args.kmer,
                threads=args.threads,exclude_beds=args.exclude_bed,
                min_region_length=args.min_region_length))
            path = path.parent
        else:
            filters = {name:getattr(args,name) for name in (
                'allele_min','allele_max','min_depth','max_depth','min_mapq','min_baseq','min_allele_reads','pair_purity')}
            filters.update(keep_details=args.keep_details, export_figures=args.export_figures)
            if args.bam_list is not None:
                from .batch import call_batch
                result = call_batch(args.bam_list,args.assets,args.out,**filters)
                print(result.summary)
                return 1 if result.failed else 0
            path = call(CallConfig(args.bam,args.assets,args.out,args.sample,**filters))
        print(path)
        return 0
    except (OSError,ValueError,RuntimeError) as exc:
        print(f'[MapPloid] ERROR: {exc}',file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print('[MapPloid] Interrupted; rerun the same command to resume.',file=sys.stderr)
        return 130
