"""Real GenMap -> samtools -> default model smoke validation; no mocked tools."""
from pathlib import Path
import csv
import json
import os
import subprocess
import sys

import numpy as np
import pysam

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'validation' / 'synthetic_e2e'


def run(*args):
    subprocess.run([sys.executable, str(ROOT/'mapploid'), *map(str,args)], check=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    reference=OUT/'reference.fa'
    rng=np.random.default_rng(482)
    sequences={f'chr{i}':''.join(rng.choice(list('ACGT'),1500)) for i in range(1,7)}
    if not reference.exists():
        reference.write_text(''.join(f'>{name}\n{seq}\n' for name,seq in sequences.items()))
    assets=OUT/'assets'
    run('prepare','--reference',reference,'--read-length',150,'--threads',2,'--out',assets)
    results=[]
    for ploidy in (2,3,4):
        bam=OUT/f'synthetic_{ploidy}x.bam'
        if not bam.exists():
            header=dict(HD=dict(VN='1.6',SO='coordinate'),SQ=[dict(SN=k,LN=len(v)) for k,v in sequences.items()])
            means={2:[.5],3:[1/3,2/3],4:[.25,.5,.75]}[ploidy]
            weights={2:[1],3:[.5,.5],4:[.35,.3,.35]}[ploidy]
            draws=np.random.default_rng(120+ploidy)
            with pysam.AlignmentFile(str(bam),'wb',header=header) as handle:
                for tid,sequence in enumerate(sequences.values()):
                    for start in range(0,len(sequence),150):
                        seq=sequence[start:start+150]
                        values=[]
                        for base in seq:
                            mean=draws.choice(means,p=weights)
                            concentration=(1-.01)/.01
                            probability=draws.beta(mean*concentration,(1-mean)*concentration)
                            alternative=int(draws.binomial(30,probability))
                            alt=next(x for x in 'ACGT' if x!=base)
                            values.append([alt]*alternative+[base]*(30-alternative))
                        for read_id in range(30):
                            read=pysam.AlignedSegment()
                            read.query_name=f'{tid}_{start}_{read_id}'
                            read.query_sequence=''.join(v[read_id] for v in values)
                            read.flag=0;read.reference_id=tid;read.reference_start=start
                            read.mapping_quality=60;read.cigar=[(0,150)]
                            read.query_qualities=pysam.qualitystring_to_array('I'*150)
                            handle.write(read)
            pysam.index(str(bam))
        destination=OUT/f'result_{ploidy}x'
        run('call','--bam',bam,'--assets',assets,'--out',destination)
        with (destination/'result.tsv').open() as handle:
            row=next(csv.DictReader(handle,delimiter='\t'))
        expected={2:'2x',3:'3x',4:'4x'}[ploidy]
        assert row['predicted_ploidy']==expected,row
        results.append(dict(truth=ploidy,**row))
    (OUT/'verification.json').write_text(json.dumps(dict(passed=True,real_genmap=True,
        real_samtools=True,default_model_parameters=True,results=results),indent=2)+'\n')
    print('SYNTHETIC_END_TO_END_OK: 2x, 3x, 4x',flush=True)


if __name__=='__main__':
    main()
