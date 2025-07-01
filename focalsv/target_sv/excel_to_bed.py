import pandas as pd 
from collections import defaultdict
df = pd.read_excel("High_confidence_callset.xlsx")
dc = defaultdict(list)
use_chrs = ['chr'+str(i) for i in range(1,23)]
for i in range(df.shape[0]):
    chrom1, pos1, chrom2, pos2, svsize, svtype = df.iloc[i,1:7]
    if chrom1 in use_chrs and chrom2 in use_chrs:
        if svtype in ['INV','DUP']:
            if svtype == 'DUP':
                if svsize >= 5000000:
                    continue
            dc[svtype].append((chrom1, min(pos1,pos2), max(pos1, pos2)))
        elif svtype == 'TRA':
            if chrom1!=chrom2:
                dc[svtype].append((chrom1,pos1,chrom2,pos2))
flank = 50000
cnt = 0
sv_cnt = 0
with open("HCC1395_SV_rich_regions.bed",'w') as f:
    for svtype, svs in dc.items():
        print(svtype, len(svs))
        sv_cnt+= len(svs)
        if svtype == 'DUP':
            
            for sv in svs:
                f.write(f"{sv[0]}\t{int(max(0,sv[1]-flank))}\t{int(sv[2]+flank)}\t{svtype}\n")
        elif svtype == 'INV':
            for sv in svs:
                f.write(f"{sv[0]}\t{int(max(0,sv[1]-flank))}\t{int(sv[1]+flank)}\t{int(max(0,sv[2]-flank))}\t{int(sv[2]+flank)}\t{svtype}\n")
        else:
            for sv in svs:
                cnt+=1
                chrom1,pos1,chrom2,pos2 = sv 
                f.write(f"{chrom1}\t{int(max(0,pos1-flank))}\t{int(pos1+flank)}\t{chrom2}\t{int(max(0,pos2-flank))}\t{int(pos2+flank)}\t{svtype}\n")
                # f.write(f"{chrom2}\t{int(max(0,pos2-flank))}\t{int(pos2+flank)}\t{svtype}\tTRA.{cnt}\n")


print("total svs: ", sv_cnt)