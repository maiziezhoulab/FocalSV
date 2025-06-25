import os

def filter(vcf, out_dir,chrs=None,dipcall=False,delly = False):
    if chrs is None:
        chrs = list()
        for chrnum in range(1,23):
            chrs.append('chr'+str(chrnum)) 


                    
    out_vcf = out_dir+'/'+vcf.replace('.vcf','').split('/')[-1]+'_DEL_INS_noXY.vcf'
    out_INS_vcf = out_dir+'/'+vcf.replace('.vcf','').split('/')[-1]+'_INS_noXY.vcf'
    out_DEL_vcf = out_dir+'/'+vcf.replace('.vcf','').split('/')[-1]+'_DEL_noXY.vcf'
    #print(vcf)
    #print(out_vcf)
    #print(out_INS_vcf)
    #print(out_DEL_vcf)
    header = []
    body = []
    with open(vcf,'r') as fin:
        for line in fin:
            if line[0]=='#':
                header.append(line)
            else:
                body.append(line)    
    if delly:
        dellyheader = "/data/maiziezhou_lab/CanLuo/Link_Reads_Project_Complete_Genomics/bin/delly_header"
        with open(dellyheader,'r') as fd:
            dh = fd.readlines()
            header = header[:-1] + dh +header[-1:]


    if 1:
        with open(out_vcf,'w') as fout:
            with open(out_INS_vcf,'w') as fins:
                with open(out_DEL_vcf,'w') as fdel:
                    fout.writelines(header)
                    fins.writelines(header)
                    fdel.writelines(header)
                    for line in body:
                        if line[0] == '#':
                            if 'contig=' in line:
                                chrinfo = line.split(",")[0].split("=")[-1]
                                if (chrinfo in chrs) or 1:
                                    fout.write(line)
                                    fins.write(line)
                                    fdel.write(line)
                            else:
                                fout.write(line)
                                fins.write(line)
                                fdel.write(line)
                        else:
                            line = line.replace("SVLEN=>", "SVLEN=") #for NanoVar
                            svinfo = line.split("\t")
                            chrinfo = svinfo[0]
                            #svtype = svinfo[7].split(";")[0]
                            #if chrinfo in chrs and svtype == 'SVTYPE=DEL':
                            if dipcall:
                                if chrinfo in chrs and (len(svinfo[3])-len(svinfo[4]))>49:
                                    fout.write(line)
                                    fdel.write(line)
                                #elif chrinfo in chrs and svtype == 'SVTYPE=INS':
                                elif chrinfo in chrs and (len(svinfo[3])-len(svinfo[4]))<-49:
                                    fout.write(line)
                                    fins.write(line)
                            else:
                                if chrinfo in chrs and 'SVTYPE=DEL' in svinfo[7]:
                                    fout.write(line)
                                    fdel.write(line)
                                #elif chrinfo in chrs and svtype == 'SVTYPE=INS':
                                elif chrinfo in chrs and 'SVTYPE=INS' in svinfo[7]:
                                    fout.write(line)
                                    fins.write(line)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('--vcf','-v',)
    parser.add_argument('--out_dir','-o_dir')
    parser.add_argument('--chrs', nargs='+')
    parser.add_argument('--dipcall',action="store_true")
    parser.add_argument('--delly',action="store_true")


    args = parser.parse_args()

    filter(args.vcf, args.out_dir, args.chrs, args.dipcall,args.delly)
