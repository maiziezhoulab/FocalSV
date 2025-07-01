from tqdm import tqdm

def load_contigs(fasta_path, logger):
    logger.info("loading "+fasta_path)
    with open(fasta_path,'r') as f:
        s = f.readlines()
    # dc = {}
    # for i in tqdm(range(0,len(s),2),desc='load contig'):
    #     name = s[i][1:-1]
    #     seq = s[i+1][:-1]
    #     dc[name]=seq
    dc = {}

    for line in tqdm(s,desc='load contig'):
        if '>' in line:
            cur_name = line[1:-1].split()[0]
        else:
            if cur_name in dc:
                dc[cur_name].append(line[:-1])
            else:
                dc[cur_name] = [line[:-1]]

    for name in tqdm(dc):
        dc[name] = ''.join(dc[name])

    logger.info("finish loading")
    return dc