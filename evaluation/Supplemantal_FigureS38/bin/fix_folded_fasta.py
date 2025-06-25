import re
from collections import defaultdict

def fix_fasta_streaming(input_path, output_path, fai, hap_tag=None,):
    fai_dict = dict()
    with open(fai, "r") as f:
        for line in f:
            if line[0] != "#":
                line = line.rstrip("\n").split("\t")
                fai_dict[line[0]] = [int(i) for i in line[1:]]
    
    with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
        in_header = False
        header_parts = []
        sequence_lines = []

        contig_name_dict = defaultdict(int)

        for line in fin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if in_header:
                    # Write previous record
                    full_header = ''.join(header_parts)

                    m = re.search(r'Region_chr(\w+)_S(-?\d+)_E(-?\d+)', full_header)
                    if m:
                        chr_, start, end = m.groups()
                        block_name = full_header.split("/")[-1].split(".")[0].replace("_hp2", "").replace("_hp1", "")

                        start = int(start)
                        if start < 0:
                            start = 0
                        end = int(end)
                        if end > fai_dict["chr"+chr_][0]:
                            end = fai_dict["chr"+chr_][0]

                        if hap_tag is not None:
                            contig_name = f"{hap_tag}_chr{chr_}_{start}_{end}_{block_name}"
                        else:
                            contig_name = f"chr{chr_}_{start}_{end}_{block_name}"

                        contig_name_dict[contig_name] += 1

                        fout.write(">"+contig_name+"-"+str(contig_name_dict[contig_name])+"\n")

                        # for seq in sequence_lines:
                        #     fout.write(seq + '\n')
                        fout.write("".join(sequence_lines)+"\n")
                    else:
                        raise ValueError(f"Cannot parse header: {full_header}")
                    
                # Start new header
                header_parts = [line[1:]]  # Remove '>'
                sequence_lines = []
                in_header = True
            elif in_header and not re.match(r'^[ACGTNacgtn]+$', line):
                # Continuation of folded header
                header_parts.append(line)
            else:
                # Sequence line
                sequence_lines.append(line)

        # Write last record
        if in_header:
            full_header = ''.join(header_parts)
            m = re.search(r'Region_chr(\w+)_S(-?\d+)_E(-?\d+)', full_header)
            if m:
                chr_, start, end = m.groups()

                start = int(start)
                if start < 0:
                    start = 0
                end = int(end)
                if end > fai_dict["chr"+chr_][0]:
                    end = fai_dict["chr"+chr_][0]
                    
                if hap_tag is not None:
                    fout.write(f">{hap_tag}_chr{chr_}_{start}_{end}\n")
                else:
                    fout.write(f">chr{chr_}_{start}_{end}\n")

                # for seq in sequence_lines:
                #     fout.write(seq + '\n')
                fout.write("".join(sequence_lines)+"\n")
            else:
                raise ValueError(f"Cannot parse header at EOF: {full_header}")

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('--input_fasta','-i',)
    parser.add_argument('--output_fasta', '-o')
    parser.add_argument('--hap_tag')
    parser.add_argument('--fai')

    args = parser.parse_args()

    fix_fasta_streaming(args.input_fasta, args.output_fasta, args.fai, args.hap_tag)
