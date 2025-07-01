from argparse import ArgumentParser
from subprocess import Popen
import os
import logging
from extract_contig_signature_CLR import extract_contig_sig_CLR
from extract_contig_signature_ONT import extract_contig_sig_ONT
from extract_contig_signature_CCS import extract_contig_sig_CCS
from extract_reads_signature import extract_reads_signature
from FP_filter_v1 import FP_filter
from remove_redundancy import remove_redundancy



def reformat_fasta(contig_file,outfile,hp):
	cnt = 0
	with open(outfile,'w') as fw:
		with open(contig_file,'r') as f:
			for line in f:
				if line[0]=='>':
					line = '>contig_%s_%d\n'%(hp,cnt)
					cnt+=1
				fw.write(line)
	return 
    
import subprocess

def extract_chrom_fasta(fasta_path, chrom, out_path):
    """
    Extract a chromosome FASTA from a genome FASTA using samtools faidx.

    Parameters:
        fasta_path (str): Path to the reference genome FASTA file.
        chrom (str): Chromosome name to extract (e.g., 'chr21').
        out_path (str): Path to write the extracted chromosome FASTA.
    """
    try:
        # Make sure the index exists
        fai_path = fasta_path + ".fai"
        subprocess.run(["samtools", "faidx", fasta_path], check=True)

        # Extract the chromosome
        with open(out_path, "w") as out_f:
            subprocess.run(["samtools", "faidx", fasta_path, chrom], check=True, stdout=out_f)

        print(f"Chromosome {chrom} written to {out_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error running samtools: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


def dippav_variant_call(data_type,
		read_bam_file,
		reference_path,
	hp1_contig_path,
	hp2_contig_path,
	output_dir,
	chr_num,
	header_file=None,
	n_thread=10,
	mem_per_thread = '1G'):

	assert data_type in ['CCS','CLR','ONT']

	# data_type = 'CLR'

	## set logger
	logging.basicConfig(
	format='%(asctime)s %(levelname)-8s %(message)s',
		level=logging.INFO,
		datefmt='%Y-%m-%d %H:%M:%S')
	global logger
	logger = logging.getLogger(" ")
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	code_dir=os.path.dirname(os.path.realpath(__file__))+'/'
	if not header_file:
		header_file = code_dir + 'header'

	outhp1 = output_dir+'/hp1.fa'
	outhp2 = output_dir+'/hp2.fa'
	contig_path = output_dir+'/assemblies.fa'
	prefix = contig_path.split('/')[-1].split('.')[0]
	bam_path = output_dir+'/'+prefix+'.sorted.bam'
	reference_path_chr = output_dir+"/ref_chr%d.fa"%chr_num


	#extract one chrom ref
	extract_chrom_fasta(reference_path, 'chr%d'%chr_num, output_dir+"/ref_chr%d.fa"%chr_num)
	
	

	reformat_fasta(hp1_contig_path ,outhp1,'hp1')
	reformat_fasta(hp2_contig_path ,outhp2,'hp2')
	
	cmd = "cat %s %s > %s"%(outhp1,outhp2,contig_path)
	Popen(cmd, shell = True).wait()


	logger.info("align contig to reference...")
	
	cmd = "minimap2 -a -x asm5 --cs -r2k -t %d \
		%s \
		%s \
			| samtools sort -@ %d -m %s > %s"%(n_thread,reference_path_chr,contig_path,n_thread, mem_per_thread,bam_path )
	logger.info(cmd)
	Popen(cmd,shell=True).wait()

	cmd = "samtools index "+bam_path
	logger.info(cmd)
	Popen(cmd,shell=True).wait()



	logger.info("Raw variant call by chromosome...")
	if data_type == 'CLR':
		extract_contig_sig_CLR(chr_number=chr_num,
					bam_path=bam_path,
					header_path=header_file,
					ref_path=reference_path_chr,
					contig_path=contig_path,
					output_dir=output_dir)
	elif data_type == 'ONT':
		extract_contig_sig_ONT(chr_number=chr_num,
					bam_path=bam_path,
					header_path=header_file,
					ref_path=reference_path_chr,
					contig_path=contig_path,
					output_dir=output_dir)
	else:
		extract_contig_sig_CCS(chr_number=chr_num,
					bam_path=bam_path,
					header_path=header_file,
					ref_path=reference_path_chr,
					contig_path=contig_path,
					output_dir=output_dir)


	logger.info("extract signatures from reads bam file...")
	extract_reads_signature(chr_number=chr_num,
				input_path= read_bam_file,
				output_dir=output_dir)


	### combine all vcf together

	cmd = "cat %s/dippav_variant_chr*vcf | grep -v \"#\" > %s/body "%(output_dir,output_dir)
	Popen(cmd,shell=True).wait()
	cmd = "cat %s %s/body > %s/dippav_raw_variant.vcf; rm %s/body "%(header_file,output_dir,output_dir,output_dir)
	Popen(cmd,shell=True).wait()



	logger.info("Filter out false positive...")
	vcf_path = output_dir+"/dippav_raw_variant.vcf"
	vcf_path_filtered = output_dir+"/dippav_variant_filtered.vcf"
	signature_dir = output_dir+'/reads_signature/'

	FP_filter(input_path=vcf_path,
		signature_dir=signature_dir,
		output_path=vcf_path_filtered
		)


	logger.info("Remove redundancy...")
	# vcf_path_noredun = output_dir+"/dippav_variant_chr%d_no_redundancy.vcf"%chr_num
	final_vcf_dir = output_dir+'/final_vcf/'

	remove_redundancy(vcf_path=vcf_path_filtered,
			output_dir=final_vcf_dir)













if __name__ == "__main__":
	parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
	parser.add_argument('--read_bam_file','-rbam')
	parser.add_argument('--hp1_contig_path','-hp1')
	parser.add_argument('--hp2_contig_path','-hp2')
	parser.add_argument('--reference_path','-ref', help = "single chrom ref")
	# parser.add_argument('--signature_dir','-sigd')
	parser.add_argument('--output_dir','-o')
	parser.add_argument('--data_type','-dtype',choices= ['Hifi','CLR','ONT'], help="CCS;CLR;ONT")
	parser.add_argument('--chr_num','-chr',type = int)
	###optional
	parser.add_argument('--header_file','-header',help='optional;if not set, will use the default header')
	parser.add_argument('--n_thread','-t', type = int, default = 10,help = 'default = 10')
	parser.add_argument('--mem_per_thread','-mempt', default = '1G',help = "Set maximum memory per thread; suffix K/M/G recognized; default = 1G")



	args = parser.parse_args()
	read_bam_file = args.read_bam_file
	hp1_contig_path = args.hp1_contig_path
	hp2_contig_path = args.hp2_contig_path
	# contig_path = args.contig_path
	reference_path = args.reference_path
	# signature_dir = args.signature_dir
	output_dir = args.output_dir
	data_type = args.data_type
	chr_num = args.chr_num
	###optional
	header_file = args.header_file
	n_thread = args.n_thread
	mem_per_thread = args.mem_per_thread

	if data_type == 'Hifi':
		data_type = 'CCS'

	# print(chr_num)
	dippav_variant_call(
		data_type,
		read_bam_file,
		reference_path,
	hp1_contig_path,
	hp2_contig_path,
	output_dir,
	chr_num,
	header_file,
	n_thread,
	mem_per_thread )




















