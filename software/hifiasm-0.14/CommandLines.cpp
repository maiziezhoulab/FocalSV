#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <sys/time.h>
#include "CommandLines.h"
#include "ketopt.h"

#define DEFAULT_OUTPUT "hifiasm.asm"

hifiasm_opt_t asm_opt;

static ko_longopt_t long_options[] = {
	{ "version",       ko_no_argument, 300 },
	{ "dbg-gfa",       ko_no_argument, 301 },
	{ "write-paf",     ko_no_argument, 302 },
	{ "write-ec",      ko_no_argument, 303 },
	{ "skip-triobin",  ko_no_argument, 304 },
	{ "max-od-ec",     ko_no_argument, 305 },
	{ "max-od-final",  ko_no_argument, 306 },
	{ "ex-list",       ko_required_argument, 307 },
	{ "ex-iter",       ko_required_argument, 308 },
    { "purge-cov",     ko_required_argument, 309 },
    { "pri-range",     ko_required_argument, 310 },
    { "high-het",      ko_no_argument, 311 },
    { "lowQ",          ko_required_argument, 312 },
	{ "min-hist-cnt",  ko_required_argument, 313 },
    { "h1",      ko_required_argument, 314 },
    { "h2",      ko_required_argument, 315 },
    { "enzyme",      ko_required_argument, 316 },
    { "b-cov",      ko_required_argument, 317 },
    { "h-cov",      ko_required_argument, 318 },
    { "m-rate",      ko_required_argument, 319 },
	{ 0, 0, 0 }
};

double Get_T(void)
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec+t.tv_usec/1000000.0;
}

void Print_H(hifiasm_opt_t* asm_opt)
{
    fprintf(stderr, "Usage: hifiasm [options] <in_1.fq> <in_2.fq> <...>\n");
    fprintf(stderr, "Options:\n");
	fprintf(stderr, "  Input/Output:\n");
    fprintf(stderr, "    -o STR      prefix of output files [%s]\n", asm_opt->output_file_name);
    fprintf(stderr, "    -i          ignore saved read correction and overlaps\n");
    fprintf(stderr, "    -t INT      number of threads [%d]\n", asm_opt->thread_num);
    fprintf(stderr, "    -z INT      length of adapters that should be removed [%d]\n", asm_opt->adapterLen);
    fprintf(stderr, "    --version   show version number\n");
	fprintf(stderr, "  Overlap/Error correction:\n");
    fprintf(stderr, "    -k INT      k-mer length (must be <64) [%d]\n", asm_opt->k_mer_length);
	fprintf(stderr, "    -w INT      minimizer window size [%d]\n", asm_opt->mz_win);
	fprintf(stderr, "    -f INT      number of bits for bloom filter; 0 to disable [%d]\n", asm_opt->bf_shift);
	fprintf(stderr, "    -D FLOAT    drop k-mers occurring >FLOAT*coverage times [%.1f]\n", asm_opt->high_factor);
	fprintf(stderr, "    -N INT      consider up to max(-D*coverage,-N) overlaps for each oriented read [%d]\n", asm_opt->max_n_chain);
    fprintf(stderr, "    -r INT      round of correction [%d]\n", asm_opt->number_of_round);
    fprintf(stderr, "  Assembly:\n");
    fprintf(stderr, "    -a INT      round of assembly cleaning [%d]\n", asm_opt->clean_round);
    fprintf(stderr, "    -m INT      pop bubbles of <INT in size in contig graphs [%lld]\n", asm_opt->large_pop_bubble_size);
    fprintf(stderr, "    -p INT      pop bubbles of <INT in size in unitig graphs [%lld]\n", asm_opt->small_pop_bubble_size);
    fprintf(stderr, "    -n INT      remove tip unitigs composed of <=INT reads [%d]\n", asm_opt->max_short_tip);
    fprintf(stderr, "    -x FLOAT    max overlap drop ratio [%.2g]\n", asm_opt->max_drop_rate);
    fprintf(stderr, "    -y FLOAT    min overlap drop ratio [%.2g]\n", asm_opt->min_drop_rate);
    fprintf(stderr, "    -u          disable post join contigs step which may improve N50\n");
    fprintf(stderr, "    --lowQ      INT\n");
    fprintf(stderr, "                output contig regions with >=INT%% inconsistency in BED format; 0 to disable [%d]\n", asm_opt->bed_inconsist_rate);
    fprintf(stderr, "    --b-cov     INT\n");
    fprintf(stderr, "                break contigs at positions with <INT-fold coverage; work with '--m-rate'; 0 to disable [%d]\n", asm_opt->b_low_cov);
    fprintf(stderr, "    --h-cov     INT\n");
    fprintf(stderr, "                break contigs at positions with >INT-fold coverage; work with '--m-rate'; -1 to disable [%d]\n", asm_opt->b_high_cov);
    fprintf(stderr, "    --m-rate    FLOAT\n");
    fprintf(stderr, "                break contigs at positions with <=FLOAT*coverage exact overlaps;\n");
    fprintf(stderr, "                only work with '--b-cov' or '--h-cov'[%.2f]\n", asm_opt->m_rate);
    
//	fprintf(stderr, "    --pri-range INT1[,INT2]\n");
//	fprintf(stderr, "                keep contigs with coverage in this range in p_ctg.gfa; -1 to disable [auto,inf]\n");

    fprintf(stderr, "  Trio-partition:\n");
    fprintf(stderr, "    -1 FILE     hap1/paternal k-mer dump generated by \"yak count\" []\n");
    fprintf(stderr, "    -2 FILE     hap2/maternal k-mer dump generated by \"yak count\" []\n");
    fprintf(stderr, "    -c INT      lower bound of the binned k-mer's frequency [%d]\n", asm_opt->min_cnt);
    fprintf(stderr, "    -d INT      upper bound of the binned k-mer's frequency [%d]\n", asm_opt->mid_cnt);
	fprintf(stderr, "    -3 FILE     list of hap1/paternal read names []\n");
	fprintf(stderr, "    -4 FILE     list of hap2/maternal read names []\n");

    fprintf(stderr, "  Purge-dups:\n");
    fprintf(stderr, "    -l INT      purge level. 0: no purging; 1: light; 2: aggressive [0 for trio; 2 for unzip]\n");
    fprintf(stderr, "    -s FLOAT    similarity threshold for duplicate haplotigs [%g]\n", 
                                     asm_opt->purge_simi_rate);
    fprintf(stderr, "    -O INT      min number of overlapped reads for duplicate haplotigs [%d]\n", 
                                     asm_opt->purge_overlap_len);
    fprintf(stderr, "    --purge-cov INT\n");
    fprintf(stderr, "                coverage upper bound of Purge-dups [auto]\n");
    fprintf(stderr, "    --high-het  enable this mode for high heterozygosity sample [experimental, not stable]\n");
    
    fprintf(stderr, "  Hi-C-partition [experimental, not stable]:\n");
    fprintf(stderr, "    --h1 FILEs   file names of Hi-C R1  [r1_1.fq,r1_2.fq,...]\n");
    fprintf(stderr, "    --h2 FILEs   file names of Hi-C R2  [r2_1.fq,r2_2.fq,...]\n");

    fprintf(stderr, "Example: ./hifiasm -o NA12878.asm -t 32 NA12878.fq.gz\n");
    fprintf(stderr, "See `man ./hifiasm.1' for detailed description of these command-line options.\n");
}

void init_opt(hifiasm_opt_t* asm_opt)
{
	memset(asm_opt, 0, sizeof(hifiasm_opt_t));
	asm_opt->flag = 0;
    asm_opt->coverage = -1;
    asm_opt->num_reads = 0;
    asm_opt->read_file_names = NULL;
    asm_opt->output_file_name = (char*)(DEFAULT_OUTPUT);
    asm_opt->required_read_name = NULL;
    asm_opt->hic_enzymes = NULL;
    asm_opt->hic_reads[0] = NULL;
    asm_opt->hic_reads[1] = NULL;
    asm_opt->thread_num = 1;
    asm_opt->k_mer_length = 51;
    asm_opt->hic_mer_length = 31;
	asm_opt->mz_win = 51;
	asm_opt->bf_shift = 37;
	asm_opt->high_factor = 5.0;
	asm_opt->max_ov_diff_ec = 0.04;
	asm_opt->max_ov_diff_final = 0.03;
	asm_opt->hom_cov = 20;
    asm_opt->het_cov = -1024;
	asm_opt->max_n_chain = 100;
	asm_opt->min_hist_kmer_cnt = 5;
    asm_opt->load_index_from_disk = 1;
    asm_opt->write_index_to_disk = 1;
    asm_opt->number_of_round = 3;
    asm_opt->adapterLen = 0;
    asm_opt->clean_round = 4;
    asm_opt->small_pop_bubble_size = 100000;
    asm_opt->large_pop_bubble_size = 10000000;
    asm_opt->min_drop_rate = 0.2;
    asm_opt->max_drop_rate = 0.8;
    asm_opt->max_hang_Len = 1000;
    asm_opt->max_hang_rate = 0.8;
    asm_opt->gap_fuzz = 1000;
    asm_opt->min_overlap_Len = 50;
    asm_opt->min_overlap_coverage = 0;
    asm_opt->max_short_tip = 3;
    asm_opt->min_cnt = 2;
    asm_opt->mid_cnt = 5;
    asm_opt->purge_level_primary = 2;
    asm_opt->purge_level_trio = 0;
    asm_opt->purge_simi_rate = 0.75;
    asm_opt->purge_simi_rate_hic = 0.85;
    asm_opt->purge_overlap_len = 1;
    asm_opt->purge_overlap_len_hic = 50;
    asm_opt->recover_atg_cov_min = -1024;
    asm_opt->recover_atg_cov_max = INT_MAX;
    asm_opt->hom_global_coverage = -1;
    asm_opt->bed_inconsist_rate = 70;
    asm_opt->hic_inconsist_rate = 30;
    ///asm_opt->bub_mer_length = 3;
    asm_opt->bub_mer_length = 1000000;
    asm_opt->b_low_cov = 0;
    asm_opt->b_high_cov = -1;
    asm_opt->m_rate = 0.75;
}

void destory_enzyme(enzyme* f)
{
    int i;
    if(f != NULL)
    {
        for (i = 0; i < f->n; i++)
        {
            free(f->a[i]);
        }
        free(f->a);
        free(f->l);
        free(f);
    }
}

void destory_opt(hifiasm_opt_t* asm_opt)
{
    if(asm_opt->read_file_names != NULL) free(asm_opt->read_file_names);
    if(asm_opt->hic_enzymes != NULL) destory_enzyme(asm_opt->hic_enzymes);
    if(asm_opt->hic_reads[0] != NULL) destory_enzyme(asm_opt->hic_reads[0]);
    if(asm_opt->hic_reads[1] != NULL) destory_enzyme(asm_opt->hic_reads[1]);
}

void ha_opt_reset_to_round(hifiasm_opt_t* asm_opt, int round)
{
    asm_opt->num_bases = 0;
    asm_opt->num_corrected_bases = 0;
    asm_opt->num_recorrected_bases = 0;
	asm_opt->mem_buf = 0;
    asm_opt->roundID = round;
}

void ha_opt_update_cov(hifiasm_opt_t *opt, int hom_cov)
{
	int max_n_chain = (int)(hom_cov * opt->high_factor + .499);
	opt->hom_cov = hom_cov;
	if (opt->max_n_chain < max_n_chain)
		opt->max_n_chain = max_n_chain;
	fprintf(stderr, "[M::%s] updated max_n_chain to %d\n", __func__, opt->max_n_chain);
}

static int check_file(char* name, const char* opt)
{
    if(!name)
    {
        fprintf(stderr, "[ERROR] file does not exist (-%s)\n", opt);
        return 0;
    } 
    FILE* is_exist = NULL;
    is_exist = fopen(name,"r");
    if(!is_exist)
    {
        fprintf(stderr, "[ERROR] %s does not exist (-%s)\n", name, opt);
        return 0;
    } 

    fclose(is_exist);
    return 1;
}

static int check_hic_reads(enzyme* f, const char* opt)
{
    int i;
    for (i = 0; i < f->n; i++)
    {
        if(check_file(f->a[i], opt) == 0) return 0;
    }
    return 1;
}

int check_option(hifiasm_opt_t* asm_opt)
{
    if(asm_opt->read_file_names == NULL || asm_opt->num_reads == 0)
    {
        fprintf(stderr, "[ERROR] missing input: please specify a read file\n");
        return 0;
    }

    if(asm_opt->output_file_name == NULL)
    {
        fprintf(stderr, "[ERROR] missing output: please specify the output name (-o)\n");
        return 0;
    }

    if(asm_opt->thread_num < 1)
    {
        fprintf(stderr, "[ERROR] the number of threads must be > 0 (-t)\n");
        return 0;
    }


    if(asm_opt->number_of_round < 1)
    {
        fprintf(stderr, "[ERROR] the number of rounds for correction must be > 0 (-r)\n");
        return 0;
    }

    if(asm_opt->clean_round < 1)
    {
        fprintf(stderr, "[ERROR] the number of rounds for assembly cleaning must be > 0 (-a)\n");
        return 0;
    }

    if(asm_opt->adapterLen < 0)
    {
        fprintf(stderr, "[ERROR] the length of removed adapters must be >= 0 (-z)\n");
        return 0;
    }


    if(asm_opt->k_mer_length >= 64)
    {
        fprintf(stderr, "[ERROR] the length of k_mer must be < 64 (-k)\n");
        return 0;
    }


    if(asm_opt->max_drop_rate < 0 || asm_opt->max_drop_rate >= 1 ) 
    {
        fprintf(stderr, "[ERROR] max overlap drop ratio must be [0.0, 1.0) (-x)\n");
        return 0;
    }

    
    if(asm_opt->min_drop_rate < 0 || asm_opt->min_drop_rate >= 1)
    {
        fprintf(stderr, "[ERROR] min overlap drop ratio must be [0.0, 1.0) (-y)\n");
        return 0;
    }

    if(asm_opt->max_drop_rate <= asm_opt->min_drop_rate)
    {
        fprintf(stderr, "[ERROR] min overlap drop ratio must be less than max overlap drop ratio (-x/-y)\n");
        return 0;
    }

    if(asm_opt->small_pop_bubble_size < 0)
    {
        fprintf(stderr, "[ERROR] the size of popped small bubbles must be >= 0 (-p)\n");
        return 0;
    }

    if(asm_opt->large_pop_bubble_size < 0)
    {
        fprintf(stderr, "[ERROR] the size of popped large bubbles must be >= 0 (-m)\n");
        return 0;
    }

    if(asm_opt->max_hang_Len < 0)
    {
        fprintf(stderr, "[ERROR] max_hang_Len must be >= 0\n");
        return 0;
    }

    if(asm_opt->max_hang_rate < 0)
    {
        fprintf(stderr, "[ERROR] max_hang_rate must be >= 0\n");
        return 0;
    }

    if(asm_opt->gap_fuzz < 0)
    {
        fprintf(stderr, "[ERROR] gap_fuzz must be >= 0\n");
        return 0;
    }

    if(asm_opt->min_overlap_Len < 0)
    {
        fprintf(stderr, "[ERROR] min_overlap_Len must be >= 0\n");
        return 0;
    }

    if(asm_opt->min_overlap_coverage < 0)
    {
        fprintf(stderr, "[ERROR] min_overlap_coverage must be >= 0\n");
        return 0;
    }

	if (asm_opt->max_ov_diff_ec < asm_opt->max_ov_diff_final) {
		fprintf(stderr, "[ERROR] max_ov_diff_ec shouldn't be smaller than max_ov_diff_final\n");
		return 0;
	}

	if (asm_opt->max_ov_diff_ec < HA_MIN_OV_DIFF) {
		fprintf(stderr, "[ERROR] max_ov_diff_ec shouldn't be smaller than %g\n", HA_MIN_OV_DIFF);
		return 0;
	}

    if(asm_opt->max_short_tip < 0)
    {
        fprintf(stderr, "[ERROR] the length of removal tips must be >= 0 (-n)\n");
        return 0;
    }

    if(asm_opt->purge_level_primary < 0 || asm_opt->purge_level_primary > 2)
    {
        fprintf(stderr, "[ERROR] the level of purge-dup should be [0, 2] (-l)\n");
        return 0;
    }

    if(ha_opt_triobin(asm_opt) && ((asm_opt->purge_level_trio < 0 || asm_opt->purge_level_trio > 1)))
    {
        fprintf(stderr, "[ERROR] the level of purge-dup for trio should be [0, 1] (-l)\n");
        return 0;
    }

    if(asm_opt->hom_global_coverage < 0 && asm_opt->hom_global_coverage != -1)
    {
        fprintf(stderr, "[ERROR] purge duplication coverage threshold should be >= 0 (--purge-cov)\n");
        return 0;
    }

    if(asm_opt->bed_inconsist_rate < 0 || asm_opt->bed_inconsist_rate > 100)
    {
        fprintf(stderr, "[ERROR] inconsistency rate should be [0, 100] (--lowQ)\n");
        return 0;
    }


    if(asm_opt->fn_bin_yak[0] != NULL && check_file(asm_opt->fn_bin_yak[0], "YAK1") == 0) return 0;
    if(asm_opt->fn_bin_yak[1] != NULL && check_file(asm_opt->fn_bin_yak[1], "YAK2") == 0) return 0;
    if(asm_opt->fn_bin_list[0] != NULL && check_file(asm_opt->fn_bin_list[0], "LIST1") == 0) return 0;
    if(asm_opt->fn_bin_list[1] != NULL && check_file(asm_opt->fn_bin_list[1], "LIST2") == 0) return 0;
    if(asm_opt->required_read_name != NULL && check_file(asm_opt->required_read_name, "b") == 0) return 0;
    
    if(asm_opt->hic_reads[0] != NULL && check_hic_reads(asm_opt->hic_reads[0], "HIC1") == 0) return 0;
    if(asm_opt->hic_reads[1] != NULL && check_hic_reads(asm_opt->hic_reads[1], "HIC2") == 0) return 0;
    if(asm_opt->hic_reads[0] != NULL && asm_opt->hic_reads[1] == NULL)
    {
        fprintf(stderr, "[ERROR] lack r2 of HiC reads (--h2)\n");
        return 0;
    }
    if(asm_opt->hic_reads[1] != NULL && asm_opt->hic_reads[0] == NULL)
    {
        fprintf(stderr, "[ERROR] lack r1 of HiC reads (--h1)\n");
        return 0;
    }

    if(asm_opt->hic_reads[0] != NULL && asm_opt->hic_reads[1] != NULL && 
                            asm_opt->hic_reads[0]->n != asm_opt->hic_reads[1]->n)
    {
        fprintf(stderr, "[ERROR] wrong r1 and r2 of HiC reads (--h1 && --h2)\n");
        return 0;
    }

    if(asm_opt->hic_enzymes != NULL && asm_opt->hic_enzymes->n == 0)
    {
        fprintf(stderr, "[ERROR] wrong HiC enzymes (--enzyme)\n");
        return 0;
    }

    if(asm_opt->hic_reads[0] != NULL && asm_opt->hic_reads[0]->n == 0)
    {
        fprintf(stderr, "[ERROR] wrong r1 of HiC reads (--h1)\n");
        return 0;
    }

    if(asm_opt->hic_reads[1] != NULL && asm_opt->hic_reads[1]->n == 0)
    {
        fprintf(stderr, "[ERROR] wrong r2 of HiC reads (--h2)\n");
        return 0;
    }

    if(asm_opt->b_low_cov < 0)
    {
        fprintf(stderr, "[ERROR] must >= 0 (--b-cov)\n");
        return 0;
    }

    if(asm_opt->b_high_cov != -1 && asm_opt->b_high_cov < 0)
    {
        fprintf(stderr, "[ERROR] must >= 0 (--h-cov)\n");
        return 0;
    }

    if(asm_opt->m_rate < 0)
    {
        fprintf(stderr, "[ERROR] must >= 0 (--m-rate)\n");
        return 0;
    }

    if(asm_opt->b_high_cov != -1 && asm_opt->b_high_cov <= asm_opt->b_low_cov)
    {
        fprintf(stderr, "[ERROR] [--h-cov] must >= [--b-cov]\n");
        return 0;
    }

    // fprintf(stderr, "input file num: %d\n", asm_opt->num_reads);
    // fprintf(stderr, "output file: %s\n", asm_opt->output_file_name);
    // fprintf(stderr, "number of threads: %d\n", asm_opt->thread_num);
    // fprintf(stderr, "number of rounds for correction: %d\n", asm_opt->number_of_round);
    // fprintf(stderr, "number of rounds for assembly cleaning: %d\n", asm_opt->clean_round);
    // fprintf(stderr, "length of removed adapters: %d\n", asm_opt->adapterLen);
    // fprintf(stderr, "length of k_mer: %d\n", asm_opt->k_mer_length);
    // fprintf(stderr, "min overlap drop ratio: %.2g\n", asm_opt->min_drop_rate);
    // fprintf(stderr, "max overlap drop ratio: %.2g\n", asm_opt->max_drop_rate);
    // fprintf(stderr, "size of popped small bubbles: %lld\n", asm_opt->small_pop_bubble_size);
    // fprintf(stderr, "size of popped large bubbles: %lld\n", asm_opt->large_pop_bubble_size);
    // fprintf(stderr, "small removed unitig threshold: %d\n", asm_opt->max_short_tip);
    // fprintf(stderr, "small removed unitig threshold: %d\n", asm_opt->max_short_tip);
    // fprintf(stderr, "min_cnt: %d\n", asm_opt->min_cnt);
    // fprintf(stderr, "mid_cnt: %d\n", asm_opt->mid_cnt);
    // fprintf(stderr, "purge_level_primary: %d\n", asm_opt->purge_level_primary);
    // fprintf(stderr, "purge_level_trio: %d\n", asm_opt->purge_level_trio);
    // fprintf(stderr, "purge_simi_rate: %f\n", asm_opt->purge_simi_rate);
    // fprintf(stderr, "purge_overlap_len: %d\n", asm_opt->purge_overlap_len);

    return 1;
}

void get_queries(int argc, char *argv[], ketopt_t* opt, hifiasm_opt_t* asm_opt)
{
    if(opt->ind == argc)
    {
        return;
    }

    asm_opt->num_reads = argc - opt->ind;
    asm_opt->read_file_names = (char**)malloc(sizeof(char*)*asm_opt->num_reads);
    
    long long i;
    gzFile dfp;
    for (i = 0; i < asm_opt->num_reads; i++)
    {
        asm_opt->read_file_names[i] = argv[i + opt->ind];
        dfp = gzopen(asm_opt->read_file_names[i], "r");
        if (dfp == 0)
        {
            fprintf(stderr, "[ERROR] Cannot find the input read file: %s\n", 
                    asm_opt->read_file_names[i]);
		    exit(0);
        }
        gzclose(dfp);
    }
}

void get_hic_enzymes(char *argv, enzyme** x, int check_name)
{
    int i, k, pre_i, len = strlen(argv);
    (*x) = (enzyme*)calloc(1, sizeof(enzyme));
    if(len == 0)
    {
        (*x)->n = 0; (*x)->l = NULL; (*x)->a = NULL;
        return;
    }


    (*x)->n = 1;
    for (i = pre_i = 0; i < len; i++)
    {
        if(argv[i] == ',')
        {
            (*x)->n++;
            continue;
        } 
        
        if(check_name)
        {
            if(argv[i] != 'A' && argv[i] != 'C' && argv[i] != 'G' && argv[i] != 'T' &&
                argv[i] != 'a' && argv[i] != 'c' && argv[i] != 'g' && argv[i] != 't' && 
                argv[i] != 'N' && argv[i] != 'n')
            {
                (*x)->n = 0;
                (*x)->l = NULL;
                (*x)->a = NULL;
                return;
            }
        }
        
    }
    (*x)->l = (int*)calloc((*x)->n, sizeof(int));
    (*x)->a = (char**)calloc((*x)->n, sizeof(char*));

    for (i = pre_i = k = 0; i < len; i++)
    {
        if(argv[i] == ',')
        {
            (*x)->l[k] = i - pre_i;
            (*x)->a[k] = (char*)malloc(sizeof(char)*((*x)->l[k]+1));
            memcpy((*x)->a[k], argv + pre_i, (*x)->l[k]);
            (*x)->a[k][(*x)->l[k]] = '\0';
            pre_i = i + 1;
            k++;
        }
    }

    (*x)->l[k] = i - pre_i;
    (*x)->a[k] = (char*)malloc(sizeof(char)*((*x)->l[k]+1));
    memcpy((*x)->a[k], argv + pre_i, (*x)->l[k]);
    (*x)->a[k][(*x)->l[k]] = '\0';
}

int CommandLine_process(int argc, char *argv[], hifiasm_opt_t* asm_opt)
{
    ketopt_t opt = KETOPT_INIT;

    int c;

    while ((c = ketopt(&opt, argc, argv, 1, "hvt:o:k:w:m:n:r:a:b:z:x:y:p:c:d:M:P:if:D:FN:1:2:3:4:l:s:O:eu", long_options)) >= 0) {
        if (c == 'h')
        {
            Print_H(asm_opt);
            return 0;
        } 
        else if (c == 'v' || c == 300)
        {
			puts(HA_VERSION);
            return 0;
        }
		else if (c == 'f') asm_opt->bf_shift = atoi(opt.arg);
        else if (c == 't') asm_opt->thread_num = atoi(opt.arg); 
        else if (c == 'o') asm_opt->output_file_name = opt.arg;
        else if (c == 'r') asm_opt->number_of_round = atoi(opt.arg);
        else if (c == 'k') asm_opt->k_mer_length = atoi(opt.arg);
        else if (c == 'i') asm_opt->load_index_from_disk = 0; 
        else if (c == 'w') asm_opt->mz_win = atoi(opt.arg);
		else if (c == 'D') asm_opt->high_factor = atof(opt.arg);
		else if (c == 'F') asm_opt->flag |= HA_F_NO_KMER_FLT;
		else if (c == 'N') asm_opt->max_n_chain = atoi(opt.arg);
        else if (c == 'a') asm_opt->clean_round = atoi(opt.arg); 
        else if (c == 'z') asm_opt->adapterLen = atoi(opt.arg);
        else if (c == 'b') asm_opt->required_read_name = opt.arg;
        else if (c == 'c') asm_opt->min_cnt = atoi(opt.arg);
        else if (c == 'd') asm_opt->mid_cnt = atoi(opt.arg);
        else if (c == '1' || c == 'P') asm_opt->fn_bin_yak[0] = opt.arg; // -P/-M reserved for backward compatibility
        else if (c == '2' || c == 'M') asm_opt->fn_bin_yak[1] = opt.arg;
        else if (c == '3') asm_opt->fn_bin_list[0] = opt.arg;
        else if (c == '4') asm_opt->fn_bin_list[1] = opt.arg;
        else if (c == 'x') asm_opt->max_drop_rate = atof(opt.arg);
        else if (c == 'y') asm_opt->min_drop_rate = atof(opt.arg);
        else if (c == 'p') asm_opt->small_pop_bubble_size = atoll(opt.arg);
        else if (c == 'm') asm_opt->large_pop_bubble_size = atoll(opt.arg);
        else if (c == 'n') asm_opt->max_short_tip = atoll(opt.arg);
        else if (c == 'e') asm_opt->flag |= HA_F_BAN_ASSEMBLY;
        else if (c == 'u') asm_opt->flag |= HA_F_BAN_POST_JOIN;
		else if (c == 301) asm_opt->flag |= HA_F_VERBOSE_GFA;
		else if (c == 302) asm_opt->flag |= HA_F_WRITE_PAF;
		else if (c == 303) asm_opt->flag |= HA_F_WRITE_EC;
		else if (c == 304) asm_opt->flag |= HA_F_SKIP_TRIOBIN;
		else if (c == 305) asm_opt->max_ov_diff_ec = atof(opt.arg);
		else if (c == 306) asm_opt->max_ov_diff_final = atof(opt.arg);
		else if (c == 307) asm_opt->extract_list = opt.arg;
		else if (c == 308) asm_opt->extract_iter = atoi(opt.arg);
        else if (c == 309) asm_opt->hom_global_coverage = atoi(opt.arg);
        else if (c == 310)
        {
            char* s = NULL;
            asm_opt->recover_atg_cov_min = strtol(opt.arg, &s, 10);
			if (*s == ',') asm_opt->recover_atg_cov_max = strtol(s + 1, &s, 10);
            if(asm_opt->recover_atg_cov_min == -1 || asm_opt->recover_atg_cov_max == -1)
            {
                asm_opt->recover_atg_cov_min = asm_opt->recover_atg_cov_max = -1;
            }
        }
        else if (c == 311) asm_opt->flag |= HA_F_HIGH_HET;
        else if (c == 312) asm_opt->bed_inconsist_rate = atoi(opt.arg);
        else if (c == 313) asm_opt->min_hist_kmer_cnt = atoi(opt.arg);
        else if (c == 314) get_hic_enzymes(opt.arg, &(asm_opt->hic_reads[0]), 0);
        else if (c == 315) get_hic_enzymes(opt.arg, &(asm_opt->hic_reads[1]), 0);
        else if (c == 316) get_hic_enzymes(opt.arg, &(asm_opt->hic_enzymes), 1);
        else if (c == 317) asm_opt->b_low_cov = atoi(opt.arg);
        else if (c == 318) asm_opt->b_high_cov = atoi(opt.arg);
        else if (c == 319) asm_opt->m_rate = atof(opt.arg);
        else if (c == 'l')
        {   ///0: disable purge_dup; 1: purge containment; 2: purge overlap
            asm_opt->purge_level_primary = asm_opt->purge_level_trio = atoi(opt.arg);
        }
        else if (c == 's') asm_opt->purge_simi_rate = atof(opt.arg);
        else if (c == 'O') asm_opt->purge_overlap_len = atoll(opt.arg);
        else if (c == ':') 
        {
			fprintf(stderr, "[ERROR] missing option argument in \"%s\"\n", argv[opt.i - 1]);
			return 1;
		} 
        else if (c == '?') 
        {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[opt.i - 1]);
			return 1;
		}
    }

    
    if (argc == opt.ind)
    {
        Print_H(asm_opt);
        return 0;
    }

    get_queries(argc, argv, &opt, asm_opt);



    return check_option(asm_opt);
}