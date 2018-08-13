'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such
as config, options, DRMAA and the logger.
'''

from utils import safe_make_dir
from runner import run_stage
import os
import math


def java_command(jar_path, mem_in_gb, command_args):
    '''Build a string for running a java command'''
    # Bit of room between Java's max heap memory and what was requested.
    # Allows for other Java memory usage, such as stack.
    java_mem = mem_in_gb - 2
    return 'java -Xmx{mem}g -jar {jar_path} {command_args}'.format(
        jar_path=jar_path, mem=java_mem, command_args=command_args)


def run_java(state, stage, jar_path, mem, args):
    command = java_command(jar_path, mem, args)
    run_stage(state, stage, command)


class Stages(object):
    def __init__(self, state):
        self.state = state
        self.reference = self.get_options('ref_grch37')
        self.vardict_bed = self.get_options('vardict_bed')
        self.AF_THR = self.get_options('AF_THR')
        self.dbsnp_b37 = self.get_options('dbsnp_b37')
        self.interval_file = self.get_options('interval_file')
        self.other_vep = self.get_options('other_vep')
        self.brcaex = self.get_options('vep_brcaex')
        self.gnomad = self.get_options('vep_gnomad')
        self.revel = self.get_options('vep_revel')
        self.maxentscan = self.get_options('vep_maxentscan')
        self.exac = self.get_options('vep_exac')
        self.dbnsfp = self.get_options('vep_dbnsfp')
        self.dbscsnv = self.get_options('vep_dbscsnv')
        self.cadd = self.get_options('vep_cadd')
        self.locatit_bed_file = self.get_options('locatit_bedfile')
        self.gatk_jar=self.get_options('gatk_jar')

    ######### General methods ######### 
    def run_gatk(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, self.gatk_jar, mem, args)

    def get_stage_options(self, stage, *options):
        return self.state.config.get_stage_options(stage, *options)

    def get_options(self, *options):
        return self.state.config.get_options(*options)


    ######### Map methods ######### 
    # Methods to trim, align,
    # de-mutiplex and sort input 
    # fastq files
    ###############################
    def original_fastqs(self, output):
        '''Original fastq files'''
        pass
    
    def glob_gatk(self, output):
        '''grab all the gatk .g.vcf files'''
        pass
    def grab_summary_file(self, output):
        '''grab the summary file for processing'''
        pass
    
    def run_surecalltrimmer(self, inputs, outputs, sample_id):
        '''Run SurecallTrimmer on the raw reads'''
        fastq_read1_in, fastq_read2_in = inputs
        # Need to make the path to this .jar file a variable in config file
        command = 'java -Xmx{mem}g -jar /projects/vh83/local_software/agent/SurecallTrimmer_v4.0.1.jar ' \
                                   '-fq1 {fastq_read1} -fq2 {fastq_read2} -hs -out_loc ./processed_fastqs; ' \
                                   'mv ./processed_fastqs/{sample}_R1_001.fastq.gz* ./processed_fastqs/{sample}_R1_001.processed.fastq.gz; ' \
                                   'mv ./processed_fastqs/{sample}_R3_001.fastq.gz* ./processed_fastqs/{sample}_R3_001.processed.fastq.gz' \
                  .format(mem=self.state.config.get_stage_options('run_surecalltrimmer', 'mem'),
                          fastq_read1=fastq_read1_in,
                          fastq_read2=fastq_read2_in,
                          sample=sample_id)
        run_stage(self.state, 'run_surecalltrimmer', command) 

    def align_bwa(self, inputs, bam_out, sample_id):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('align_bwa', 'cores')
        read_group = '"@RG\\tID:{sample}\\tSM:{sample}\\tPU:lib1\\tPL:Illumina"' \
            .format(sample=sample_id)
        command = 'bwa mem -M -t {cores} -R {read_group} {reference} {fastq_read1} {fastq_read2} ' \
                  '| samtools view -uh  - | samtools sort -o {bam}' \
                  .format(cores=cores,
                          read_group=read_group,
                          fastq_read1=fastq_read1_in,
                          fastq_read2=fastq_read2_in,
                          reference=self.reference,
                          bam=bam_out)
        run_stage(self.state, 'align_bwa', command)

    def run_locatit(self, inputs, bam_out):
        bam_in = [ip for ip in inputs if ip.endswith('.bam')][0]
        index_file = [ip for ip in inputs if ip.endswith('R2_001.fastq.gz')][0]
        #need to make the path to this .jar file a variable in config file
        command = 'java -Xmx{mem}G -jar /projects/vh83/local_software/agent/LocatIt_v4.0.1.jar ' \
                  '-U -q 25 -m 1 -d 0 -IB -OB -b {locatit_bed_file} ' \
                  '-o {bam_out} {bam_in} {index_file}' \
                  .format(mem=self.state.config.get_stage_options('run_locatit', 'mem'),
                         locatit_bed_file=self.locatit_bed_file,
                         bam_in=bam_in,
                         bam_out=bam_out,
                         index_file=index_file)
        run_stage(self.state, 'run_locatit', command)

    def sort_bam(self, input, bam_out):
        '''sort the locatit bam files'''
        command = 'samtools sort -o {bam_out} {input}; samtools index {bam_out}'.format(input=input, bam_out=bam_out)
        run_stage(self.state, 'sort_bam', command)

    def call_haplotypecaller_gatk(self, bam_in, vcf_out):
        '''Call variants using GATK'''
        bam_in = 'alignments/' + os.path.basename(bam_in)
        cores = self.get_stage_options('call_haplotypecaller_gatk', 'cores')
        gatk_args = "-T HaplotypeCaller -R {reference} --min_base_quality_score 20 " \
                    "--emitRefConfidence GVCF " \
                    "-A AlleleBalance -A AlleleBalanceBySample " \
                    "-A ChromosomeCounts -A ClippingRankSumTest " \
                    "-A Coverage -A DepthPerAlleleBySample " \
                    "-A DepthPerSampleHC -A FisherStrand " \
                    "-A GCContent -A GenotypeSummaries " \
                    "-A HardyWeinberg -A HomopolymerRun " \
                    "-A LikelihoodRankSumTest -A LowMQ " \
                    "-A MappingQualityRankSumTest -A MappingQualityZero " \
                    "-A QualByDepth " \
                    "-A RMSMappingQuality -A ReadPosRankSumTest " \
                    "-A SampleList -A SpanningDeletions " \
                    "-A StrandBiasBySample -A StrandOddsRatio " \
                    "-A TandemRepeatAnnotator -A VariantType " \
                    "-nct {cores} " \
                    "-I {bam} -L {interval_list} -o {out}".format(reference=self.reference,
                                                                  bam=bam_in, interval_list=self.interval_file,
                                                                  out=vcf_out, cores=cores)
        self.run_gatk('call_haplotypecaller_gatk', gatk_args)

    def run_vardict(self, bam_in, vcf_out, sample_name):
        '''call variants with vardict'''
        bam_in = 'alignments/' + os.path.basename(bam_in) 
        command = 'export PATH=/home/jste0021/scripts/git_controlled/VarDict:$PATH; ' \
                  'vardict -G {reference} -f {AF_THR} -N {sample_name} -b {bam_in} -c 1 -S 2 -E 3 -g 4 {vardict_bed} | ' \
                  'teststrandbias.R | ' \
                  'var2vcf_valid_b37_chrnames.pl -N {sample_name} -E -f {AF_THR} > {vcf_out}'.format(
                             reference=self.reference,
                             AF_THR=self.AF_THR,
                             sample_name=sample_name,
                             bam_in=bam_in,
                             vardict_bed=self.vardict_bed,
                             vcf_out=vcf_out)
        run_stage(self.state, 'run_vardict', command) 

    def sort_vcfs(self, vcf_in, vcf_out):
        '''sort vardict vcf files'''
        command = 'bcftools sort -o {vcf_out} -O z {vcf_in}'.format(vcf_out=vcf_out, vcf_in=vcf_in)
        run_stage(self.state, 'sort_vcfs', command)

    def index_vcfs(self, vcf_in, vcf_out):
        '''Index vardcit vcf files'''
        command = 'bcftools index -f --tbi {vcf_in}'.format(vcf_in=vcf_in)
        run_stage(self.state, 'index_vcfs', command)


    ################## Map metrics methods ##################
    # Methods to generate metrics on the aligned bam files,
    # check the metrics, and produce a list of bam files that 
    # pass these metrics for HaplotypeCaller to call variants on.
    #########################################################
    def generate_amplicon_metrics(self, bam_in, txt_out, sample):
        '''Generate depth information for each amplicon and sample for heatmap plotting'''
        command = 'bedtools coverage -f 5E-1 -a {bed_intervals} -b {bam_in} | ' \
                  'sed "s/$/	{sample}/g" > {txt_out}'.format(bed_intervals=self.interval_file,
                                                                bam_in=bam_in,
                                                                sample=sample,
                                                                txt_out=txt_out)
        run_stage(self.state, 'generate_amplicon_metrics', command)

    def intersect_bed(self, bam_in, bam_out):
        '''intersect the bed file with the interval file '''
        command = "intersectBed -abam {bam_in} -b {interval_file} > {bam_out} ".format(
                     bam_in=bam_in, interval_file=self.interval_file, bam_out=bam_out)
        run_stage(self.state, 'intersect_bed', command)           


    def coverage_bed(self, bam_in, txt_out):
        ''' make coverage files '''
        command = "coverageBed -b {bam_in} -a {interval_file} -hist | grep all > {txt_out}".format(
                     bam_in=bam_in, interval_file=self.interval_file, txt_out=txt_out)
        run_stage(self.state, 'coverage_bed', command)
    
    def genome_reads(self, bam_in, txt_out):
        '''count reads that map to the genome'''
        command = 'samtools view -c -F4 {bam_in} > {txt_out}'.format(
                        bam_in=bam_in, txt_out=txt_out)
        run_stage(self.state, 'genome_reads', command)

    def target_reads(self, bam_in, txt_out):
        '''count reads that map to target panel'''
        command = 'samtools view -c -F4 {bam_in} > {txt_out}'.format(
                        bam_in=bam_in, txt_out=txt_out)
        run_stage(self.state, 'target_reads', command)

    def total_reads(self, bam_in, txt_out):
        '''count the total number of reads that we started with'''
        command = 'samtools view -c {bam_in} > {txt_out}'.format(
                        bam_in=bam_in, txt_out=txt_out)
        run_stage(self.state, 'total_reads', command)

    def generate_stats(self, inputs, txt_out, samplename, joint_output):
        '''run R stats script'''
        # Assigning inputfiles to correct variables based on suffix
        for inputfile in inputs:
            if inputfile.endswith('.bedtools_hist_all.txt'):
                a = inputfile
            elif inputfile.endswith('.mapped_to_genome.txt'):
                b = inputfile
            elif inputfile.endswith('.mapped_to_target.txt'):
                c = inputfile
            elif inputfile.endswith('.total_raw_reads.txt'):
                d = inputfile
        e = samplename
        command = 'touch {txt_out}; Rscript --vanilla /projects/vh83/pipelines/code/modified_summary_stat.R ' \
                  '{hist_in} {map_genome_in} {map_target_in} {raw_reads_in} {sample_name} ' \
                  '{joint_output}'.format(txt_out=txt_out, 
                                      hist_in=a, 
                                      map_genome_in=b, 
                                      map_target_in=c, 
                                      raw_reads_in=d , 
                                      sample_name=e , 
                                      joint_output=joint_output)
        run_stage(self.state, 'generate_stats', command)

    def filter_stats(self, txt_in, txt_out):
        '''filter the summary file to make a 'passed' file'''
        # Only mark samples as pass if >= 80% of target is covered at at least 10X
        # Set to 0 for now since I want to make everything pass
        awk_comm = "{if($8 >= 0){print $1\".sorted.locatit.bam\"}}"
        command = "awk '{awk_comm}' {summary_file} > {final_file}".format(
                                        awk_comm=awk_comm,
                                        summary_file=txt_in,
                                        final_file=txt_out)
        run_stage(self.state, 'filter_stats', command)

    def read_samples(self, input_pth, outputs):
        '''Reads the list of pass samples and touches files accordingly 
        in the alignments/pass_samples folder'''
        with open(input_pth, 'r') as inputf:
            pass_files = inputf.read().split('\n')
        command_l = []
        for f in pass_files:
            command_l.append("metrics/pass_samples/{}".format(f))
        command = 'touch {}'.format(' '.join(command_l))
        run_stage(self.state, 'read_samples', command)
    
    ################## Process methods ##################
    # Methods to combine gvcf files produced by the map stage, 
    # joint-genotype, process and annotate
    #####################################################
    def combine_gvcf_gatk(self, vcf_files_in, vcf_out):
        '''Combine G.VCF files for all samples using GATK'''
        safe_make_dir('processed')
        safe_make_dir('processed/gatk')
        merge_commands = []
        temp_merge_outputs = []
        for n in range(0, int(math.ceil(float(len(vcf_files_in)) / 200.0))):
            start = n * 200
            filelist = vcf_files_in[start:start + 200]
            filelist_command = ' '.join(['--variant ' + vcf for vcf in filelist])
            temp_merge_filename = vcf_out.rstrip('.vcf') + ".temp_{start}.vcf".format(start=str(start))
            gatk_args_full = "java -Xmx{mem}g -jar {jar_path} -T CombineGVCFs -R {reference} " \
                             "--disable_auto_index_creation_and_locking_when_reading_rods " \
                             "{g_vcf_files} -o {vcf_out}; ".format(reference=self.reference, 
                                                                   jar_path=self.gatk_jar, 
                                                                   mem=self.state.config.get_stage_options('combine_gvcf_gatk', 'mem'), 
                                                                   g_vcf_files=filelist_command, 
                                                                   vcf_out=temp_merge_filename)
            merge_commands.append(gatk_args_full)
            temp_merge_outputs.append(temp_merge_filename)

        final_merge_vcfs = ' '.join(['--variant ' + vcf for vcf in temp_merge_outputs])
        gatk_args_full_final = "java -Xmx{mem}g -jar {jar_path} -T CombineGVCFs -R {reference} " \
                               "--disable_auto_index_creation_and_locking_when_reading_rods " \
                               "{g_vcf_files} -o {vcf_out}".format(reference=self.reference, 
                                                                   jar_path=self.gatk_jar, 
                                                                   mem=self.state.config.get_stage_options('combine_gvcf_gatk', 'mem'), 
                                                                   g_vcf_files=final_merge_vcfs, 
                                                                   vcf_out=vcf_out)

        merge_commands.append(gatk_args_full_final)
        final_command = ''.join(merge_commands)
        run_stage(self.state, 'combine_gvcf_gatk', final_command)

    def glob_vardict(self, output):
        '''grab all the vardict vcf files'''
        pass    

    def concatenate_vcfs(self, vcf_files_in, vcf_out):
        safe_make_dir('processed') 
        safe_make_dir('processed/vardict')
        
        merge_commands = []
        temp_merge_outputs = []
        for n in range(0, int(math.ceil(float(len(vcf_files_in)) / 200.0))):
            start = n * 200
            filelist = vcf_files_in[start:start + 200]
            filelist_command = ' '.join([vcf for vcf in filelist])
            temp_merge_filename = vcf_out.rstrip('.vcf') + ".temp_{start}.vcf".format(start=str(start))
            command1 = 'bcftools merge -O z -o {vcf_out} {join_vcf_files} && bcftools index -t -f {vcf_out}; '.format(vcf_out=temp_merge_filename, join_vcf_files=filelist_command)
            merge_commands.append(command1)
            temp_merge_outputs.append(temp_merge_filename)

        final_merge_vcfs = ' '.join([vcf for vcf in temp_merge_outputs])
        command2 = 'bcftools merge -O z -o {vcf_out} {join_vcf_files} '.format(vcf_out=vcf_out, join_vcf_files=final_merge_vcfs)

        merge_commands.append(command2)
        final_command = ''.join(merge_commands)
        run_stage(self.state, 'concatenate_vcfs', final_command)

    def genotype_gvcf_gatk(self, combined_vcf_in, vcf_out):
        '''Genotype G.VCF files using GATK'''
        cores = self.get_stage_options('genotype_gvcf_gatk', 'cores')
        gatk_args = "-T GenotypeGVCFs -R {reference} " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "--dbsnp {dbsnp} " \
                    "--num_threads {cores} --variant {combined_vcf} --out {vcf_out}" \
                    .format(reference=self.reference, dbsnp=self.dbsnp_b37,
                            cores=cores, combined_vcf=combined_vcf_in, vcf_out=vcf_out)
        self.run_gatk('genotype_gvcf_gatk', gatk_args)


    def genotype_filter_gatk(self, vcf_in, vcf_out):
        '''Apply GT filters to the genotyped VCF'''
        gatk_args = "-T VariantFiltration -R {reference} " \
                    "-V {vcf_in} " \
                    "-o {vcf_out} " \
                    "-G_filter \"g.isHetNonRef() == 1\" " \
                    "-G_filterName \"HetNonRef\" " \
                    "-G_filter \"g.isHet() == 1 && g.isHetNonRef() != 1 && " \
                    "1.0 * AD[vc.getAlleleIndex(g.getAlleles().1)] / (DP * 1.0) < 0.20\" " \
                    "-G_filterName \"AltFreqLow\" " \
                    "-G_filter \"DP < 20.0\" " \
                    "-G_filterName \"LowDP\"".format(reference=self.reference,
                                                    vcf_in=vcf_in,
                                                    vcf_out=vcf_out)
        self.run_gatk('genotype_filter_gatk', gatk_args)


    def vt_decompose_normalise(self, vcf_in, vcf_out):
        '''Decompose multiallelic sites and normalise representations'''
        command = "vt decompose -s {vcf_in} | vt normalize -r {reference} -o " \
                  "{vcf_out} -".format(reference=self.reference,
                                       vcf_in=vcf_in,
                                       vcf_out=vcf_out)
        run_stage(self.state, 'vt_decompose_normalise', command)

    def variant_annotator_gatk(self, vcf_in, vcf_out):
        '''Annotate G.VCF files using GATK'''
        cores = self.get_stage_options('variant_annotator_gatk', 'cores')
        gatk_args = "-T VariantAnnotator -R {reference} " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "-A AlleleBalance -A AlleleBalanceBySample " \
                    "-A ChromosomeCounts -A ClippingRankSumTest " \
                    "-A Coverage -A DepthPerAlleleBySample " \
                    "-A DepthPerSampleHC -A FisherStrand " \
                    "-A GCContent -A GenotypeSummaries " \
                    "-A HardyWeinberg -A HomopolymerRun " \
                    "-A LikelihoodRankSumTest " \
                    "-A MappingQualityRankSumTest -A MappingQualityZero " \
                    "-A QualByDepth " \
                    "-A RMSMappingQuality -A ReadPosRankSumTest " \
                    "-A SampleList -A SpanningDeletions " \
                    "-A StrandBiasBySample -A StrandOddsRatio " \
                    "-A TandemRepeatAnnotator -A VariantType " \
                    "--num_threads {cores} --variant {vcf_in} --out {vcf_out}" \
                    .format(reference=self.reference, 
                           cores=cores,
                           vcf_in=vcf_in,
                           vcf_out=vcf_out)
        self.run_gatk('variant_annotator_gatk', gatk_args)


    def gatk_filter(self, vcf_in, vcf_out):
        '''Filtering variants (separate filters for SNPs and indels)'''
        gatk_args = "-T VariantFiltration " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "-R {reference} " \
                    "-l ERROR " \
                    "--filterExpression \"QUAL < 30.0\" --filterName GNRL_VeryLowQual " \
                    "--filterExpression \"QD < 2.0\" --filterName GNRL_LowQD " \
                    "--filterExpression \"DP < 10\" --filterName GNRL_LowCoverage " \
                    "--filterExpression \"ReadPosRankSum < -20.0\" " \
                    "--filterName GNRL_ReadPosRankSum " \
                    "--filterExpression \"OLD_MULTIALLELIC =~ '.+'\" " \
                    "--filterName MultiAllelicSite " \
                    "--filterExpression \"VariantType == 'SNP' && MQ < 30.0\" " \
                    "--filterName SNP_LowMappingQual " \
                    "--filterExpression \"VariantType == 'SNP' && SOR > 3.0\" " \
                    "--filterName SNP_StrandBias " \
                    "--filterExpression \"VariantType == 'SNP' && MQRankSum < -12.5\" " \
                    "--filterName SNP_MQRankSum " \
                    "--filterExpression \"VariantType == 'SNP' && ReadPosRankSum < -8.0\" " \
                    "--filterName SNP_ReadPosRankSum " \
                    "--variant {vcf_in} " \
                    "-o {vcf_out}".format(vcf_in=vcf_in,
                                          vcf_out=vcf_out,
                                          reference=self.reference)
        self.run_gatk('gatk_filter', gatk_args)

    def apply_vep(self, inputs, vcf_out):
        '''Apply VEP'''
        vcf_in = inputs
        cores = self.get_stage_options('apply_vep', 'cores')
        vep_command = "vep --cache --dir_cache {other_vep} " \
                      "--assembly GRCh37 --refseq --offline " \
                      "--fasta {reference} " \
                      "--sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs --format vcf " \
                      "--vcf --force_overwrite --flag_pick --no_stats " \
                      "--custom {brcaexpath},brcaex,vcf,exact,0,Clinical_significance_ENIGMA," \
                      "Comment_on_clinical_significance_ENIGMA,Date_last_evaluated_ENIGMA," \
                      "Pathogenicity_expert,HGVS_cDNA,HGVS_Protein,BIC_Nomenclature " \
                      "--custom {gnomadpath},gnomAD,vcf,exact,0,AF_NFE,AN_NFE " \
                      "--custom {revelpath},RVL,vcf,exact,0,REVEL_SCORE " \
                      "--plugin MaxEntScan,{maxentscanpath} " \
                      "--plugin ExAC,{exacpath},AC,AN " \
                      "--plugin dbNSFP,{dbnsfppath},REVEL_score,REVEL_rankscore " \
                      "--plugin dbscSNV,{dbscsnvpath} " \
                      "--plugin CADD,{caddpath} " \
                      "--fork {cores} " \
                      "-i {vcf_in} " \
                      "-o {vcf_out}".format(other_vep=self.other_vep, 
                                            cores=cores, 
                                            vcf_out=vcf_out, 
                                            vcf_in=vcf_in, 
                                            reference=self.reference, 
                                            brcaexpath=self.brcaex, 
                                            gnomadpath=self.gnomad, 
                                            revelpath=self.revel, 
                                            maxentscanpath=self.maxentscan, 
                                            exacpath=self.exac, 
                                            dbnsfppath=self.dbnsfp, 
                                            dbscsnvpath=self.dbscsnv, 
                                            caddpath=self.cadd)
        run_stage(self.state, 'apply_vep', vep_command)

    def index_vcfs(self, vcf_in, vcf_out):                                                                             
        '''Index undr_rover vcf files'''                                                                               
        command = 'bcftools index -f --tbi {vcf_in}'.format(vcf_in=vcf_in)                                             
        run_stage(self.state, 'index_vcfs', command)      
