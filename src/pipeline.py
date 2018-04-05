'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from, regex
from stages import Stages
import glob

def make_pipeline_map(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='haloplexpipe')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Stages are dependent on the state
    stages = Stages(state)

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)

    pipeline.transform(
        task_func=stages.run_surecalltrimmer,
        name='run_surecalltrimmer',
        input=output_from('original_fastqs'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+)_R1.fastq.gz'),
        add_inputs=add_inputs('{path[0]}/{sample[0]}_R3.fastq.gz'),
        extras=['{sample[0]}'],
        # output only needs to know about one file to track progress of the pipeline, but the second certainly exists after this step.
        output='processed_fastqs/{sample[0]}_R1.processed.fastq.gz')
    

    # Align paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('run_surecalltrimmer'),
        filter=formatter('processed_fastqs/(?P<sample>[a-zA-Z0-9_-]+)_R1.processed.fastq.gz'),
        add_inputs=add_inputs('processed_fastqs/{sample[0]}_R3.processed.fastq.gz'),
        extras=['{sample[0]}'],
        output='alignments/{sample[0]}.bam')

    # Run locatit from agilent.  this should produce sorted bam files, so no sorting needed at the next step
    pipeline.collate(
        task_func=stages.run_locatit,
        name='run_locatit',
        input=output_from('align_bwa', 'original_fastqs'),
        filter=regex(r'.+/(.+_L\d\d\d).+'),
        output=r'alignments/\1.locatit.bam')
    
    pipeline.transform(
        task_func=stages.sort_bam,
        name='sort_bam',
        input=output_from('run_locatit'),
        filter=suffix('.locatit.bam'),
        output='.sorted.locatit.bam')        

    # # # # # Metrics stages # # # # #
    # generate mapping metrics (post locatit)
    pipeline.transform(
        task_func=stages.generate_amplicon_metrics,
        name='generate_amplicon_metrics',
        input=output_from('sort_bam'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).sorted.locatit.bam'),
        output='alignments/metrics/{sample[0]}.amplicon-metrics.txt',
        extras=['{sample[0]}'])

    # Intersect the bam file with the region of interest
    pipeline.transform(
        task_func=stages.intersect_bed,
        name='intersect_bed',
        input=output_from('sort_bam'),
        filter=suffix('.sorted.locatit.bam'),
        output='.intersectbed.bam')

    # Calculate coverage metrics from the intersected bam file
    pipeline.transform(
        task_func=stages.coverage_bed,
        name='coverage_bed',
        input=output_from('intersect_bed'),
        filter=suffix('.intersectbed.bam'),
        output='.bedtools_hist_all.txt')

    # Count the number of mapped reads
    pipeline.transform(
        task_func=stages.genome_reads,
        name='genome_reads',
        input=output_from('sort_bam'),
        filter=suffix('.sorted.locatit.bam'),
        output='.mapped_to_genome.txt')

    # Count the number of on-target reads
    pipeline.transform(
        task_func=stages.target_reads,
        name='target_reads',
        input=output_from('intersect_bed'),
        filter=suffix('.intersectbed.bam'),
        output='.mapped_to_target.txt')

    # Count the number of total reads
    pipeline.transform(
        task_func=stages.total_reads,
        name='total_reads',
        input=output_from('sort_bam'),
        filter=suffix('.sorted.locatit.bam'),
        output='.total_raw_reads.txt')

    # Generate summary metrics from the stats files produces
    pipeline.collate(
        task_func=stages.generate_stats,
        name='generate_stats',
        input=output_from('coverage_bed', 'genome_reads', 'target_reads', 'total_reads'), 
        #filter=regex(r'.+/(.+BS\d{4,6}.+S\d+)\..+\.txt'),
        filter=regex(r'.+/(.+)\.(bedtools_hist_all|mapped_to_genome|mapped_to_target|total_raw_reads)\.txt'),
        output=r'alignments/metrics/all_sample.summary.\1.txt',
        extras=[r'\1', 'alignments/metrics/all_sample.summary.txt'])
    # # # # # Metrics stages end # # # # #

    # # # # # Checking metrics and calling # # # # #
    # Originate to set the location of the metrics summary file
    (pipeline.originate(
        task_func=stages.grab_summary_file,
        name='grab_summary_file',
        output='alignments/metrics/all_sample.summary.txt')
            .follows('generate_stats'))

    # Awk command to produce a list of bam files passing filters
    pipeline.transform(
        task_func=stages.filter_stats,
        name='filter_stats',
        input=output_from('grab_summary_file'),
        filter=suffix('.summary.txt'),
        output='.passed.summary.txt')

    # Touch passed bams to the pass_samples folder and pass the glob of that folder to HaplotypeCaller
    pipeline.subdivide(
        name='passed_filter_files', 
        task_func=stages.read_samples,
        input=output_from('filter_stats'),
        filter=formatter(),
        output="alignments/pass_samples/*.bam")

    # Call variants using GATK
    (pipeline.transform(
        task_func=stages.call_haplotypecaller_gatk,
        name='call_haplotypecaller_gatk',
        input=output_from('passed_filter_files'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9-_]+).sorted.locatit.bam'),
        output='variants/gatk/{sample[0]}.g.vcf')
        .follows('sort_bam'))

    return(pipeline)


def make_pipeline_process(state):
    #originate process pipeline state
    
    # Define empty pipeline
    pipeline = Pipeline(name='haloplexpipe')
    # Get a list of paths to all the directories to be combined for variant calling
    run_directories = state.config.get_option('runs')
    #grab files from each of the processed directories in "runs"
    gatk_files = []
    for directory in run_directories:
        gatk_files.extend(glob.glob(directory + '/variants/gatk/*.g.vcf'))

    stages = Stages(state)

    #dummy stage to take the globbed outputs of each run that is to be processed
    pipeline.originate(
        task_func=stages.glob_gatk,
        name='glob_gatk',
        output=gatk_files)

    
    # Combine G.VCF files for all samples using GATK
    pipeline.merge(
        task_func=stages.combine_gvcf_gatk,
        name='combine_gvcf_gatk',
        input=output_from('glob_gatk'),
        output='processed/ALL.combined.vcf')

    # Genotype G.VCF files using GATK
    pipeline.transform(
        task_func=stages.genotype_gvcf_gatk,
        name='genotype_gvcf_gatk',
        input=output_from('combine_gvcf_gatk'),
        filter=suffix('.combined.vcf'),
        output='.raw.vcf')

    # Apply GT filters to genotyped vcf
    pipeline.transform(
        task_func=stages.genotype_filter_gatk,
        name='genotype_filter_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.raw.vcf'),
        output='.raw.gt-filter.vcf')

    # Decompose and normalise multiallelic sites
    pipeline.transform(
        task_func=stages.vt_decompose_normalise,
        name='vt_decompose_normalise',
        input=output_from('genotype_filter_gatk'),
        filter=suffix('.raw.gt-filter.vcf'),
        output='.raw.gt-filter.decomp.norm.vcf')

    # Annotate VCF file using GATK
    pipeline.transform(
        task_func=stages.variant_annotator_gatk,
        name='variant_annotator_gatk',
        input=output_from('vt_decompose_normalise'),
        filter=suffix('.raw.gt-filter.decomp.norm.vcf'),
        output='.raw.gt-filter.decomp.norm.annotate.vcf')

    # Filter vcf
    pipeline.transform(
        task_func=stages.gatk_filter,
        name='gatk_filter',
        input=output_from('variant_annotator_gatk'),
        filter=suffix('.raw.gt-filter.decomp.norm.annotate.vcf'),
        output='.raw.gt-filter.decomp.norm.annotate.filter.vcf')

    #Apply VEP
    pipeline.transform(
        task_func=stages.apply_vep,
        name='apply_vep',
        input=output_from('gatk_filter'),
        filter=suffix('.raw.gt-filter.decomp.norm.annotate.filter.vcf'),
        output='.raw.gt-filter.decomp.norm.annotate.filter.vep.vcf')

    return pipeline
