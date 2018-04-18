# haloplex pipe

## A bioinformatics pipeline for variant calling for haloplex data.

This version written for use on Monash University M3.massive cluster by Jason Steen (jason.steen@monash.edu) , Jarred Burke and Tu Nguyen.
Based on the original verision by Khalid Mahmood (https://github.com/khalidm/hiplexpipe).

## License

See LICENSE.txt in source repository.

## Installation

#### Installation example on Monash clusters

We recommend using a virtual environment. 

```
/usr/local/python/2.7.12_static/bin/virtualenv --system-site-packages pipeline_env
module load drmaa
export DRMAA_LIBRARY_PATH=/usr/local/drmaa/1.0.7/lib/libdrmaa.so
source pipeline_env/bin/activate
pip install numpy scipy
pip install git+https://github.com/SoutheyLab/haloplex
haloplexpipe --mode map --config m3_example_pipeline.config --use_threads --log_file pipeline.log --jobs 100 --verbose 3 --just_print
```

The following lines need to be added to your .bash_profile so that VEP will function correctly.

```
PERL5LIB=$PERL5LIB:/projects/vh83/reference/annotation_databases/VEP_Modules/.vep/Plugins
export PERL5LIB
```

###Pipeline Overview

The pipeline runs in two steps. 

The first step (--mode map) will take a folder of fastq files, process them, map them to a genome, deconvolute on UMI's, generate stats and run haplotype caller on any files that meet the QC criteria.
The second step will combine all the haplotype caller files, joint genotype them, filter and annotate the vcf, and return a single combined annotated vcf for all samples.


###More detailed instructions

Create run directory

Create virtual env.  on M3 you need to use a static python binary for this to work properly. 

Load drmaa modules to allow submission to the slurm queue

Activate env

install required python modules

install haloplexpipe directly from github

copy example.config file from to run directory

confirm all example.config parameters are correct. most important are the ref genome and the two interval files.

create ./fastqs directory in run directory

create symlinks to all files required for the the project into ./fastqs directory

run the pipeline in --mode map 

once map mode is complete, run in --mode process to call variants






