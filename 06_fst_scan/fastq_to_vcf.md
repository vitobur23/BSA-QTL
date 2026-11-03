# Serpentine adaptation - A. arenosa and A. lyrata FastQ to VCF (mapped to arenosa)

*Milos Duchoslav, 2026-01*

**Contents**

- [Preparations](#preparations)
  - [Original preparations](#original-preparations)
  - [Preparations before adding more samples in May 2026](#preparations-before-adding-more-samples-in-may-2026)
- [Running Jörn's pipeline](#running-jörns-pipeline)
  - [Beginning](#beginning)
  - [Up to VCF](#up-to-vcf)
  - [Problem with plotting statistics for merging of VCFs](#problem-with-plotting-statistics-for-merging-of-vcfs)
  - [Further filtration of VCF](#further-filtration-of-vcf)
  - [Adding new samples (May 2026)](#adding-new-samples-may-2026)
- [Statistics of VCF](#statistics-of-vcf)
  - [VCF - ts/tv plots](#vcf---tstv-plots)
- [Filtering bad samples](#filtering-bad-samples)
  - [Missingness per sample](#missingness-per-sample)
  - [Samples to remove](#samples-to-remove)
- [VCF filtration](#vcf-filtration)
  - [Steps 1-3](#steps-1-3)
- [VCF filtration 2 - branch with one species only](#vcf-filtration-2---branch-with-one-species-only)
  - [Prepare list of samples from one species only](#prepare-list-of-samples-from-one-species-only)
  - [Steps 5 - apply depthmask](#steps-5---apply-depthmask)


Here I make a VCF from A. arenosa and A. lyrata samples mapped to A. arenosa genome.

This is a second run of fastq to VCF, the first was using the same samples, but mapping them to A. lyrata genome.


Input data:

Sequencing data from serpentine - non-serpentine population pairs of Arabidopsis lyrata and Arabidopsis arenosa. The data, with the exception of the newly sequenced samples, are included in the [*Arabidopsis lyrata* and *Arabidopsis arenosa* resource](https://lyrata.mpipz.mpg.de/), published here:

> Glushkevich, Anna, Laura Steinmann, Nikita Tikhomirov, et al. 2026. “A Global Genomic Resource for Outcrossing Arabidopsis Lyrata and Arabidopsis Arenosa.” In bioRxiv. BioRxiv, March 3. <https://doi.org/10.64898/2026.03.02.709016>.


## Preparations

### Original preparations

#### New folder
`/path/to/arenosa_lyrata_serpentine/fastq_to_vcf`
The settings files were copied from lyrataarenosa runs and adapted.

#### Make folder for settings
```bash
mkdir -p /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa
```

I prepared the yaml file with config in a way that the following will be shared with arenosalyrata run:
- prerun scripts
- sample list
- sample read file
- trimmed fastq files

By this, I don't need to run so many steps that will be the same anyway.


#### Download reference genome

Download reference genome and annotation (annotation will be used later for extraction of 4-fold degenerated sites).

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf

mkdir arenosa_genome_grini
cd arenosa_genome_grini

wget https://github.com/PaulGrini/Arabidopsis_assemblies/raw/refs/heads/main/01_arenosa_assembly/Arabidopsis_arenosa_genome.softmasked.fna.gz
wget https://github.com/PaulGrini/Arabidopsis_assemblies/raw/refs/heads/main/01_arenosa_assembly/Arabidopsis_arenosa_genome.annotation.gtf.gz

gunzip Arabidopsis_arenosa_genome.softmasked.fna.gz
```

#### Intervals for settings for arenosa reference genome

using `seqtk`

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/arenosa_genome_grini

module load seqtk
# info about columns
seqtk comp
# info about sequences
seqtk comp Arabidopsis_arenosa_genome.softmasked.fna | head -n 50 | column -t
# sequence names
seqtk comp Arabidopsis_arenosa_genome.softmasked.fna | head -n 50 | cut -f 1

# preparation of list of unscaffolded contigs
seqtk comp Arabidopsis_arenosa_genome.softmasked.fna | tail -n +9 | cut -f 1 > /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa_unscaffolded.list

# sanity checks
# number of unscaffolded contigs
wc -l /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa_unscaffolded.list
# number of sequences (should be number of unscaffolded contigs + 8)
grep -c ">" Arabidopsis_arenosa_genome.softmasked.fna
```

### Preparations before adding more samples in May 2026

Changes in samples

- Added newly sequenced samples
    - GAU - new serpentine population
    - GLE - possible better pair for GUL, so far only 2 samples sequenced to check where it is in pop structure
    - Plants from "RNAseq serpentine experiment 2023" (BOR, VLA, GUL, ING)
- Added older samples
    - BAL - non-serpentine pop as a pair for GAU
- Removed samples
    - WEK population - It has very low sequencing depth and is therefore not usable.

#### Check searching for samples
```bash
find /path/to/ERC/raw_data -iname *7003*R1*.fastq.gz | sort
find /path/to/ERC/raw_data -iname *LITE12-2-7778*R1*.fastq.gz | sort
find /path/to/ERC/raw_data -iname *AA605A*R1*.fastq.gz | sort
```

#### Back-up of previous version of VCFs

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf
mkdir version_2026-02
mv -t version_2026-02 lyrataarenosa.merged.vcf.gz* sample_list.txt sample_list_without_bad_samples_1.txt testing_subsets

cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_filtered
mkdir version_2026-02
mv -t version_2026-02 *

cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_without_bad_samples_1
mkdir version_2026-02
mv -t version_2026-02 *

cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/depth_mask
mkdir version_2026-02
mv -t version_2026-02 *

cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/report
mkdir version_2026-02
mv -t version_2026-02 MergeSubVCFs filterbigt filterbisnp filterinvariants mask merge_bam_stats tstv_plots
cd version_2026-02/MergeSubVCFs/
rm *.pickle

```


## Running Jörn's pipeline

### Beginning

#### Run in oven (sample read file and up to VCF in one script)
This was stopped as it was using unnecessary steps.
```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub -q oven lyrataarenosa_run_snakemake_in_oven_01_sample_read_file_and_bipassed.bt.vcf.gz.sh
```

#### Unlocking directory

```bash
# Snakemake activation (version 7.32.4 installed on brno2)
source /path/to/SW/mambaforge/bin/activate Snakemake # version 7.32.4, it does not work with newer
# Go to workflow folder
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/polyploid_variant_calling/workflow

# Unlocking directory after I stopped Snakemake run
snakemake --unlock --configfile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa.yaml --profile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/cluster_profile_config/snakemake/snakemake_metacentrum
```



### Up to VCF

```bash

# Snakemake activation (version 7.32.4 installed on brno2)
source /path/to/SW/mambaforge/bin/activate Snakemake # version 7.32.4, it does not work with newer
# Go to workflow folder
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/polyploid_variant_calling/workflow

# Dry run
snakemake -n /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_filtered/lyrataarenosa.bipassed.bt.vcf.gz /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/report/merge_bam_stats/lyrataarenosa.bam.table.tsv --configfile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa.yaml --profile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/cluster_profile_config/snakemake/snakemake_metacentrum -j 250 --max-status-checks-per-second 1 --keep-incomplete --latency-wait 60 --rerun-incomplete --retries 1 --quiet --rerun-triggers {mtime,input}
 
 # Real run command (run in script in oven node)
# with "--retries 1" to re-run problematic jobs
# adding bam.table.tsv as another output
snakemake /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_filtered/lyrataarenosa.bipassed.bt.vcf.gz /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/report/merge_bam_stats/lyrataarenosa.bam.table.tsv --configfile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa.yaml --profile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/cluster_profile_config/snakemake/snakemake_metacentrum -j 250 --max-status-checks-per-second 1 --keep-incomplete --latency-wait 60 --rerun-incomplete --retries 1 --rerun-triggers {mtime,input}

# run in oven
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub -q oven lyrataarenosa_run_snakemake_in_oven_02_bipassed.bt.vcf.gz.sh
```


hapcaller was stopping for unscaffolded. I renamed lyrataarenosa_unscaffolded.lst to lyrataarenosa_unscaffolded.list and changed it in the lyrataarenosa_intervals.lst.


#### Move log files into specific folders

There were too many log files in one folder.

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/log

# make folders based on names of files preceding underscore
ls | sed -E "s/_.*$//" | uniq | xargs -I{} mkdir -p "{}"

# put the log files into these folders
for NAME in *.log*
do
FOLD=$(echo $NAME | sed -E "s/_.*$//")
mv -i $NAME ${FOLD}/
done

```

### Problem with plotting statistics for merging of VCFs

The job MergeSubVCFsbcftools is failing because of Python plotting with this error:

```
Traceback (most recent call last):
  File "/scratch/job_17179588.pbs-m1/MergeSubVCFs_arenosalyrata/combine_parse_bcftools_pickles.py", line 289, in <module>
    plt.savefig(output+"_GT_counts_biallelic.pdf",dpi=200)
  File "/path/to/SW/mambaforge/envs/polyintro_software/lib/python3.12/site-packages/matplotlib/pyplot.py", line 1135, in savefig
    fig.canvas.draw_idle()  # Need this if 'transparent=True', to reset colors.
    ^^^^^^^^^^^^^^^^^^^^^^
  File "/path/to/SW/mambaforge/envs/polyintro_software/lib/python3.12/site-packages/matplotlib/backend_bases.py", line 1893, in draw_idle
    self.draw(*args, **kwargs)
  File "/path/to/SW/mambaforge/envs/polyintro_software/lib/python3.12/site-packages/matplotlib/backends/backend_agg.py", line 383, in draw
    self.renderer = self.get_renderer()
                    ^^^^^^^^^^^^^^^^^^^
  File "/path/to/SW/mambaforge/envs/polyintro_software/lib/python3.12/site-packages/matplotlib/backends/backend_agg.py", line 398, in get_renderer
    self.renderer = RendererAgg(w, h, self.figure.dpi)
                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/path/to/SW/mambaforge/envs/polyintro_software/lib/python3.12/site-packages/matplotlib/backends/backend_agg.py", line 70, in __init__
    self._renderer = _RendererAgg(int(width), int(height), dpi)
                     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ValueError: Image size of 82500x500 pixels is too large. It must be less than 2^16 in each direction.
```

I tried to solve it by modifying the script `combine_parse_bcftools_pickles.py`. I put there maximal values for figure dimensions.

Testing the script:


```bash

# Start interactive job
qsub -I -l select=1:ncpus=1:mem=32gb:scratch_local=500gb:brno=True -l walltime=2:00:00

source /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_arenosalyrata/prerun_scripts/4_filter_bcftools.sh

cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/report/MergeSubVCFs/

cp /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa_intervals.lst .

sub_intervals=$(awk -F/ '{print $NF}' <<< /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa_intervals.lst)

awk '{print "lyrataarenosa_"$1".vcf.gz"}' $sub_intervals > input_files.list
awk '{print "lyrataarenosa_"$1".vcfstats.pickle"}' $sub_intervals > input_pickles.list

python3 combine_parse_bcftools_pickles.py --histogram_bins 50 --output lyrataarenosa --biallelic --invariants --multiallelic --pickle_list input_pickles.list


```

### Further filtration of VCF

```bash

# Snakemake activation (version 7.32.4 installed on brno2)
source /path/to/SW/mambaforge/bin/activate Snakemake # version 7.32.4, it does not work with newer
# Go to workflow folder
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/polyploid_variant_calling/workflow

# Dry run
snakemake -n /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_filtered/lyrataarenosa.merged.masked.bt.vcf.gz --configfile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa.yaml --profile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/cluster_profile_config/snakemake/snakemake_metacentrum -j 250 --max-status-checks-per-second 1 --keep-incomplete --latency-wait 60 --rerun-incomplete --retries 1 --quiet --rerun-triggers {mtime,input}
 
 # Real run command (run in script in oven node)
# with "--retries 1" to re-run problematic jobs
snakemake /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_filtered/lyrataarenosa.merged.masked.bt.vcf.gz --configfile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa.yaml --profile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/cluster_profile_config/snakemake/snakemake_metacentrum -j 250 --max-status-checks-per-second 1 --keep-incomplete --latency-wait 60 --rerun-incomplete --retries 1 --rerun-triggers {mtime,input}

# run in oven
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub -q oven lyrataarenosa_run_snakemake_in_oven_03_merged.masked.bt.vcf.gz.sh
```

### Adding new samples (May 2026)

- File `fastq_to_vcf\settings_arenosalyrata\arenosalyrata_samples.tsv` was updated.
- Some VCFs, reports etc. were backed-up (in subfolders called `version_2026-02`) - see section [Preparations before adding more samples in May 2026](#Preparations-before-adding-more-samples-in-May-2026)

#### Making new sample_read_file.tsv

This file is shared with the version with mapping to lyrata.

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub -q oven run_snakemake_in_oven_01_samples_read_file.sh
```

#### Running Snakemake pipeline up to filtered VCF

```bash

# Snakemake activation (version 7.32.4 installed on brno2)
source /path/to/SW/mambaforge/bin/activate Snakemake # version 7.32.4, it does not work with newer
# Go to workflow folder
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/polyploid_variant_calling/workflow

# Dry run
snakemake -n /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_filtered/lyrataarenosa.merged.masked.bt.vcf.gz --configfile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa.yaml --profile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/cluster_profile_config/snakemake/snakemake_metacentrum -j 250 --max-status-checks-per-second 1 --keep-incomplete --latency-wait 60 --rerun-incomplete --retries 1 --quiet --rerun-triggers {mtime,input}
 
 # Real run command (run in script in oven node)
# with "--retries 1" to re-run problematic jobs
snakemake /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_filtered/lyrataarenosa.merged.masked.bt.vcf.gz --configfile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa.yaml --profile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/cluster_profile_config/snakemake/snakemake_metacentrum -j 250 --max-status-checks-per-second 1 --keep-incomplete --latency-wait 60 --rerun-incomplete --retries 1 --rerun-triggers {mtime,input}

# run in oven
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub -q oven lyrataarenosa_run_snakemake_in_oven_03_merged.masked.bt.vcf.gz.sh
```

The joint genotyping job took up to 95 h and 235 GB of memory (for scaffold 1, which is longest).


#### Producing BAM statistics


```bash

# Snakemake activation (version 7.32.4 installed on brno2)
source /path/to/SW/mambaforge/bin/activate Snakemake # version 7.32.4, it does not work with newer
# Go to workflow folder
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/polyploid_variant_calling/workflow

# Dry run
snakemake -n /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/report/merge_bam_stats/lyrataarenosa.bam.table.tsv --configfile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa.yaml --profile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/cluster_profile_config/snakemake/snakemake_metacentrum -j 250 --max-status-checks-per-second 1 --keep-incomplete --latency-wait 60 --rerun-incomplete --retries 1 --quiet --rerun-triggers {mtime,input}
 
 # Real run command (run in script in oven node)
# with "--retries 1" to re-run problematic jobs
# adding bam.table.tsv as another output
snakemake /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/report/merge_bam_stats/lyrataarenosa.bam.table.tsv --configfile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/settings_lyrataarenosa/lyrataarenosa.yaml --profile /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/cluster_profile_config/snakemake/snakemake_metacentrum -j 250 --max-status-checks-per-second 1 --keep-incomplete --latency-wait 60 --rerun-incomplete --retries 1 --rerun-triggers {mtime,input}

# run in oven
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub -q oven lyrataarenosa_run_snakemake_in_oven_04_bam.table.tsv.sh
```





## Statistics of VCF


### VCF - ts/tv plots

Script for plotting: <https://github.com/jgerchen/polyploid_popgen/tree/main/tstv_plot>

#### Testing

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf/

mkdir testing_subsets

module load bcftools

# extract small part of vcf
time bcftools view --threads 1 -r scaffold_1:0-100000 -o testing_subsets/lyrataarenosa.merged_first_100kbp.vcf.gz -O z lyrataarenosa.merged.vcf.gz


cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf/testing_subsets/

# this module is probably loaded automatically on Zenith frontend
module load python/3.7.7-intel-19.0.4-mgiwa7z

# download Joern's script
wget https://raw.githubusercontent.com/jgerchen/polyploid_popgen/refs/heads/main/tstv_plot/plot_tstv.py

# run script
python3 plot_tstv.py --vcf lyrataarenosa.merged_first_100kbp.vcf.gz

```

#### Preparation for bigger run

```bash
## download script

mkdir -p /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/scripts/tstv_plot

cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/scripts/tstv_plot

# download Joern's script
wget https://raw.githubusercontent.com/jgerchen/polyploid_popgen/refs/heads/main/tstv_plot/plot_tstv.py
wget https://raw.githubusercontent.com/jgerchen/polyploid_popgen/refs/heads/main/tstv_plot/README.md


## make folder for results

mkdir -p /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/report/tstv_plots

cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/report/tstv_plots

```

```bash
### Script for Metacentrum
## Plotting ts/tv ratio based on various variables to find the right filtering values for VCF

#PBS -N tstv_plots
#PBS -l select=1:ncpus=1:mem=64gb:scratch_local=1000gb
#PBS -l walltime=24:00:00 
#PBS -m ae

# input vcf
VCF="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_filtered/lyrataarenosa.bisel.bt.vcf.gz"
# folder for output
OUTDIR="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/report/tstv_plots"
# script for plotting
plot_tstv="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/scripts/tstv_plot/plot_tstv.py"


# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# loading modules needed
module load python/3.7.7-intel-19.0.4-mgiwa7z
module load py-numpy/1.19.0-intel-19.0.4-oq5edni
module load py-matplotlib/3.2.2-intel-19.0.4-ktqvh4e

echo "Modules loaded" | ts '[%Y-%m-%d %H:%M:%S]'

# Make subfolder for results based on VCF name
new_folder=$(basename $VCF .vcf.gz)

cd $OUTDIR

mkdir -p $new_folder

cd $new_folder

echo "Folder for results prepared" | ts '[%Y-%m-%d %H:%M:%S]'

# run Joern's script fro plotting
python3 $plot_tstv --vcf $VCF

echo "Python script finished" | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch

# Resources: 26 min, 99 % CPU time, 26 GB memory
```


```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub lyrataarenosa_tstv_01_bisel.sh
```


## Filtering bad samples

### Missingness per sample

Calculate missingness per sample in `lyrataarenosa.bigt.dp.bt.vcf.gz`. This is biallelic VCF after filtration for depth of genotypes but before filtration for missingness per site.

I will get the statistics from `lyrataarenosa.bigt.dp.merged.tsv`.

```r
setwd("/path/to/arenosa_lyrata_serpentine")

vcf.stats <- read.table(file = "fastq_to_vcf/lyrataarenosa_results/report/filterbigt/lyrataarenosa.bigt.dp.merged.tsv", sep = "\t", header = T)
# vcf.stats <- read.table(file = "fastq_to_vcf/lyrataarenosa_results/report/filterbigt/lyrataarenosa.bigt.dp.m.merged.tsv", sep = "\t", header = T)
# vcf.stats <- read.table(file = "fastq_to_vcf/lyrataarenosa_results/report/filterbigt/lyrataarenosa.bigt.merged.bt.tsv", sep = "\t", header = T)

summary(vcf.stats)
head(vcf.stats)

unique(vcf.stats$statistic_type)

gt.count <- vcf.stats[vcf.stats$statistic_type == "GT_count", ]

# Do the genotype counts sum to the same number for all saples?
gt.sum <- tapply(X = gt.count$statistic, INDEX = gt.count$sample, FUN = sum)
gt.sum

# only missing genotypes
miss.count <- gt.count[grep(pattern = ".", x = gt.count$var_or_bin_start, fixed = T), ]
# sum per sample
miss.sum <- tapply(X = miss.count$statistic, INDEX = miss.count$sample, FUN = sum)

# are the names of samples same?
all(names(gt.sum) == names(miss.sum))

# frequency of missing genotypes per sample
f.miss <- miss.sum/gt.sum
f.miss

pop <- substr(x = names(f.miss), start = 1, 3)

### plot
# pdf ("basic_popgen/arenosalyrata_output/arenosalyrata_missing data.pdf", width=14, height=7, onefile = T)

barplot(sort(f.miss, decreasing = T)[1:20], las = 2, 
        main = "Missing data per sample (samples with most missing data)")
abline(h = 0.4, col = 2)

# missing data per pop
mis.pop <- tapply(X = f.miss, INDEX = pop, FUN = mean)
# barplot(sort(mis.pop, decreasing = T), las = 2, main = "Missing data per population")

# Plot showing all samples and pop means
plot.default(x = jitter(as.integer(as.factor(pop))), y = f.miss, xaxt = "n", type = "n",
             xlab = "Population", ylab = "Proportion of missing data", main = "Proportion of missing data per sample",
             lwd =  2)
points(x = 1:nlevels(as.factor(pop)), y = mis.pop, pch = 19, col = 2, cex = 1.5)
text(x = jitter(as.integer(as.factor(pop))), y = f.miss, 
     labels = as.integer(gsub(pattern = "^....|..$", replacement = "", x = names(f.miss))))
abline(v = 2:nlevels(as.factor(pop)) - 0.5, col = "grey")
axis(side = 1, at = 1:nlevels(as.factor(pop)), labels = levels(as.factor(pop)), las = 2)

# boxplot(f.miss ~ as.factor(pop), las = 2)
# dev.off()

# samples to discard
f.miss[f.miss > 0.35]
```

Samples to discard (old VCF, February 2026, threshold 0.4)
```
 GUL-20ta  HOC-06ta  MIL-05da  WEK-01ta  WEK-04ta  WEK-05ta  WEK-08ta 
0.7171502 0.4004579 0.5047806 0.4017787 0.4088638 0.5958083 0.5960117
```

Samples to discard (new VCF, May 2026, threshold 0.35)

```
 GUL-20ta  HOC-06ta  LOM-11dl  MIL-05da  PEK-15tl 
0.7191103 0.3982357 0.3911670 0.5065613 0.3692828 
```

### Samples to remove

To remove

- whole population WEK
    - I don't need it and it has a lot of missing data (low sequencing depth)
- samples GUL-20ta, HOC-06ta, LOM-11dl,  MIL-05da,  PEK-15tl
    - Too much missing data
- sample PEK-02tl
    - Strange position in PCA and NeighborNet


```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf

mkdir -p ../vcf_without_bad_samples_1

module load bcftools

# list samples
bcftools query -l lyrataarenosa.merged.vcf.gz > sample_list.txt

grep -E -v 'WEK|GUL-20ta|HOC-06ta|LOM-11dl|MIL-05da|PEK-15tl|PEK-02tl' sample_list.txt > sample_list_without_bad_samples_1.txt

wc -l sample_list.txt
wc -l sample_list_without_bad_samples_1.txt

# testing the command
bcftools view -r scaffold_1:0-10000 -S sample_list_without_bad_samples_1.txt -o ../vcf_without_bad_samples_1/lyrataarenosa.merged.vcf.gz -O z lyrataarenosa.merged.vcf.gz

# check samples
bcftools query -l ../vcf_without_bad_samples_1/lyrataarenosa.merged.vcf.gz
```

#### Script for samples removal

```bash
### Script for Metacentrum
## Removing bad samples from VCF

#PBS -N lyrataarenosa_remove_bad_samples_from_vcf_1
#PBS -l select=1:ncpus=3:mem=64gb:scratch_local=1000gb:brno=True
#PBS -l walltime=24:00:00 
#PBS -m ae

# input vcf
vcf_in="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf/lyrataarenosa.merged.vcf.gz"
# output vcf
vcf_out="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_without_bad_samples_1/lyrataarenosa.merged.vcf.gz"
# samples to be retained
sample_list="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf/sample_list_without_bad_samples_1.txt"


# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# loading modules needed
module load bcftools

echo "Modules loaded" | ts '[%Y-%m-%d %H:%M:%S]'

# Make folder for output if it does not exist
mkdir -p $(dirname $vcf_out)

echo "Folder for results prepared" | ts '[%Y-%m-%d %H:%M:%S]'

# filter vcf
# --threads 1 adds one thread for compression (so 2 CPUs are needed)
bcftools view --threads 2 -S $sample_list -o $vcf_out -O z $vcf_in || { echo >&2 "Error in bcftools!"; exit 2; }
# index new vcf
tabix $vcf_out

echo "Filtering finished" | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch

# Resources: 8 h, 72 % CPU, 100 % memory
# Version with 2 CPUs was much slower.
```

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub lyrataarenosa_remove_bad_samples_from_vcf_1.sh
```


## VCF filtration

- Filtering v.02
    1. Put tag in filter column to originally multiallelic sites (in order to be able to filter them later if needed).
    2. Turn to no-call genotypes with low depth (FMT/DP), remove unobserved alleles and remove invariants
    	- This might slightly lower number of alleles per site
    	- Removing invariants reduces the file size a they are probably not needed
    	- The resulting VCF might be good for heatmaps
    3. Remove invariants, MNPs and indels, keep only SNPs (the others are not used by Piawka anyway)
    	- Try Piawka also after this step
    	- Is there some interesting variability in regions where only some populations have data?
    4. Filter sites with too much missing data (gen_max_missing: 0.5)
    5. Make and apply depthmask
    	- This is probably the best VCF for Piawka.
    

### Steps 1-3

#### Script for filtration

```bash
### Script for Metacentrum
## Filtration of VCF, steps 01-03

#PBS -N lyrataarenosa_filtering_vcf_01_03
#PBS -l select=1:ncpus=3:mem=64gb:scratch_local=1000gb:brno=True
#PBS -l walltime=24:00:00 
#PBS -m ae

gen_min_depth=3 # if 3 is set, it means that 4 and higher will be retained
vcf_start="lyrataarenosa.merged.vcf.gz"
input_dir="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_without_bad_samples_1/"
output_dir="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_without_bad_samples_1/"

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# loading modules needed
module load bcftools

echo "Modules loaded" | ts '[%Y-%m-%d %H:%M:%S]'

# Make folder for output if it does not exist
mkdir -p $output_dir

echo "Folder for results prepared" | ts '[%Y-%m-%d %H:%M:%S]'

cd $SCRATCHDIR


# VCF filename without suffix
vcf_base=$(basename $vcf_start .vcf.gz)

### filter vcf
# --threads 2 adds two threads for compression (so 3 CPUs are needed)

## Filtering step 01 - Tag sites that were originally multiallelic
time bcftools filter --threads 2 -m+ -s'MultiAllelic' -e'N_ALT>1' -o ${vcf_base}_01_tagMultiAllelic.vcf.gz ${input_dir}${vcf_base}.vcf.gz || { echo >&2 "Error in bcftools!"; exit 2; }

echo "Filtering step 01 done" | ts '[%Y-%m-%d %H:%M:%S]'

## Filtering step 02 - Turn to no-call genotypes with low depth (FMT/DP)
# Fix header - without that the next step does not work
# Extract the header
bcftools view -h ${vcf_base}_01_tagMultiAllelic.vcf.gz > ${vcf_base}_01_tagMultiAllelic.header
# Change the header
sed -i 's/ID=PL,Number=G/ID=PL,Number=./g' ${vcf_base}_01_tagMultiAllelic.header
# Filtering pipeline
time bcftools reheader -h ${vcf_base}_01_tagMultiAllelic.header ${vcf_base}_01_tagMultiAllelic.vcf.gz | # Fix the header (put there the changed one)
	bcftools filter -i "FMT/DP>$gen_min_depth" --set-GTs . --threads 2 -Ou | # Turn to no-call genotypes with low depth (FMT/DP)
	bcftools +fill-AN-AC -Ou | # Recalculate AN and AC
	bcftools view -a -Ou | # Remove alleles that are not in the genotypes
	bcftools view -m2 --threads 2 -o ${output_dir}${vcf_base}_02_dp.vcf.gz || # Remove invariants to decrease file size
	{ echo >&2 "Error in bcftools!"; exit 3; }
	
tabix ${output_dir}${vcf_base}_02_dp.vcf.gz

echo "Filtering step 02 done" | ts '[%Y-%m-%d %H:%M:%S]'

# Filtering 03 - Keep only SNPs - remove indels, MNPs and sites overlapping deletions (*)
# It might be possible to obtain more SNPs by turning genotypes with * alleles (overlapping deletions) and indels to no-call. But I did not find a simple way how to do that.
time bcftools view -m2 -e 'ALT="*"' --exclude-types indels,mnps --threads 2 -o ${output_dir}${vcf_base}_03_dp.snp.vcf.gz ${output_dir}${vcf_base}_02_dp.vcf.gz  || { echo >&2 "Error in bcftools!"; exit 4; }

tabix ${output_dir}${vcf_base}_03_dp.snp.vcf.gz

echo "Filtering step 03 done" | ts '[%Y-%m-%d %H:%M:%S]'
echo ""

# Checking VCF produced in step 3
echo "Values in FILTER column after step 03 (${output_dir}${vcf_base}_03_dp.snp.vcf.gz):"
# Check the values in FILTER column
bcftools query -f '%FILTER\n' ${output_dir}${vcf_base}_03_dp.snp.vcf.gz | sort | uniq -c | sort -nr

echo ""
echo "Checking of VCF after step 03 done" | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch

# Resources: 19 h, 66 % CPU, 100 % memory
```

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub lyrataarenosa_filtering_vcf_01_03.sh
```

#### Checking VCF produced

```
Values in FILTER column after step 03 (/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_without_bad_samples_1/lyrataarenosa.merged_03_dp.snp.vcf.gz):
21147690 PASS
4551414 MultiAllelic

```
```r
n.pass = 21147690
n.multi = 4551414
# fraction of multiallelic or originally multiallelic sites
n.multi/(n.multi + n.pass)
# total sites
n.pass + n.multi
```

Multiallelic sites (or originally multiallelic) make 17.7 % of the VCF with SNPs only.



## VCF filtration 2 - branch with one species only

### Prepare list of samples from one species only

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_without_bad_samples_1

mkdir -p ../vcf_one_species_2

module load bcftools

# list samples
bcftools query -l lyrataarenosa.merged_03_dp.snp.vcf.gz > sample_list_without_bad_samples_1.txt

# list populations
sed 's/-.*//' sample_list_without_bad_samples_1.txt | uniq
sed 's/-.*//' sample_list_without_bad_samples_1.txt | uniq | paste -s -d"|"

# keep only populations of one species
# remove lyrata populations and GLE and BDO populations which was included just for checking neigbhor-net
grep -E -v 'AHM|BRA|KHM|LOM|NKM|PEK|SPN|SUN|TCM|TEM|GLE|BDO' sample_list_without_bad_samples_1.txt > ../vcf_one_species_2/sample_list_one_species_2.txt

wc -l sample_list_without_bad_samples_1.txt
wc -l ../vcf_one_species_2/sample_list_one_species_2.txt

```

#### Script for samples removal, missingness filtration (step 4) and depthmask calculation

1. Remove samples from the other species 
    - also remove alleles that are not in the genotypes and then invariants
2. Step 04 of filtration: Remove sites with high missingness.
3. Calculate depthmask
    - Make mask for sites with excessive depth using [Jörn's script "make_depth_mask.py"](https://github.com/jgerchen/polyploid_popgen/tree/main/depth_mask). These sites can be e.g. repetitive sequences collapsed in the reference genome.

```bash
### Script for Metacentrum
## Removing samples of the other species from VCF, filtering for missingness and calculating depthmask

#PBS -N lyrataarenosa_one_species_2_samples_removal_missingness_and_depthmask
#PBS -l select=1:ncpus=3:mem=64gb:scratch_local=1000gb:brno=True
#PBS -l walltime=24:00:00 
#PBS -m ae

# input vcf
vcf_in="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_without_bad_samples_1/lyrataarenosa.merged_03_dp.snp.vcf.gz"
# base for vcf name
vcf_base="lyrataarenosa.one_sp"
output_dir="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/"
# samples to be retained
sample_list="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/sample_list_one_species_2.txt"

## settings for missingness filtering
# Maximal fraction of samples with missing genotype
gen_max_missing=0.5

## Settings for depthmask calculation
# This is the '-m, --min_depth' parameter for Python script ("Minimum depth below which a site is not considered for calculation of mean depth and Stderror"). I set it the same as for filtering - these sites should be already filtered so everything is used for calculation.
gen_min_depth=3 # meaning that 4 and higher will be retained
# Proportion of samples with excess depth required for site to be included in depth mask
depthmask_prop_excess_depth=0.5
# number of sites for depthmask script will be calculated

output_dir_depthmask="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/depthmask/"
# Joern's script for depth mask calculation (https://github.com/jgerchen/polyploid_popgen/tree/main/depth_mask)
depth_mask_script="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/polyploid_variant_calling/workflow/scripts/make_depth_mask.py"


# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# loading modules needed
module load bcftools

echo "Modules loaded" | ts '[%Y-%m-%d %H:%M:%S]'

# Make folder for output if it does not exist
mkdir -p ${output_dir}

echo "Folder for results prepared" | ts '[%Y-%m-%d %H:%M:%S]'

## Samples removal

# filter vcf
# --threads 1 adds one thread for compression (so 2 CPUs are needed)
bcftools view --threads 2 -S $sample_list -Ou $vcf_in | # Remove unwanted samples
	bcftools +fill-AN-AC -Ou | # Recalculate AN and AC
	bcftools view -a -Ou | # Remove alleles that are not in the genotypes
	bcftools view -m2 --threads 2 -o ${output_dir}${vcf_base}_03_dp.snp.vcf.gz || # Remove invariants
	{ echo >&2 "Error in bcftools!"; exit 2; }
# index new vcf
tabix ${output_dir}${vcf_base}_03_dp.snp.vcf.gz

echo "Removal of samples done" | ts '[%Y-%m-%d %H:%M:%S]'
echo ""

# Checking VCF produced
echo "Values in FILTER column after samples removal (${output_dir}${vcf_base}_03_dp.snp.vcf.gz):"
# Check the values in FILTER column
bcftools query -f '%FILTER\n' ${output_dir}${vcf_base}_03_dp.snp.vcf.gz | sort | uniq -c | sort -nr

echo ""
echo "Checking of VCF after samples removal done" | ts '[%Y-%m-%d %H:%M:%S]'


## Filtering step 04 - remove sites with fraction of samples with missing genotypes higher then $gen_max_missing

bcftools view -i "F_MISSING<=$gen_max_missing" --threads 2 -o ${output_dir}${vcf_base}_04_dp.snp.m.vcf.gz ${output_dir}${vcf_base}_03_dp.snp.vcf.gz || { echo >&2 "Error in bcftools!"; exit 2; }

echo "Filtering step 04 done" | ts '[%Y-%m-%d %H:%M:%S]'

tabix ${output_dir}${vcf_base}_04_dp.snp.m.vcf.gz

echo "Indexing done" | ts '[%Y-%m-%d %H:%M:%S]'

echo ""

# Checking VCF produced in step 4
echo "Values in FILTER column after step 04 (${output_dir}${vcf_base}_04_dp.snp.m.vcf.gz):"
# Check the values in FILTER column
bcftools query -f '%FILTER\n' ${output_dir}${vcf_base}_04_dp.snp.m.vcf.gz | sort | uniq -c | sort -nr

echo ""
echo "Checking of VCF after step 04 done" | ts '[%Y-%m-%d %H:%M:%S]'

echo "Filtering finished" | ts '[%Y-%m-%d %H:%M:%S]'


## Depthmask calculation

# unload bcftools to be safer that the right python is used (maybe it is not needed)
module unload bcftools

# Make folder for output if it does not exist
mkdir -p $output_dir_depthmask

echo "Folder for depthmask prepared" | ts '[%Y-%m-%d %H:%M:%S]'

# calculate number of sites (use, if it is not set directly)
n_loci=$(zcat ${output_dir}${vcf_base}_04_dp.snp.m.vcf.gz | grep -c "^[^#]")

echo "Number of sites in VCF: $n_loci"

# calculate number of samples
n_samples=$(echo $(zgrep -m 1 "^#CHROM" ${output_dir}${vcf_base}_04_dp.snp.m.vcf.gz | head -n 1 | wc -w) - 9 | bc)
echo "Number of samples in VCF: $n_samples"
# calculate min number of samples (from fraction) with excess depth for site to be masked
depthmask_n_samples_excess_depth=$(echo $n_samples '*' $depthmask_prop_excess_depth | bc)
echo "Number of samples with excess depth for site to be masked: $depthmask_n_samples_excess_depth"

### Depth mask calculation
python3 $depth_mask_script -v ${output_dir}${vcf_base}_04_dp.snp.m.vcf.gz -o ${output_dir_depthmask}${vcf_base}_04_dp.snp.m.depthmask.tsv -c ${output_dir_depthmask}${vcf_base}_04_dp.snp.m.depthmask_counts_out.tsv -p ${output_dir_depthmask}${vcf_base}_04_dp.snp.m.depthmask_hist_out.tsv -n $depthmask_n_samples_excess_depth -l $n_loci -m $gen_min_depth > ${output_dir_depthmask}${vcf_base}_04_dp.snp.m.depthmask_script_output.txt

# Conversion to bed
awk '{printf "%s\t%d\t%d\n" ,$1,($2 - 1),$2}' ${output_dir_depthmask}${vcf_base}_04_dp.snp.m.depthmask.tsv > ${output_dir_depthmask}${vcf_base}_04_dp.snp.m.depthmask.bed


echo "Depthmask created." | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch

# Resources: 3 h, 59 % CPU, 100 % memory
```

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub lyrataarenosa_one_species_2_samples_removal_missingness_and_depthmask.sh
```

#### Output of the script
```
[2026-06-11 13:38:06] Modules loaded
[2026-06-11 13:38:06] Folder for results prepared
[2026-06-11 14:38:59] Removal of samples done

Values in FILTER column after samples removal (/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/lyrataarenosa.one_sp_03_dp.snp.vcf.gz):
18166856 PASS
4371738 MultiAllelic

[2026-06-11 14:48:10] Checking of VCF after samples removal done
[2026-06-11 15:32:24] Filtering step 04 done
[2026-06-11 15:44:05] Indexing done

Values in FILTER column after step 04 (/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/lyrataarenosa.one_sp_04_dp.snp.m.vcf.gz):
16823282 PASS
4023676 MultiAllelic

[2026-06-11 15:52:46] Checking of VCF after step 04 done
[2026-06-11 15:52:46] Filtering finished
[2026-06-11 15:52:47] Folder for depthmask prepared
Number of sites in VCF: 20846958
Number of samples in VCF: 209
Number of samples with excess depth for site to be masked: 104.5
[2026-06-11 16:56:25] Depthmask created.

```

#### Checking histogram produced by depthmask script

Column 1: How many samples are above threshold of depth.
Column 2: How many sites have this many samples above threshold of depth.

```r
setwd("/path/to/arenosa_lyrata_serpentine/")
if(!exists(x = "old.par")) old.par<-par(no.readonly = T)

mask.hist <- read.table(file = "fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/depthmask/lyrataarenosa.one_sp_04_dp.snp.m.depthmask_hist_out.tsv", header = F, sep = "\t")
colnames(mask.hist) <- c("n_samples", "n_sites")

plot(x = mask.hist[, 1], y = mask.hist[, 2], type = "l", ylab = "# sites above cut-off", xlab = "# samples above cut-off")
plot(x = mask.hist[, 1], y = mask.hist[, 2], type = "b", log = "y", ylab = "# sites above cut-off", xlab = "# samples above cut-off")

# summing sites that will be masked for various thresholds
mask.hist$summed_n_sites <- NA
for(i in nrow(mask.hist):1) {
	mask.hist$summed_n_sites[i] <- sum(mask.hist$n_sites[i:nrow(mask.hist)])
}

plot(x = mask.hist[, 1], y = mask.hist[, 3], type = "b", ylab = "# sites above cut-off (summed)", xlab = "# samples above cut-off", xlim = rev(range(mask.hist[, 1])))

# plot to be exported
pdf(file = paste0("fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/depthmask/depthmask_cutoff_plot", ".pdf"),
    paper = "a4r")
plot(x = mask.hist[, 1], y = mask.hist[, 3], type = "b", log = "y", ylab = "# sites above cut-off (summed)", xlab = "# samples above cut-off", xlim = rev(range(mask.hist[, 1])))
# setting for making depthmask
# Proportion of samples with excess depth required for site to be included in depth mask
depthmask_prop_excess_depth=0.5
abline(v = max(mask.hist$n_samples)*depthmask_prop_excess_depth, col = 2)
abline(h = mask.hist$summed_n_sites[max(mask.hist$n_samples)*depthmask_prop_excess_depth + 1], col = 2)
mtext(paste0("Number of samples with excess depth required for site to be included in depthmask: ", max(mask.hist$n_samples)*depthmask_prop_excess_depth),
			line = 1)
mtext(paste0("Number of sites to be masked: ", mask.hist$summed_n_sites[max(mask.hist$n_samples)*depthmask_prop_excess_depth + 1], 
						 " out of ", sum(mask.hist$n_sites), 
			" (", round(mask.hist$summed_n_sites[max(mask.hist$n_samples)*depthmask_prop_excess_depth + 1]/sum(mask.hist$n_sites)*100, 2), " %)"))
dev.off()
```


#### Merging sites to bigger regions in depthmask

Script for depthmask is calculating the depth for every site, not for regions. I will use bedtools to merge nearby sites to regions. I have chosen to merge sites which are maximally 500 bp apart.

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/depthmask/

module load bedtools2

# number of sites
wc -l lyrataarenosa.one_sp_04_dp.snp.m.depthmask.bed # 147077
# checking number of masked regions for various merging settings
bedtools merge -d 0 -i lyrataarenosa.one_sp_04_dp.snp.m.depthmask.bed | wc -l # 80935
bedtools merge -d 100 -i lyrataarenosa.one_sp_04_dp.snp.m.depthmask.bed | wc -l # 2431
bedtools merge -d 200 -i lyrataarenosa.one_sp_04_dp.snp.m.depthmask.bed | wc -l # 1697
bedtools merge -d 500 -i lyrataarenosa.one_sp_04_dp.snp.m.depthmask.bed | wc -l # 1284
bedtools merge -d 1000 -i lyrataarenosa.one_sp_04_dp.snp.m.depthmask.bed | wc -l # 1084

# merging regions
bedtools merge -d 500 -i lyrataarenosa.one_sp_04_dp.snp.m.depthmask.bed > lyrataarenosa.one_sp_04_dp.snp.m.depthmask_merged.bed
```

#### Comparison of depthmasks for one species and for both species

```r
setwd("/path/to/arenosa_lyrata_serpentine/")
if(!exists(x = "old.par")) old.par<-par(no.readonly = T)

## depthmasks without merged regions

# depthmask calculated with data from both species
dm.both <- read.table(file = "fastq_to_vcf/lyrataarenosa_results/vcf_without_bad_samples_1/depthmask/lyrataarenosa.merged_04_dp.snp.m.depthmask.bed", header = F, sep = "\t")
# depthmask calculated with data from one species only
dm.one <- read.table(file = "fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/depthmask/lyrataarenosa.one_sp_04_dp.snp.m.depthmask.bed", header = F, sep = "\t")

colnames(dm.both) <- c("chrom", "start", "end")
colnames(dm.one) <- c("chrom", "start", "end")

# get all scaffolds occuring in at least one of the bed files
lev.chrom.pre <- unique(c(dm.both$chrom, dm.one$chrom))
lev.chrom <- lev.chrom.pre[order(as.integer(sub(pattern = "scaffold_", replacement = "", x = lev.chrom.pre)))]
# set the levels of the factors
dm.both$chrom <- factor(dm.both$chrom, levels = lev.chrom)
dm.one$chrom <- factor(dm.one$chrom, levels = lev.chrom)

summary(dm.both)

# plot the sites
par(mar = c(8, 4, 4, 2) + 0.1)
plot.default(x = as.integer(dm.both$chrom) - 0.2, y = dm.both$end, xaxt = "n", pch = 3)
axis(side = 1, at = seq_along(lev.chrom), labels = lev.chrom, las = 2)
points(x = as.integer(dm.one$chrom) + 0.2, y = dm.one$end, xaxt = "n", pch = 3, col = 2)
par(old.par)


## depthmasks with merged regions

# depthmask calculated with data from both species
dm.both <- read.table(file = "fastq_to_vcf/lyrataarenosa_results/vcf_without_bad_samples_1/depthmask/lyrataarenosa.merged_04_dp.snp.m.depthmask_merged.bed", header = F, sep = "\t")
# depthmask calculated with data from one species only
dm.one <- read.table(file = "fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/depthmask/lyrataarenosa.one_sp_04_dp.snp.m.depthmask_merged.bed", header = F, sep = "\t")

colnames(dm.both) <- c("chrom", "start", "end")
colnames(dm.one) <- c("chrom", "start", "end")

# get all scaffolds occuring in at least one of the bed files
lev.chrom.pre <- unique(c(dm.both$chrom, dm.one$chrom))
lev.chrom <- lev.chrom.pre[order(as.integer(sub(pattern = "scaffold_", replacement = "", x = lev.chrom.pre)))]
# set the levels of the factors
dm.both$chrom <- factor(dm.both$chrom, levels = lev.chrom)
dm.one$chrom <- factor(dm.one$chrom, levels = lev.chrom)

summary(dm.both)

# plot the regions

pdf(file = paste0("fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/depthmask/depthmask_comparison_both_or_one_species", ".pdf"), width = 10, height = 20)
par(mar = c(8, 4, 4, 2) + 0.1)
plot.default(x = as.integer(dm.both$chrom) - 0.2, y = dm.both$end, xaxt = "n", pch = 3, type = "n", main = "Regions masked by depthmask")
axis(side = 1, at = seq_along(lev.chrom), labels = lev.chrom, las = 2)
rect(xleft = as.integer(dm.both$chrom) - 0.3, xright = as.integer(dm.both$chrom) - 0, ybottom = dm.both$start, ytop = dm.both$end, col = 1, border = 1)
rect(xleft = as.integer(dm.one$chrom) + 0.3, xright = as.integer(dm.one$chrom) + 0, ybottom = dm.one$start, ytop = dm.one$end, col = 2, border = 2)
legend("topright", legend = c("data from both species", "data from one species only"), fill = c(1, 2))
par(old.par)
dev.off()

```

Result:

The depthmasks differ in many regions. It is probably better to calculate them for each species separately.

I expected that for arenosa it will be more similar to the depthmask calculated for all samples, because the most of the samples is from arenosa. However, the depthmask calculated for arenosa also differs.

### Steps 5 - apply depthmask


```bash
### Script for Metacentrum
## Filtration of VCF, step 05

#PBS -N lyrataarenosa_one_species_2_filtering_vcf_05_apply_depthmask
#PBS -l select=1:ncpus=3:mem=64gb:scratch_local=1000gb:brno=True
#PBS -l walltime=24:00:00 
#PBS -m ae

vcf_base="lyrataarenosa.one_sp"
# To be compatible with previous script, the input file ("${output_dir}${vcf_base}_03_dp.snp.vcf.gz") will be taken from the $output_dir and the output file will be put in the same folder.
# input_dir="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_without_bad_samples_1/"
output_dir="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/"
depthmask="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/depthmask/lyrataarenosa.one_sp_04_dp.snp.m.depthmask_merged.bed"

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# loading modules needed
module load bcftools

echo "Modules loaded" | ts '[%Y-%m-%d %H:%M:%S]'

# Make folder for output if it does not exist
mkdir -p $output_dir

echo "Folder for results prepared" | ts '[%Y-%m-%d %H:%M:%S]'

cd $SCRATCHDIR



### filter vcf
# --threads 2 adds two threads for compression (so 3 CPUs are needed)

## Filtering step 05 - apply depthmask
bcftools view --threads 2 -T ^$depthmask -o ${output_dir}${vcf_base}_05_dp.snp.m.dm.vcf.gz ${output_dir}${vcf_base}_04_dp.snp.m.vcf.gz || 
	{ echo >&2 "Error in bcftools!"; exit 2; }
	
tabix ${output_dir}${vcf_base}_05_dp.snp.m.dm.vcf.gz

echo "Filtering step 05 done" | ts '[%Y-%m-%d %H:%M:%S]'

echo "Values in FILTER column of the resulting vcf:"
# Check the values in FILTER column
bcftools query -f '%FILTER\n' ${output_dir}${vcf_base}_05_dp.snp.m.dm.vcf.gz | sort | uniq -c | sort -nr

echo "Filtering completed" | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch

# Resources: 5 h 31 min, 64 % CPU, 100 % memory
```

```bash
cd /path/to/arenosa_lyrata_serpentine/fastq_to_vcf/metacentrum_scripts/
qsub lyrataarenosa_one_species_2_filtering_vcf_05_apply_depthmask.sh
```

#### Output of the script

```
[2026-06-25 16:24:30] Filtering step 05 done
Values in FILTER column of the resulting vcf:
16705937 PASS
3966439 MultiAllelic
```


