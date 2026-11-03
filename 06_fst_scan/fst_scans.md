# Serpentine adaptation - A. arenosa and A. lyrata (mapped to arenosa) - Fst scans

*Milos Duchoslav, 2026-03*

**Contents**

- [Preparations](#preparations)
  - [New folder](#new-folder)
  - [Preparation of bed file with windows for Piawka](#preparation-of-bed-file-with-windows-for-piawka)
  - [Group file for Piawka](#group-file-for-piawka)
  - [Lists of samples for each population pair](#lists-of-samples-for-each-population-pair)
- [Subsetting vcf](#subsetting-vcf)
  - [Script for subsetting vcf](#script-for-subsetting-vcf)
- [Running Piawka on pairs of populations](#running-piawka-on-pairs-of-populations)
  - [Script for running Piawka on pairs](#script-for-running-piawka-on-pairs)
- [Analysis in R](#analysis-in-r)
  - [Preparation of data](#preparation-of-data)
  - [Filter windows and assign genes](#filter-windows-and-assign-genes)



## Preparations

### New folder
```bash
mkdir -p /path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input
```

### Preparation of bed file with windows for Piawka

Use `piawka win` and VCF header to get bed file with 1kb windows as input for Piawka.

```bash
cd /path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input

# load bcftools
module load bcftools
# add piawka to path
export PATH="/path/to/SW/piawka:${PATH}"

# vcf (info about scaffold sizes in header will be used to make windows)
vcf="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/lyrataarenosa.one_sp_05_dp.snp.m.dm.vcf.gz"
# make windows
piawka win -v $vcf -w 1000 | head
piawka win -v $vcf -w 1000 > windows_1kb.bed

# filter the bed file for the testing
head -n 50 windows_1kb.bed > windows_1kb_fragment.bed
```

### Group file for Piawka

```bash
vcf="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/lyrataarenosa.one_sp_05_dp.snp.m.dm.vcf.gz"

cd /path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input

module load bcftools

# list samples
bcftools query -l $vcf > sample_list_one_species_2.txt

# add pops as groups to second column
paste sample_list_one_species_2.txt <(sed -E 's/-.*$//' sample_list_one_species_2.txt) > groups_2_pre_merging.txt

# merge some lyrata Norwegian pops
# exclude FEL-10 as it is probably mixture of two samples from different populations
# remove PEK-02, which is an outlier
cat groups_2_pre_merging.txt | 
	sed 's/\tNKM$/\tSNK/' | 
	sed 's/\tSUN$/\tSNK/' | 
	sed 's/\tSPN$/\tLSP/' | 
	sed 's/\tLOM$/\tLSP/' | 
	sed 's/^FEL-10ta\tFEL$/FEL-10ta\texclude/' | 
	sed 's/^PEK-02tl\tPEK$/PEK-02tl\texclude/' | 
	cat > groups_2.txt

```

### Lists of samples for each population pair

```bash
## Make a list of pairs just for lyrata and arenosa
cd /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/

head -n9 pop_pairs_1.tsv > pop_pairs_1_arenosa.tsv
tail -n4 pop_pairs_1.tsv > pop_pairs_1_lyrata.tsv


## prepare lists of samples for each pair

cd /path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input

sample_groups="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input/groups_2.txt"

while IFS=$'\t' read -r pop1 pop2 || [[ -n "$pop1" || -n "$pop2" ]]; do
	echo "$pop1 + $pop2"
	grep -P "\t$pop1|\t$pop2" $sample_groups | cut -f1 > samples_list_${pop1}_${pop2}.txt
done  < /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv
# done  < <(head -n 1 /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv)
```



#### Testing version 3 - devide vcf to pairs

Extract just pop pairs from vcf and filter out invariant sites.

```bash
cd /path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/tests

# load bcftools module to get all the dependencies
module load bcftools
# add piawka to path
export PATH="/path/to/SW/piawka:${PATH}"

cd /path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/tests

# get just a piece of vcf for testing
vcf="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/lyrataarenosa.one_sp_05_dp.snp.m.dm.vcf.gz"
bcftools view -r 'scaffold_1:0-50000' -o lyrataarenosa.one_sp_05_dp.snp.m.dm_50kb.vcf.gz $vcf
tabix lyrataarenosa.one_sp_05_dp.snp.m.dm_50kb.vcf.gz


## groups for samples
cd /path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input
less groups_2.txt
	
	
# Make a list of pairs just for lyrata and arenosa
cd /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/

head -n9 pop_pairs_1.tsv > pop_pairs_1_arenosa.tsv
tail -n4 pop_pairs_1.tsv > pop_pairs_1_lyrata.tsv

cd /path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input



## prepare lists of samples for each pair
sample_groups="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input/groups_2.txt"

while IFS=$'\t' read -r pop1 pop2 || [[ -n "$pop1" || -n "$pop2" ]]; do
	echo "$pop1 + $pop2"
	grep -P "\t$pop1|\t$pop2" $sample_groups | cut -f1 > samples_list_${pop1}_${pop2}.txt
done  < /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv
# done  < <(head -n 1 /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv)


## extract pair from vcf
cd /path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/tests

module load bcftools

input_vcf="lyrataarenosa.one_sp_05_dp.snp.m.dm_50kb.vcf.gz"
base_input_vcf=$(basename $input_vcf .vcf.gz)

while IFS=$'\t' read -r pop1 pop2 || [[ -n "$pop1" || -n "$pop2" ]]; do
	echo "$pop1 + $pop2"
	samples_list="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input/samples_list_${pop1}_${pop2}.txt"
	# subset the vcf
	bcftools view -S $samples_list -O u --threads 2 $input_vcf | # filter samples
	bcftools +fill-AN-AC -Ou | # Recalculate AN and AC
	bcftools view -a -Ou | # Remove alleles that are not in the genotypes
	bcftools view -m2 --threads 2 -o ${base_input_vcf}_${pop1}_${pop2}.vcf.gz # Remove invariants
	# make index
	bcftools index --tbi ${base_input_vcf}_${pop1}_${pop2}.vcf.gz
done  < /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv
# done  < <(head -n 1 /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv)


## calculate Fst
# load bcftools module to get all the dependencies
module load bcftools
# add piawka to path
export PATH="/path/to/SW/piawka:${PATH}"

sample_groups="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input/groups_2.txt"
bed="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input/windows_1kb_fragment.bed"
miss=0.5 # maximum share of missing genotypes per group for filtering by Piawka
lines=10 # minimum SNPs per window for filtering
jobs=1

output="Fst_scan_test_5_persite_filt"
pop1=VLA
pop2=BOR


while IFS=$'\t' read -r pop1 pop2 || [[ -n "$pop1" || -n "$pop2" ]]; do
	echo "$pop1 + $pop2"
	vcf="${base_input_vcf}_${pop1}_${pop2}.vcf.gz"
	
output="Fst_scan_test_5_persite_filt_sum_filt"
time piawka calc --persite --jobs $jobs --vcf $vcf --groups <(grep -P "\t$pop1|\t$pop2" $sample_groups) --bed $bed --stats fst,fstwc,afd,miss,lines,segr,maf,dxy --mult --dependencies | piawka filt -e "miss <= $miss" | piawka sum --bed $bed | piawka filt -e "lines >= $lines" | bgzip > ${output}_${pop1}_${pop2}.bed.gz
zgrep -c 'fstwc' ${output}_${pop1}_${pop2}.bed.gz

# done  < /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv
done  < <(head -n 1 /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv)



```
## Subsetting vcf

Separate vcf on population pairs and remove sites that are invariant in this subset. That improves calculation of number of SNPs per window.

When Piawka will be improved to produce also something like "between-groups lines" or "between-groups segr", this will not be needed anymore.

### Script for subsetting vcf
```bash
### Script for Metacentrum
## Subsetting of VCF for scans

#PBS -N lyrataarenosa_subsetting_vcf
#PBS -l select=1:ncpus=3:mem=64gb:scratch_local=1000gb:brno=True
#PBS -l walltime=24:00:00 
#PBS -m ae

input_vcf="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/lyrataarenosa.one_sp_05_dp.snp.m.dm.vcf.gz"
output_dir="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input"
# other needed inputs:
# - samples lists ("/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input/samples_list_${pop1}_${pop2}.txt")
# - list of population pairs (/path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv)

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# loading modules needed
module load bcftools

echo "Modules loaded" | ts '[%Y-%m-%d %H:%M:%S]'

# Make folder for output if it does not exist
mkdir -p $output_dir

echo "Folder for results prepared" | ts '[%Y-%m-%d %H:%M:%S]'

cd $output_dir


# VCF filename without path and suffix
base_input_vcf=$(basename $input_vcf .vcf.gz)

### filter vcf
# --threads 2 adds two threads for compression (so 3 CPUs are needed)

echo "Start filtering" | ts '[%Y-%m-%d %H:%M:%S]'

while IFS=$'\t' read -r pop1 pop2 || [[ -n "$pop1" || -n "$pop2" ]]; do
	echo "Processing $pop1 + $pop2" | ts '[%Y-%m-%d %H:%M:%S]'
	samples_list="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input/samples_list_${pop1}_${pop2}.txt"
	# subset the vcf
	bcftools view -S $samples_list -O u --threads 2 $input_vcf | # filter samples
		bcftools +fill-AN-AC -Ou | # Recalculate AN and AC
		bcftools view -a -Ou | # Remove alleles that are not in the genotypes
		bcftools view -m2 --threads 2 -o ${base_input_vcf}_${pop1}_${pop2}.vcf.gz # Remove invariants
	# make index
	bcftools index --tbi ${base_input_vcf}_${pop1}_${pop2}.vcf.gz
done  < /path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv

echo "Filtering done" | ts '[%Y-%m-%d %H:%M:%S]'


# clean the SCRATCH directory
clean_scratch

# Resources: 3,5 h, 57 % CPU, 48 GB memory
```

```bash
cd /path/to/arenosa_lyrata_serpentine/Fst_scans/metacentrum_scripts/
qsub lyrataarenosa_subsetting_vcf.sh
```

## Running Piawka on pairs of populations

### Script for running Piawka on pairs

```bash
### Script for Metacentrum

#PBS -N lyrataarenosa_piawka_Fst_scan_on_pairs
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=50gb
#PBS -l walltime=24:00:00 
#PBS -m ae

# define variables
# directory for output
output_dir="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_output/"

# original vcf (before subsetting)
input_vcf="/path/to/arenosa_lyrata_serpentine/fastq_to_vcf/lyrataarenosa_results/vcf_one_species_2/lyrataarenosa.one_sp_05_dp.snp.m.dm.vcf.gz"
# directory with subsetted vcf
vcf_dir="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input/"
sample_groups="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input/groups_2.txt"
species_pairs="/path/to/arenosa_lyrata_serpentine/basic_popgen/metadata/pop_pairs_1_arenosa.tsv"
bed="/path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_input/windows_1kb.bed"
output="Fst_scan_02"


miss=0.5 # maximum share of missing genotypes per group for filtering by Piawka
lines=10 # minimum SNPs per window for filtering
jobs=8

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $PBS_O_WORKDIR/jobs_info.txt

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# redirect /tmp to scratch
export TMPDIR=$SCRATCHDIR

echo "Starting"  | ts '[%Y-%m-%d %H:%M:%S]'

vcf_base=$(basename $input_vcf .vcf.gz)

## loading bcftools and piawka
# load bcftools module to get all the dependencies
module load bcftools
# add bcftools to path
# export PATH=/path/to/SW/bcftools-1.21/bin:$PATH
echo ""
echo "bcftools version:"  | ts '[%Y-%m-%d %H:%M:%S]'
bcftools --version
# add piawka to path
export PATH="/path/to/SW/piawka:${PATH}"
# check piawka version
echo ""
echo "Piawka version:"  | ts '[%Y-%m-%d %H:%M:%S]'
piawka
echo ""
echo "Piawka commit:"  | ts '[%Y-%m-%d %H:%M:%S]'
git -C /path/to/SW/piawka rev-parse HEAD
echo ""

echo "SW loaded."  | ts '[%Y-%m-%d %H:%M:%S]'

cd $output_dir

## calculate Fst
echo "Starting calculation."  | ts '[%Y-%m-%d %H:%M:%S]'

# loop over population pairs
while IFS=$'\t' read -r pop1 pop2 || [[ -n "$pop1" || -n "$pop2" ]]; do
	echo "Processing $pop1 + $pop2" | ts '[%Y-%m-%d %H:%M:%S]'
	# vcf for one pair
	vcf="${vcf_dir}${vcf_base}_${pop1}_${pop2}.vcf.gz"
	# piawka calculation
	piawka calc --persite --jobs $jobs --vcf $vcf --groups <(grep -P "\t$pop1|\t$pop2" $sample_groups) --bed $bed \
						--stats fst,fstwc,afd,miss,lines,segr --mult --dependencies --quiet | # calculate values per site
		piawka filt -e "miss <= $miss" | # filter out individual sites that have too much missing values
		piawka sum --bed $bed | # sum the values per windows
		piawka filt -e "lines >= $lines" | # filter out windows with too few SNPs (this is slightly problematic step, see below)
		bgzip > ${output}_${pop1}_${pop2}.bed.gz  || { echo >&2 "Piawka failed!"; exit 3; } # compress
done  < $species_pairs
# done  < <(head -n 1 $species_pairs) # for testing

echo "Scans done."  | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch

# Resources: 
# With 16 CPUs: 1 h, 100 % (32 GB) memory, 94 % CPU
```
##### Filtration within piawka calculation

There are these steps:

1. Calculation of Fst and other variables per individual site (SNP).
2. Missingness filtration - sites with more than 50 % of missing data within population are removed. If the other population in pair has for this site enough data, the information of within-population variables (lines, segr, miss) stays there. The between-population variables are removed for such site.
3. Summing up - the values are summed-up for 1kb windows. The number of lines for pop2 includes also the lines where information for pop1 was removed due to high missingness and vice versa.
4. Number of lines filtration - between-population variables from windows with less than 10 lines in at least one population are removed. Due to behaviour described above, the real number of SNPs used to calculate Fst can be lower. This can happen in case that different lines were removed from each population due to missingness, as the lines in the other population were included in calculation of number of lines.

```bash
cd /path/to/arenosa_lyrata_serpentine/Fst_scans/metacentrum_scripts/
qsub lyrataarenosa_piawka_Fst_scan_on_pairs.sh
```

#### Checking results

```bash
cd /path/to/arenosa_lyrata_serpentine/Fst_scans/lyrataarenosa/piawka_output/


```

## Analysis in R

### Preparation of data

```r
setwd("/path/to/arenosa_lyrata_serpentine/")
old.par<-par(no.readonly = T)

library(data.table)
library(qqman) # Needed for manhattan plot
# library(VennDiagram)
library(UpSetR)
library(SuperExactTest)

### setting of species
species <- c("Arabidopsis arenosa", "AA")

## pairs and names
# read table with pairs
pop.pairs <- read.table(file = "basic_popgen/metadata/pop_pairs_1.tsv",
												 header = F, sep = "\t", quote = "", comment.char = "")
colnames(pop.pairs) <- c("N", "S")
# Name the pairs
pop.pairs$pair <- c(paste0("AA", 1:9), paste0("AL", 1:4))
pop.pairs$pair2 <- paste0(pop.pairs$pair, " (", pop.pairs$N, "-", pop.pairs$S, ")")

pop.pairs.species <- pop.pairs[grep(pattern = species[2], x = pop.pairs$pair), ]


# scan files
scans <- list.files(path = "Fst_scans/lyrataarenosa/piawka_output/")
scans.f <- scans[grep(pattern = "Fst_scan_02_.*\\.bed.gz", x = scans)]
scan.pairs <- gsub(pattern = "Fst_scan_02_|.bed.gz", replacement = "", x = scans.f)
scans.species <- scans.f[match(x = paste(pop.pairs.species$N, pop.pairs.species$S, sep = "_"), table = scan.pairs)]
scans.species.f <- scans.species[!is.na(scans.species)]

# prepare list for scans
pi.res.1 <- as.list(vector(length = length(scans.species.f)))
names(pi.res.1) <- pop.pairs.species$pair[!is.na(scans.species)]


i=1
i=4

# reading files
for(i in seq_along(scans.species.f)) {
	scan.file <- scans.species.f[i]
	# read data from scans
	# Note: Package 'R.utils' is needed to be installed (but not loaded) to read .gz file
	one.pair <- fread(file = paste0("Fst_scans/lyrataarenosa/piawka_output/", scan.file), header=T)
	# remove "#" from colnames
	colnames(one.pair) <- sub(pattern = "#", replacement = "", x = colnames(one.pair))
	# save in list
	pi.res.1[[i]] <- one.pair
	
	
}

# checks
one.pair <- pi.res.1[[1]]
one.pair[one.pair$stat == "lines" & one.pair$pop1 == "VLA", ]
dim(one.pair[one.pair$stat == "lines" & one.pair$pop1 == "VLA", ])
dim(one.pair[one.pair$stat == "lines" & one.pair$pop1 == "BOR", ])

table(one.pair$stat, one.pair$pop1)


### change to wide format
# also add center of the window

pi.res.2 <- as.list(vector(length = length(scans.species.f)))
names(pi.res.2) <- names(pi.res.1)

for(i in seq_along(scans.species.f)) {
	one.pair <- pi.res.1[[i]]
	# separate scaffolds and windows location to separate columns, add middle of the window
	# one.pair$scaff <- sub(pattern = ":.*$", replacement = "", x = one.pair$locus)
	# one.pair$wStart <- as.integer(gsub(pattern = "^.*:|-.*$", replacement = "", x = one.pair$locus))
	# one.pair$wEnd <- as.integer(gsub(pattern = "^.*-", replacement = "", x = one.pair$locus))
	
	
	### change to wide format
	w.one.pair.1 <- unique(one.pair[, 1:4])
	w.one.pair.1$wCenter <- as.integer(rowMeans(w.one.pair.1[, c("start", "end")]))
	w.one.pair.1$chrNum <- as.integer(gsub(pattern = "scaffold_", replacement = "", x = w.one.pair.1$chr))
	# variable = "fst"
	
	## add variables for pair
	if(exists("w.one.pair.2")) remove(w.one.pair.2) # needs to be removed if it exist for loop to work safely
	
	for(variable in c("fst", "fstwc", "afd")) {
		if(!exists("w.one.pair.2")) {
			w.one.pair.2 <- merge.data.table(x = w.one.pair.1, y = one.pair[one.pair$stat == variable, c(1:6, 8)], by = c("chr", "start", "end", "locus"), all.x = T)
		}
		else {
			w.one.pair.2 <- merge.data.table(x = w.one.pair.2, y = one.pair[one.pair$stat == variable, c(1:4, 8)], by = c("chr", "start", "end", "locus"), all.x = T)
		}
		colnames(w.one.pair.2)[colnames(w.one.pair.2) == "value"] <- variable
	}
	
	# get pop names
	pop.1 <- w.one.pair.2$pop1[!is.na(w.one.pair.2$pop1)][1]
	pop.2 <- w.one.pair.2$pop2[!is.na(w.one.pair.2$pop1)][1]
	
	## add variables for single pops
	w.one.pair.3 <- w.one.pair.2
	for(variable in c("lines", "miss", "segr")) {
		w.one.pair.3 <- merge.data.table(x = w.one.pair.3, y = one.pair[one.pair$stat == variable & one.pair$pop1 == pop.1, c(1:4, 8)], 
																		 by = c("chr", "start", "end", "locus"), all.x = T)
		colnames(w.one.pair.3)[colnames(w.one.pair.3) == "value"] <- paste0(variable, ".pop1")
		w.one.pair.3 <- merge.data.table(x = w.one.pair.3, y = one.pair[one.pair$stat == variable & one.pair$pop1 == pop.2, c(1:4, 8)], 
																		 by = c("chr", "start", "end", "locus"), all.x = T)
		colnames(w.one.pair.3)[colnames(w.one.pair.3) == "value"] <- paste0(variable, ".pop2")
	}
	
	
	# save in the list
	pi.res.2[[i]] <- w.one.pair.3
}

# checks
one.pair <- pi.res.2[[1]]
one.pair <- pi.res.2[[4]]

plot(one.pair$fst, one.pair$fstwc)

hist(one.pair$lines.pop1)
hist(one.pair$lines.pop2)

plot(one.pair$lines.pop1, one.pair$lines.pop2)

hist(one.pair$miss.pop1)
hist(one.pair$miss.pop2)

plot(one.pair$miss.pop1, one.pair$miss.pop2)

plot(one.pair$miss.pop1, one.pair$fst)
plot(one.pair$miss.pop1, one.pair$fstwc)

plot(rowMeans(one.pair[, c("miss.pop1", "miss.pop2")]), one.pair$fst)
plot(rowMeans(one.pair[, c("miss.pop1", "miss.pop2")]), one.pair$fstwc)

plot(rowMeans(one.pair[, c("lines.pop1", "lines.pop2")]), one.pair$fst)
plot(rowMeans(one.pair[, c("lines.pop1", "lines.pop2")]), one.pair$fstwc)

plot(one.pair$afd, one.pair$fst)
plot(one.pair$afd, one.pair$fstwc)

plot(one.pair$segr.pop1, one.pair$segr.pop2)
```

#### Preparation genome of annotation

```r
## read the gff file with annotation
gff.1 <- read.table(file = "/path/to/resources/arenosa_Paul_Grini/Arabidopsis_arenosa_genome.annotation.gtf", 
							 header = F, sep = "\t", quote = "" ,
							 comment.char = "#",
							 #, nrows = 1500
)
# gff.1 <- read.table(file = "/path/to/resources/ncbi_lyrata_reference_101/GCF_000004255.2_v.1.0_genomic.gff", 
# 										header = F, sep = "\t", comment.char = "#", quote = "" 
#                      #, nrows = 1500
#                     )
dim(gff.1)
summary(gff.1)
head(gff.1)
gff.1[1:50, ]

# Types of features
levels(as.factor(gff.1$V3))
table(as.factor(gff.1$V3))

# scaffold names
levels(as.factor(gff.1$V1))

# filter just for genes and pseudogenes
genes.gff <- gff.1[grep(pattern = "gene", x = gff.1$V3), ]
dim(genes.gff)
table(as.factor(genes.gff$V3))

# gene IDs to new column
genes.gff$gene_ID <- genes.gff$V9
# genes.gff$gene_ID <- gsub(pattern = "ID=gene-|;.*$", replacement = "", x = genes.gff$V9)

# change to data.table
genes.gff <- as.data.table(genes.gff)


# prepare annotation for fast overlap (it will also sort the table)
setkey(x = genes.gff, V1, V4, V5)



### add surroundings of genes (probable regulatory regions)
# check, if all the genes have direction
table(genes.gff[, 7])

## upstream region
# copy gff
upstr.gff <- genes.gff

# Set the start of the upstream region of gene:
# - either to 2 kb before the start of the gene 
# - or to the end of the previous gene.
# - If there are some overlapping genes, set the region to 2 kb.
upstr.gff[, 4] <- mapply(# function to be used in mapply
													FUN = function(chr_this, chr_previous, start_of_this, stop_of_previous) {
																	if(chr_this == chr_previous & start_of_this > stop_of_previous) { # test if the previous gene lies on the same chromosome 
																																																		# and if the genes are not overlapping
																			max((start_of_this - 2000), (stop_of_previous + 1))
																		} else
																			max((start_of_this - 2000), (0 + 1))
																	}, 
												 # parameters
												 chr_this = unlist(genes.gff[, 1]), # chromosome of the current gene
												 chr_previous = c("none", unlist(genes.gff[-nrow(genes.gff), 1])), # chromosome of the previous gene (the vector is shifted by one position)
												 start_of_this = unlist(genes.gff[, 4]), # start of the current gene
												 stop_of_previous = c(0, unlist(genes.gff[-nrow(genes.gff), 5])) # stop of the previous gene (again, the vector is shifted by one position)
												 )
# end of upstream region
upstr.gff[, 5] <- genes.gff[, 4] - 1 # end of the upstream region is start of the gene minus 1
# is it really upstream or downstream of the gene? (set it based on the direction)
upstr.gff[, 3] <- ifelse(test = genes.gff[, 7] == "+", yes = "upstream_region", no = "downstream_region")

## downstream region
# copy gff
downstr.gff <- genes.gff

# Set the end of the downstream region of gene:
# - either to 2 kb after the end of the gene 
# - or to the start of the next gene.
# - If there are some overlapping genes, set the region to 2 kb.
downstr.gff[, 5] <- mapply(# function to be used in mapply
													FUN = function(chr_this, chr_next, end_of_this, start_of_next) {
																	if(chr_this == chr_next & end_of_this < start_of_next) { # test if the next gene lies on the same chromosome 
																																																		# and if the genes are not overlapping
																			min((end_of_this + 2000), (start_of_next - 1))
																		} else
																			(end_of_this + 2000)
																	}, 
												 # parameters
												 chr_this = unlist(genes.gff[, 1]), # chromosome of the current gene
												 chr_next = c(unlist(genes.gff[-1, 1]), "none"), # chromosome of the next gene (the vector is shifted by one position)
												 end_of_this = unlist(genes.gff[, 5]), # end of the current gene
												 start_of_next = c(unlist(genes.gff[-1, 4]), 0) # start of the next gene (again, the vector is shifted by one position)
												 )
# start of downstream region
downstr.gff[, 4] <- genes.gff[, 5] + 1 # start of downstream region is end of the gene plus 1
# is it really upstream or downstream of the gene? (set it based on the direction)
downstr.gff[, 3] <- ifelse(test = genes.gff[, 7] == "+", yes = "downstream_region", no = "upstream_region")

## merge the gff with gene body, upstream and downstream region
genes.gff <- rbind(genes.gff, upstr.gff, downstr.gff)

# remove lines where the end is 0 (this happens in cases when gene begins exactly at the beginning of scaffold - at position 1)
genes.gff <- genes.gff[!genes.gff$V5 == 0, ]

# prepare annotation for fast overlap (it will also sort the table)
setkey(x = genes.gff, V1, V4, V5)

# checks
all(genes.gff[, 4] < genes.gff[, 5])
all(genes.gff[, 4] <= genes.gff[, 5])
genes.gff[as.logical(! genes.gff[, 4] < genes.gff[, 5]), ]
genes.gff[genes.gff$gene_ID %in% c("jg29852", "jg29853"), ]
genes.gff[genes.gff$gene_ID %in% c("jg15560", "jg15559"), ]
genes.gff[genes.gff$V1 %in% c("scaffold_105", "scaffold_106"), ]

# lenght of the upstream and downstream region
hist((genes.gff$V5 - genes.gff$V4)[genes.gff$V3 %in% c("upstream_region", "downstream_region")])
# number of regions with length of 2 kb
sum((genes.gff$V5 - genes.gff$V4)[genes.gff$V3 %in% c("upstream_region", "downstream_region")] == 1999)
# number of regions with length less than 2 kb
sum((genes.gff$V5 - genes.gff$V4)[genes.gff$V3 %in% c("upstream_region", "downstream_region")] < 1999)


## Functional annotation and thaliana orthologues

fun.annot <- readRDS(file = "/path/to/Brassicaceae_orthology/brassicaceae_3/functional_annotation/4_full_rds/Arabidopsis_arenosa_At_orthologues_and_functional_annotations_full.rds")
# fun.annot <- readRDS(file = "/path/to/resources/orthofinder/brassicaceae_2/R_analysis/Arabidopsis_lyrata_NCBI_At_orthologues_and_functional_annotations_full.rds")
dim(fun.annot)
fun.annot[1:10, 1:10]

# # Gene - protein translation table
# tr.gene.prot <- read.table(file = "/path/to/resources/ncbi_lyrata_reference_101/lyrata_NCBI_101_annotation_-_protein_id-mRNA_id-gene_id_translation_table.tsv", header = T, sep = "\t")
# head(tr.gene.prot)
# dim(tr.gene.prot)
# sort(table(tr.gene.prot$gene), decreasing = T)[1:30]
# 
# # add gene ID to fun.annot table
# fun.annot$Arabidopsis_lyrata_NCBI_gene <- tr.gene.prot$gene[match(x = fun.annot$Arabidopsis_lyrata_NCBI, table = tr.gene.prot$protein)]
# sort(table(fun.annot$Arabidopsis_lyrata_NCBI_gene), decreasing = T)[1:30]
# fun.annot[1:20, c("Arabidopsis_lyrata_NCBI", "Arabidopsis_lyrata_NCBI_gene")]
```

### Filter windows and assign genes

The results were already filtered for at least 10 lines and max. 0.5 missingness.

**Background sets of genes**

All genes/windows with number of SNPs over threshold.

Used later for calculation of super exact test and for calculation of fraction of genes and other features repeated and not repeated.

```r
# set threshold for lines of SNPs per window
lines_threshold <- 10
# lines_threshold <- 20
# lines_threshold <- 30
# lines_threshold <- 50

# set quantile for top windows
topq <- 0.95

# set threshold for missingness
miss_thr <- 0.5

# save Manhattan plots?
save.manhattan = F

## prepare empty lists
# selected candidates
full.candidates <- list()
candidate.proteins <- list()
# whole dataset (all windows/proteins with enough SNPs)
bckgr.full <- list()
bckgr.proteins <- list()

i=1
i=2

# make folder for results if it does not exist
dir.create("Fst_scans/lyrataarenosa/results/Fst_scan_02/figures", recursive = T)
dir.create("Fst_scans/lyrataarenosa/results/Fst_scan_02/gene_lists", recursive = T)

# loop to get candidates for each population pair
for(i in seq_along(pi.res.2)) {
	one.pair <- pi.res.2[[i]]
	
	if(save.manhattan) {
		# Manhattan plots with increasing treshold for number of SNPs
		pdf(file = paste0("Fst_scans/lyrataarenosa/results/Fst_scan_02/figures/manhattan_Fstwc_", 
											names(pi.res.2)[i], "_various_linesMins.pdf"), 
				width=26, height=14, pointsize=20, onefile = T)
		for(j in seq(from = 0, to = 150, by = 10)) {
			g.data.pre <- one.pair[one.pair$lines.pop1 >= j & one.pair$lines.pop2 >= j, ]
			g.data <- g.data.pre[g.data.pre$fstwc > 0.2, ] # show only the higher points (to not have all to low points in the plot)
			# function manhattan is from qqman package
			manhattan(x = g.data, chr = "chrNum", bp = "wCenter", p = "fstwc", snp = "wCenter", logp = F,
								ylab = paste0("Fst ", names(pi.res.2)[i]), 
								suggestiveline=as.numeric(quantile(g.data.pre$fstwc, 0.99, na.rm = T)),
								main = paste("Min. number of SNPs per window >=", j))
			mtext(paste("Number of windows:", nrow(g.data.pre)))
			mtext(paste("Blue line is showing 99% quantile."), col = "blue", line = -1)
		}
		dev.off()
	}
	
	# filter for minimal number of SNPs and missingness
	one.pair.f <- one.pair[one.pair$lines.pop1 >= lines_threshold & one.pair$lines.pop2 >= lines_threshold & 
												 	one.pair$miss.pop1 <= miss_thr & one.pair$miss.pop2 <= miss_thr, ]
	# dim(one.pair)
	# dim(one.pair.f)
	
	# calculate percentile
	one.pair.f$fstwc.percentile <- rank(one.pair.f$fstwc, ties.method= "min", na.last = F)/length(one.pair.f$fstwc) * 100
	
	## Assign genes to windows
	# annotate windows by genes
	all.windows <- foverlaps(x = one.pair.f, y = genes.gff, 
												 by.x = c("chr", "start", "end"), by.y = c("V1", "V4", "V5"), type = "any")
	

	# filter for 1% highest Fst (or other based on topq)
	top.windows <- all.windows[all.windows$fstwc.percentile >= topq * 100, ]
	dim(top.windows)
	dim(all.windows)
	dim(top.windows)/dim(all.windows) # the numbers does not fully correspond to percentile as there are some windows multiple times in "all.windows" and "top.windows" objects (covering gene, upstream or downstream region)
	
	
	# ## Assign genes to windows
	# # annotate windows by genes
	# top.genes <- foverlaps(x = top.windows, y = genes.gff, 
	# 											 by.x = c("chr", "start", "end"), by.y = c("V1", "V4", "V5"), type = "any")
	
	# # checks
	# table(as.factor(top.genes$V3))
	# # Number of windows without genes
	# sum(is.na(top.genes$V3))
	# # Numbers of top windows in one gene
	# sort(table(top.genes$gene_ID), decreasing = T)
	
	# Add functional annotation of genes
	all.fun <- merge(x = all.windows, y = fun.annot, by.x = "gene_ID", by.y = "Arabidopsis_arenosa", all.x = T)
	# top.fun <- merge(x = top.genes, y = fun.annot, by.x = "gene_ID", by.y = "Arabidopsis_lyrata_NCBI_gene", all.x = T)
	
	# Add functional annotation of genes to candidates
	top.fun <- merge(x = top.windows, y = fun.annot, by.x = "gene_ID", by.y = "Arabidopsis_arenosa", all.x = T)
	# top.fun <- merge(x = top.genes, y = fun.annot, by.x = "gene_ID", by.y = "Arabidopsis_lyrata_NCBI_gene", all.x = T)
	
	# sort by Fst
	all.fun <- all.fun[order(all.fun$fstwc, decreasing = T), ]
	top.fun <- top.fun[order(top.fun$fstwc, decreasing = T), ]
	
	# check the top genes
	# colnames(top.fun)
	# View(top.fun[1:30, c(1:2, 18:27, 27, 57:58, 68)])
	# View(top.fun[!is.na(top.fun$gene_ID)][1:30, c(1:2, 18:27, 27, 57:58, 68)])
	
	# export results
	write.table(x = top.fun[, -"InterProScan_pathways"], 
							file =  paste0("Fst_scans/lyrataarenosa/results/Fst_scan_02/gene_lists/candidates_", 
										names(pi.res.2)[i], "_linesMin_", lines_threshold, "_topq_", topq, "_miss_", miss_thr, ".tsv"), 
							sep = "\t", row.names = F, quote = F
							)
	
	# save results to list
	bckgr.full[[i]] <- all.fun
	names(bckgr.full)[i] <- names(pi.res.2)[i]
	full.candidates[[i]] <- top.fun
	names(full.candidates)[i] <- names(pi.res.2)[i]
	
	## extract just the gene/protein names
	all.proteins <- unique(all.fun$gene_ID[!is.na(all.fun$gene_ID)])
	#  save results to list
	bckgr.proteins[[i]] <- all.proteins
	names(bckgr.proteins)[i] <- names(pi.res.2)[i]
	
	top.proteins <- unique(top.fun$gene_ID[!is.na(top.fun$gene_ID)])
	# just.proteins <- unique(top.fun$Arabidopsis_lyrata_NCBI[!is.na(top.fun$Arabidopsis_lyrata_NCBI)])
	#  save results to list
	candidate.proteins[[i]] <- top.proteins
	names(candidate.proteins)[i] <- names(pi.res.2)[i]
	
	print(paste(names(pi.res.2)[i], "done"))
}

top.fun <- full.candidates$AA1

# export the candidate windows/genes of the first pair (non-serpentine VLA vs serpentine BOR)
	write.table(x = top.fun[, 1:27],
							file = "Fst_scans/lyrataarenosa/results/Fst_scan_02/gene_lists/fst_candidates.tsv",
							sep = "\t", row.names = F, quote = F
							)

# plot(rowMeans(top.fun[, c("miss.pop1", "miss.pop2")]), top.fun$fstwc, col = c(2, 1)[is.na(top.fun$Arabidopsis_thaliana) + 1])
# 
# plot(rowMeans(top.fun[, c("lines.pop1", "lines.pop2")]), top.fun$fstwc, col = c(2, 1)[is.na(top.fun$Arabidopsis_thaliana) + 1])
# 
# plot(top.fun$afd, top.fun$fstwc, , col = c(2, 1)[is.na(top.fun$Arabidopsis_thaliana) + 1])

# numbers of candidates
sapply(X = full.candidates, FUN = nrow)
sapply(X = candidate.proteins, FUN = length)

# numbers for background
sapply(X = bckgr.full, FUN = nrow)
sapply(X = bckgr.proteins, FUN = length)
```
