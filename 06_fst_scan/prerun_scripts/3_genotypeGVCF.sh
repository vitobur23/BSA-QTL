# module add picard-2.8.1
set +eu
source /path/to/SW/mambaforge/bin/activate gatk4
set -eu
GATK4=/path/to/SW/downloaded_gatk4/gatk-4.5.0.0/gatk
# module load bcftools

# It seems that the Java version in Conda is not working properly (at least sometimes), it is needed to load the Java module (2026-01-27)
module load openjdk
