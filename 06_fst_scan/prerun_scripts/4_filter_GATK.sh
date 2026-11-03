set +eu
source /path/to/SW/mambaforge/bin/activate gatk4
set -eu
GATK4=/path/to/SW/downloaded_gatk4/gatk-4.5.0.0/gatk
module load picard-2.8.1
module load bedops-2.4.39
module load openjdk