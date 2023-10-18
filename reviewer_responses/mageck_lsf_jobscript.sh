#BSUB -q normal
#BSUB -J HAP1_MAGeCK
#BSUB -R "select[mem>5000] rusage[mem=5000] span[hosts=1]"
#BSUB -M 5000
#BSUB -oo /lustre/scratch124/casm/team113/users/vo1/hap1_wt_vs_lig4_daSilva_crispr_screen/MAGeCK/mageck_stdout.o
#BSUB -eo /lustre/scratch124/casm/team113/users/vo1/hap1_wt_vs_lig4_daSilva_crispr_screen/MAGeCK/mageck_stderr.e

# Load MAGeCK module which uses Singularity container
module load /software/team113/modules/modulefiles/mageck/0.5.9.3

# Run MAGeCK comparing HAP1 WT to Plasmid
mageck test \
	--norm-method 'none' \
	--remove-zero 'none' \
	-t 'WT_rep1,WT_rep2,WT_rep3' \
	-c 'plasmid' \
	-k '/lustre/scratch124/casm/team113/users/vo1/hap1_wt_vs_lig4_daSilva_crispr_screen/daSilva_normalised_counts_for_mageck.tsv' \
	-n '/lustre/scratch124/casm/team113/users/vo1/hap1_wt_vs_lig4_daSilva_crispr_screen/MAGeCK/WT_vs_plasmid'

# Run MAGeCK comparing HAP1 LIG4 to Plasmid
mageck test \
	--norm-method 'none' \
	--remove-zero 'none' \
	-t 'LIG4_rep1,LIG4_rep2,LIG4_rep3' \
	-c 'plasmid' \
	-k '/lustre/scratch124/casm/team113/users/vo1/hap1_wt_vs_lig4_daSilva_crispr_screen/daSilva_normalised_counts_for_mageck.tsv' \
	-n '/lustre/scratch124/casm/team113/users/vo1/hap1_wt_vs_lig4_daSilva_crispr_screen/MAGeCK/LIG4_vs_plasmid'

# Run MAGeCK comparing HAP1 LIG4 to HAP1 WT
mageck test \
	--norm-method 'none' \
	--remove-zero 'none' \
	-t 'LIG4_rep1,LIG4_rep2,LIG4_rep3' \
	-c 'WT_rep1,WT_rep2,WT_rep3' \
	-k '/lustre/scratch124/casm/team113/users/vo1/hap1_wt_vs_lig4_daSilva_crispr_screen/daSilva_normalised_counts_for_mageck.tsv' \
	-n '/lustre/scratch124/casm/team113/users/vo1/hap1_wt_vs_lig4_daSilva_crispr_screen/MAGeCK/LIG4_vs_WT'

