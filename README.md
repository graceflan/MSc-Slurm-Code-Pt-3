# MSc-Slurm-Code-Pt-3
Part 3: IQ Tree to Astral   
Note: Remove genes found in overlost.txt again. 
You will need IQTREE, raxml-ng and Astral.

## IQ Tree
Quick advice note: put *name*_alM_r_cT_f.fasta in a new file with your gene name file.

```
#!/bin/bash
#
#SBATCH --chdir=/home/DIR/HP_out/alignments/alM_r_o/IQtree
#SBATCH --job-name=iqtree
#SBATCH --partition=short       
#SBATCH --array=1-*number of genes*%25
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=email@kew.org
#SBATCH --mail-type=END,FAIL

echo $SLURM_ARRAY_TASK_ID

name=$(awk -v lineid=$SLURM_ARRAY_TASK_ID 'NR==lineid{print;exit}' /home/gflanaga/DIR/HP_out/alignments/edited/alM_r_o/IQtree/test-gene-names-new.txt)

echo $name

/home/DIR/apps/iqtree-1.6.12-Linux/bin/iqtree -s "$name".fasta -m MF -AICc -mset JC69,HKY85,GTR,K80 -nt AUTO -ntmax 2
```

Upload IQTREE
```

scp iqtree-script.sh DIR/HP_out/alignments/edited/alM_r_o/IQtree
```

Run IQTREE
```
sbatch /home/DIR/HP_out/alignments/edited/alM_r_o/IQtree/iqtree-script.sh
```














