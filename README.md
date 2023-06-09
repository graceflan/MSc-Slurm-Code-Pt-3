# MSc-Slurm-Code-Pt-3
Part 3: IQ Tree to Astral   
Note: Remove genes found in overlost.txt again. 
You will need IQTREE, raxml-ng and Astral.

## IQ Tree
Quick advice note: put *name*_alM_r_cT_f.fasta in a new file with your gene name file.

```
#!/bin/bash
#
#SBATCH --chdir=/home/DIR/HP_out/alignments/edited/alM_r_o/IQtree
#SBATCH --job-name=iqtree
#SBATCH --partition=short       
#SBATCH --array=1-*number of genes*
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=email@kew.org
#SBATCH --mail-type=END,FAIL

echo $SLURM_ARRAY_TASK_ID

name=$(awk -v lineid=$SLURM_ARRAY_TASK_ID 'NR==lineid{print;exit}' /home/DIR/HP_out/alignments/edited/alM_r_o/IQtree/test-gene-names-new.txt)

echo $name

/home/DIR/apps/iqtree-1.6.12-Linux/bin/iqtree -s "$name"_alM_r_cT_f.fasta -m MF -AICc -nt AUTO -ntmax 2

```

Upload IQTREE
```

scp iqtree-script.sh DIR/HP_out/alignments/edited/alM_r_o/IQtree
```

Run IQTREE
```
sbatch /home/DIR/HP_out/alignments/edited/alM_r_o/IQtree/iqtree-script.sh
```

## RAxML

### Get the best model for each gene from the output of IQtree
Extract the best model from all IQtree output files, from the folder containing the IQtree outputs, run (in an interractive slurm window, not in a slurm script):
```
ls *.log > file_names.txt
while read f; do grep "Best-fit model" $f; done < file_names.txt >> models.txt
paste -d "\t" file_names.txt models.txt > All_models.txt
rm file_names.txt
rm models.txt
```
Edit the model list to keep only the alignment file name and the model
```
sed -e 's/.fasta.log//g' -e 's/Best-fit model: //g' -e 's/ chosen according to AICc//g' All_models.txt > All_models_2.txt
```

Check what models were found and if some names need to be edited to be recognised by RAxML-NG
```
cut -f2 All_models_2.txt > All_models_2_models.txt
```
```
sort -u All_models_2_models.txt > All_models_2_models_su.txt
```
Download the All_models_2_models_su.txt file and open it in excel to see what models there are
Look if the model name needs to be edited and write down what should be the new model name  (see Sidonie's spreadsheet + here https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model and here: https://github.com/ddarriba/modeltest/wiki/Models-of-Evolution)

Run an edited version of the command below to edit the model names in the original list (in an interactive slurm window):
```
sed -e 's/GTR+/GTR+/g' -e 's/HKY+/HKY+/g' -e 's/K2P+/K80+/g' -e 's/K3P+/K81+/g' -e 's/K3Pu+/K81uf+/g' -e 's/SYM+/SYM+/g' -e 's/TIM2e+/TIM2+/g' -e 's/TIM2+/TIM2uf+/g' -e 's/TIM3e+/TIM3+/g' -e 's/TIM3+/TIM3uf+/g' -e 's/TIMe+/TIM1+/g' -e 's/TIM+/TIM1uf+/g' -e 's/TNe+/TN93ef+/g' -e 's/TN+/TN93+/g' -e 's/TPM2+/TPM2+/g' -e 's/TPM2u+/TPM2uf+/g' -e 's/TPM3+/TPM3+/g' -e 's/TPM3u+/TPM3uf+/g' -e 's/TVMe+/TVMef+/g' -e 's/TVM+/TVM+/g' All_models_2.txt > All_models_3.txt
```

You may need to edit that *_3.txt list locally by hand for the models such as K2P and K3P that do not have a "+" or so after the model name, because making a sed command for them like above may lead to unexpected results. Save the edited file as All_models_4.txt and reupload and make sure that end of lines are still just \n


Generate a list of models with their edited names, and a list of corresponding gene names (in an interactive slurm window):
```
cut -f2 All_models_4.txt > All_models_4_models.txt
cut -f1 All_models_4.txt > All_models_4_names.txt
```
The first model in All_models_4_models.txt is the model selected for the first gene in All_models_4_names.txt and so on...

### Run raxml-ng (because standard bootstrap is faster than in IQtree and more reliable or interpretable than their ultra-fast bootstrap, and more accurate and quick than RAxML)

```
#!/bin/bash
#
#SBATCH --chdir=/home/DIR/HP_out/alignments/edited/alM_r_o/IQtree
#SBATCH --job-name=raxml
#SBATCH --partition=short
#SBATCH --array=1-265
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --mail-user=email@kew.org
#SBATCH --mail-type=END,FAIL

echo $SLURM_ARRAY_TASK_ID

name=$(awk -v lineid=$SLURM_ARRAY_TASK_ID 'NR==lineid{print;exit}' /home/DIR/HP_out/alignments/edited/alM_r_o/IQtree/All_models_4_names.txt)

echo $name

model=$(awk -v lineid=$SLURM_ARRAY_TASK_ID 'NR==lineid{print;exit}' /home/DIR/HP_out/alignments/edited/alM_r_o/IQtree/All_models_4_models.txt)

echo $model

rm "$name".raxml.*
rm RF6_$name*

# raxmlHPC-PTHREADS -T 8 -m GTRGAMMA -f a -p 2345 -x 2345 -# 500 -k -s "$name".fasta -n "$name"_tree
#mv RAxML_*"$name"_tree trees

raxml-ng --all --msa "$name".fasta --model "$model" --bs-trees 500 --threads auto{8} --tree pars{30},rand{30} --seed 1 --prefix $name
# Check that only one tree (global optimum) instead of many (local optima)
raxml-ng --rfdist --tree $name.raxml.mlTrees --prefix RF6_$name
# Check convergence of BP trees
raxml-ng --bsconverge --bs-trees $name.raxml.bootstraps --prefix $name --seed 1 --threads auto{8} --bs-cutoff 0.03
```

Upload
```
scp RAxML-script.sh DIR/HP_out/alignments/edited/alM_r_o/IQtree
```

Run
```
sbatch /home/DIR/HP_out/alignments/edited/alM_r_o/IQtree/RAxML-script.sh
```

### SPECIES TREE
I recommend making a directory in HP_out called CDS and moving all your raxml files there to keep the data from your final steps together.	
Create a file called Tree_files_names_for_Astral.txt and write "all_trees" in it	
Check that the names in your gene name list correspond to the genes with which to do the species tree (so exclude failed raxml genes from this list for instance)

## ASTRAL
Download ASTRAL, Phyutility, and Newick.

```
#!/bin/bash
#
#SBATCH --chdir=/home/DIR/HP_out/CDS
#SBATCH --job-name=astral
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G 
#SBATCH --mail-user=email@kew.org
#SBATCH --mail-type=END,FAIL


#cd /home/DIR/HP_out/CDS
#while read name
#do cat "$name"_alM_r_cT_f.raxml.support >> /home/DIR/HP_out/CDS/all_trees.tre && echo "" >> /home/DIR/HP_out/CDS/all_trees.tre
#done < /home/DIR/HP_out/CDS/test-gene-names-new.txt


cd /home/DIR/HP_out/CDS
while read name;
do	
	/home/DIR/apps/conda/bin/nw_ed "$name".tre 'i & b<=10' o > "$name"_BP10.tre
	java -jar /home/DIR/apps/ASTRAL-master/Astral/astral.5.7.8.jar -i all_trees_BP10.tre -t 2 -o "$name"_BP10_SpeciesTree_annotQ.trel
	java -jar /home/DIR/apps/ASTRAL-master/Astral/astral.5.7.8.jar -i all_trees_BP10.tre -t 0 -o "$name"_BP10_SpeciesTree.tre
	/home/DIR/apps/phyx-master/src/pxrr -t all_trees_BP10_SpeciesTree.tre -g *outgroup* > all_trees_BP10_SpeciesTree_rooted.tre
	/home/DIR/apps/phyx-master/src/pxrr -t all_trees_BP10_SpeciesTree_annotQ.tre -g *outgroup* -o all_trees_BP10_SpeciesTree_annotQ_rooted.tre
	sed 's/\;\n/\;\r\n/' "$name"_BP10_SpeciesTree_rooted.tre > all_trees_BP10_SpeciesTree_rooted2.tre
	sed 's/\;\n/\;\r\n/' "$name"_BP10_SpeciesTree_annotQ_rooted.tre > all_trees_BP10_SpeciesTree_annotQ_rooted2.tre
	done < /home/DIR/HP_out/CDS/Tree_files_names_for_Astral.txt
```
If it only makes all_trees.tre and doesn't continue to the second section check:	
- If there are extra spaces between lines in all_trees.tre	
- All the directories are correct	
- If it is still not working replace the '$name' sections with all_trees or what you picked	

Download the species trees generated by ASTRAL	
	
YAY ALL DONE :)		
Well at least for slurm. Now move to R to turn the tre files into a phylogenetic tree.
