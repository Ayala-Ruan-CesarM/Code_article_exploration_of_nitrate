# First:  download reads with SRA-toolkit
```
IW="SRR14424469"
PW="SRR14424468"

woriking_dir="home" #Any other dir could be used
mkdir -p IW
mkdir -p PW
prefetch -O $woriking_dir/IW/ $IW
fasterq-dump -s -e 20 -O $woriking_dir/ -o IW $woriking_dir/$IW/*.sra

prefetch -O $woriking_dir/IW $IW
fasterq-dump -s -e 20 -O $woriking_dir/ -o PW $woriking_dir/$PW/*.sra

```
# Second: quality check with fastp
```
cd $woriking_dir/IW/
mkdir -p fastp
fastp -i IW_1.fastq -o fastp/IW_1_qc.fastq -I IW_2.fastq -O fastp/IW_2_qc.fastq -w 20 -h report.html
cd $woriking_dir/PW/
fastp -i PW_1.fastq -o fastp/PW_1_qc.fastq -I PW_2.fastq -O fastp/PW_2_qc.fastq -w 20 -h report.html
```
# Third: taxonomical assigment 
Here, we need the path to the .k2d index files needed for kraken.
```
DB="/mnt/d/DataBase/" # This is the path to download db_index
kraken2 --db "$DB" --threads --output PW_taxonomy.txt --report PW_taxonomy.report --paired fastp/PW_1_qc.fastq fastp/PW_2_qc.fastq
```
Now, read correction with bracken
```
bracken -d $DB" -i PW_taxonomy.report -o PW_taxonomy.bracken -w PW_taxonomy.bracken.report -r 150 -l G
```
The "PW_taxonomy.bracken" file will be use latter.
This same step should be applied to IW fastq reads in order to obtained the file "IW_taxonomy.bracke" that would be use latter

# Fourth: De-novo assembly and gene predciton 
## Reads merged 
```
mkdir -p Assembly
cat fastp/IW_1_qc.fastq fastp/PW_1_qc.fastq > Assembly/Merge_1_Reads.fastq ;
cat fastp/IW_2_qc.fastq fastp/PW_2_qc.fastq > Assembly/Merge_2_Reads.fastq  
```
## Spades assembly 
We're going to use merge reads from previous step
```
spades --meta -1 Assembly/Merge_1_Reads.fastq -2 Assembly/Merge_2_Reads.fastq -o Emsamble/ -t 30
```
## Prokka gene prediction
From Spades output dir "Ensamble" we're going to use the file "contigs.fasta" as input for prokka
Here we use --prefix to have all output files with the same basename and the same locus tag "OF" (Oil Field)
```
 metaprokka --outdir Prokka_annotation/ --force --prefix OILFIELD --locustag OF --mincontiglen 200 --compliant --cpus 30 Emsamble/contigs.fasta
```
# Fifth: gene extractions
On this step we're going to extract the genes related to Nitrate reduction and Sulfate reduction from all metagenome predicted genes
To do so for loop would be used:
```
#!/bin/bash 
gene_list=("narG" "narL" "narX" "narH" "narI" "narJ" "sat_" "dsrA" "dsrB") 
for gene in "${gene_list[@]}"; do
    grep "$gene" OILFILED.TSV | awk '{print $1}' | uniq > IDs_$gene.txt
done

# Unfourtnaly for aprA and aprB a change on the command was requiered
grep "aprA" OILFIELD.tsv | grep "Adenylylsulfate" | awk '{print $1}' | uniq > IDs_aprA.txt
grep "aprB" OILFIELD.tsv | grep "Adenylylsulfate" | awk '{print $1}' | uniq > IDs_aprB.txt

```

Now that we have all gene IDs that match the desired gene according to it "OF" locustag
we're going to extract the DNA sequence of each gene.

```
list=("IDs_narG.txt" "IDs_narL.txt" "IDs_narX.txt" "IDs_narH.txt" "IDs_narI.txt" "IDs_narJ.txt" "IDs_sat.txt" "IDs_aprA.txt" "IDs_aprB.txt" "IDs_dsrA.txt" "IDs_dsrB.txt")
for file in "${list[@]}"; do 
    basename=${filename#"IDs_"}
    basename=${basename%.txt}
    grep -f $file OILFIELD.ffn > "$basename"_genes.fasta
done
```
# Sixth: taxonomical annotation 
Then we're going to do the taxonomical assignment to each gene using kraken2
Two separete bash cycles are used to keep results separeted.

```
echo "Nitrate genes taxonomical assignment"
outdir_nitrate="Nitrate_receptors_annotations"
list=("narL_genes.fasta" "narX_genes.fasta" "narH_genes.fasta" "narI_genes.fasta" "narJ_genes.fasta")
for file in "${list[@]}"; do
    base="${file%_genes.fasta}"  
    kraken2 --db /mnt/d/DataBase/ --threads 40 --classified-out $outdir_nitrate/"$base"_classified.tsv \
    --output $outdir_nitrate/"$base"_output.txt --report $base/narG_report.txt \
    --unclassified-out $outdir_nitrate/"$base"_unclassifed.txt --use-names $file
done 

echo "Sulfate genes taxonomical assignment"
outdir_sulfate="Sulfate_taxonomical_annotations"
list2=("sat_genes.fasta" "aprA_genes.fasta" "aprB_genes.fasta" "dsrA_genes.fasta" "dsrB_genes.fasta")
for file in "${list2[@]}"; do
    base="${file%_genes.fasta}"  
    kraken2 --db /mnt/d/DataBase/ --threads 40 --classified-out $outdir_sulfate/"$base"_classified.tsv \
    --output $outdir_sulfate/"$base"_output.txt --report "$base"_report.txt \
    --unclassified-out $outdir_sulfate/"$base"_unclassifed.txt --use-names $file
```
# Seventh: data tables generation and quality check
Then, we're going to format the output results into a tabular format. 

```
echo "Nitrate table preparation"
cd "$outdir_nitrate"/ 
list=("narL_output.txt" "narX_output.txt" "narH_output.txt" "narI_output.txt" "narJ_output.txt")
for file in "${list[@]}"; do
    basename="${file%_output.txt}"
    awk '{print $3 }' "$file" | sort | uniq -c | sed 's/^[[:space:]]*//' | awk '{count=$1; $1=""; sub(/^ /, ""); print $0 "\t" count}' > "$basename"_per_specie.txt
    (echo -e "Taxa\t$basename" && cat "$basename"_per_specie.txt) > tmp ; mv tmp "$basename"_per_specie.txt
done

echo "Sulfate table preparation"
cd "$outdir_sulfate"/ 
list3=("sat_output.txt" "aprA_output.txt" "aprB_output.txt" "dsrA_output.txt" "dsrB_output.txt")
for file in "${lis3[@]}"; do
    basename="${file%_output.txt}"
    awk '{print $3 }' "$file" | sort | uniq -c | sed 's/^[[:space:]]*//' | awk '{count=$1; $1=""; sub(/^ /, ""); print $0 "\t" count}' > "$basename"_per_specie.txt
    (echo -e "Taxa\t$basename" && cat "$basename"_per_specie.txt) > tmp ; mv tmp "$basename"_per_specie.txt
done
```

Lastly, we're going to merge all file into two singular dataframes, one for each metabolic pathway.
To do so, it is only requiered to have all "$basename"_per_specie.txt files on same director and run the python script provided
Additionally, an automatic filtering step is introduce here within the same python script.
The reason behiend it is that we need to make sure that all genes are coming from taxa present only on the producton water (PW)
That's why on step third step the PW_taxonomy.bracken was created.

```
python merge_filter.py
```
There are fout files generated as output from this script: 
"Sulfate_locus_per_taxa.txt"
"Nitrate_locus_per_taxa.txt"
"Nitrate_locus_filtered.txt" 
"Sulfate_locus_filtered.txt"

However, the "Nitrate_locus_filtered.txt" it is not correct to use it as the script only trys to compare gene annotations at 
the gene level and in the "Nitrate_locus_filtered.txt" file there are genes that couldn't be annotated at gene level.
So, instead one have to manually compare  the "Nitrate_locus_per_taxa.txt" file against the "PW_taxonomy.bracen" and remove
all the rows that are not presented on "PW_taxonomy.bracen" file. 
To that last filterd file we recommend call it "Nitrate_locus_per_taxa_filtered_manual.txt" so that plot generation script does not break.

# Eight: Gene presence plot 
"Nitrate_locus_per_taxa_filtered_manual.txt" and "Sulfate_locus_filtered.txt" have to be on same folder to plot the results
used as figure 4 on the article. 
```
python gene_presence_plot.py
```
# Tenth: Swarm plot
In order to recreate the Figure 1. presented on the manuscript that displays the results from the Posgate C medium
we have to execute the following:

```
python swarmplot.py
```

The previous script already have the results obtained and the treatments assgined.
# Eleventh: Microbial Network construction
We're going to use the output generated from bracken at third step (PW_taxonomy.bracken and IW_taxonomy.bracken)
They must to be place on the same directory
## mOTU table generation
The output from this script is a tab separeted table.
```
python bracken_to_mOTU.py -i PW_taxonomy.bracken and IW_taxonomy.bracken -o PWIW_G_mOTU.tsv
```
## Samples simulation and distribution selection
With the following script we have the ability to simulate technical replicates for microbial samples, filter low-abundance OTUs, 
and perform PCA, PCoA (Bray-Curtis), and UMAP analysis.
The option that were used on this work are:
```
python replicates_simulation.py -h 
  --input INPUT         Path to the input TSV mOTU table file (e.g., PWIW_motu.tsv).
  --output_simulated_data OUTPUT_SIMULATED_DATA
                        Path to save the new mOTU table with simulated replicates.
  --num_replicates NUM_REPLICATES
                        Number of technical replicates to generate per sample (default: 3).
  --distribution {normal,poisson,negative_binomial}
                        Statistical distribution for generating variation ('normal', 'poisson', 'negative_binomial'). Default: 'normal'.
  --variation_factor VARIATION_FACTOR
                        Degree of variation parameter (float): - 'normal': Standard deviation (e.g., 0.1 for ~10% variation). - 'poisson': Not
                        applicable (variance tied to mean). - 'negative_binomial': Dispersion parameter (alpha), higher value means more
                        dispersion. Default: 0.1.
  --random_seed RANDOM_SEED
                        Random seed for reproducibility (default: 42).
```
The script was ran once for each distribution ('normal', 'poisson', 'negative_binomial') with --num_replicates 200, --variation_factor 0.05 and --random_seed 57.

For example: 

```
python replicates_simulation.py --input PWIW_G_mOTU.tsv \
 --output PWIW_simulated_neg_binom_05_rsed57_n200.tsv \
 --distribution negative_binomial --num_replicates 200 \
 --variation_factor 0.05 --random_seed 57
```

The output file generated "PWIW_simulated_neg_binom_05_rsed57_n200.tsv" is going to be use for the network generation

## Network construction
THe output from previous simulation will be used for the network generation.
The following R packages are requiered to be installed:
tidyverse (Wickham et al. (2019)),  NetCoMi (Peschel et al. (2020)), igraph (Antonov M et al. (2023)) and phyloseq (McMurdie PJ & Holmes S. (2013))

Make sure that change the path accordingly on the script:
```
file_paths <- list.files(path = "Simulations_genera_level/", #It also can be full path
                         pattern = "PWIW_simulated_neg_binom_05_rsed57_n200.tsv", # name of the
                         full.names = TRUE)
```
Execution:
```
Rscript microbial_network.R
```
Althought a graphical interface instead of terminal is recommended

# References
Wickham et al., (2019). Welcome to the Tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
Peschel, S., Müller, C. L., Von Mutius, E., Boulesteix, A., & Depner, M. (2020). NetCoMi: network construction and comparison for microbiome data in R. Briefings in Bioinformatics, 22(4). https://doi.org/10.1093/bib/bbaa290
Antonov M, Csárdi G, Horvát S, Müller K, Nepusz T, Noom D, Salmon M, Traag V, Welles BF, Zanini F (2023). “igraph enables fast and robust network analysis across programming languages.” arXiv preprint arXiv:2311.10260. doi:10.48550/arXiv.2311.10260.
McMurdie PJ, Holmes S (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE 8(4): e61217. https://doi.org/10.1371/journal.pone.0061217

