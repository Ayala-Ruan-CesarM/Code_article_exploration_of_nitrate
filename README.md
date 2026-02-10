# First,  download reads with SRA-toolkit
```
IW="SRR14424469"
PW="SRR14424468"

woriking_dir="home"
mkdir -p IW
mkdir -p PW
prefetch -O $woriking_dir/IW/ $IW
fasterq-dump -s -e 20 -O $woriking_dir/ -o IW $woriking_dir/$IW/*.sra

prefetch -O $woriking_dir/IW $IW
fasterq-dump -s -e 20 -O $woriking_dir/ -o PW $woriking_dir/$PW/*.sra

```
# Second, quality check with fastp
```
cd $woriking_dir/IW/
mkdir -p fastp
fastp -i IW_1.fastq -o fastp/IW_1_qc.fastq -I IW_2.fastq -O fastp/IW_2_qc.fastq -w 20 -h report.html
cd $woriking_dir/PW/
fastp -i PW_1.fastq -o fastp/PW_1_qc.fastq -I PW_2.fastq -O fastp/PW_2_qc.fastq -w 20 -h report.html
```
# Third, taxonomical assigment 
Here, we need the path to the .k2d index files needed for kraken.
```
DB="/path/to/db_index"
kraken2 --db "$DB" --threads --output PW_taxonomy.txt --report PW_taxonomy.report --paired fastp/PW_1_qc.fastq fastp/PW_2_qc.fastq
```
Now, read correction with bracken
```
bracken -d $DB" -i PW_taxonomy.report -o PW_taxonomy.brack -w PW_taxonomy.brack.report -r 150 -l G
```
# Fourth, De-novo assembly and gene predciton 
Reads merged 
```
mkdir -p Assembly
cat fastp/IW_1_qc.fastq fastp/PW_1_qc.fastq > Assembly/Merge_1_Reads.fastq ;
cat fastp/IW_2_qc.fastq fastp/PW_2_qc.fastq > Assembly/Merge_2_Reads.fastq  
```
Spades assembly 
```
spades --meta 
```

