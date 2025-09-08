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
```
