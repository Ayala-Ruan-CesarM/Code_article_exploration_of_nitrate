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
