# First, reads download

''
# Data ID
IW="SRR14424469"
PW="SRR14424468"

woriking_dir="home"
#1 SRA-Toolkit for sequences retrieval

prefetch -O $woriking_dir/IW -t $IW
fasterq-dump $IW.sra

prefetch -O $woriking_dir/IW -t $IW
fasterq-dump $PW.sra
''
