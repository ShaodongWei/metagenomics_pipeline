#!/bin/bash

################################################################################################################
############################################## help page #######################################################
Help()
{
   # Display Help
   echo "Add your parameters here."
   echo "options:"
   echo "-i|--input   Path to the raw fastq files directory"
   echo "-o|--output    Path to the output directory"
   echo "-t|--threads   Number of threads to use"   
   echo "-a|--assembly    Use "metaspades" [default] or "megahit" tool for assembling. "megahit" require less computation time and memory, but with worse assembling quality."
   echo "-b|--binning   If use "-b" or "--binning", binning will be done."
   echo "-H|--humann    If use "-H" or "--humann", functional profiling with humann will be done."   
   echo "--gzip   If use "--gzip", clean_fastq files will be gzipped."
}

################################################################################################################
############################################## set parameters ##################################################
## these parameters are for pipeline developers 
db=/home/yuanzh/shaodong/database # databases for kraken2, metaphlan4, checkm, 
sfw=/home/yuanzh/shaodong/softwares # softwares includign salmon, bowtie2, diamond, etc. 
keepbothreads="true" # !!! we use true as the default, though trimmomatic use false as default. If forward_only_surviving_pct is high, means our shotgun insert size is smaller than 150bp, so we have many reads forward and reverse are 100% same, so by default (false) trimmomatic only keep the forward, check the keepbothreads
mem=50 #in GB
parameter_a="metaspades"
parameter_b="false"
parameter_H="false"
parameter_gz="false"

##
VALID_ARGS=$(getopt -o hi:o:t:a:b::H:: --long input:,output:,threads:,assembly:,binning::,humann::,gzip:: -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -h|--help)
     Help
     exit
     ;; 
    -i|--input)
    parameter_i=$2
    shift 2 
    ;;
    -o|--output)
    parameter_o=$2
    shift 2
    ;;
    -t|--threads)
    parameter_t=$2
    shift 2
    ;;
    -a|--assembly)
    parameter_a=$2
    shift 2
    ;;
    -b|--binning) 
    parameter_b="true"
    shift 2
    ;;
    -H|--humann)
    parameter_H="true"
    shift 2
    ;;
    --gzip)
    parameter_gz="true"
    shift 2
    ;;
    --) 
    shift;
    break;;
  esac
done

############################################## check parameters 
if [[ ! -z $parameter_i ]]; then echo "Your input fastq is in $parameter_i"; else echo -e "\nYour input fastq is missing !!!"; echo "use -h for help"; exit;fi
if [[ ! -d $parameter_i ]]; then echo -e "\n$parameter_i does not exist !!!"; exit;fi
if [ ! -z $parameter_o ]; then echo "Your out directory is $parameter_o"; else echo -e "\nYour output directory is missing !!!"; echo "use -h for help"; exit;fi
if [[ ! -d $parameter_o ]]; then echo -e "\n$parameter_o does not exist !!!"; exit;fi
if [[ ! -z $parameter_t ]]; then echo "You will use $parameter_t threads"; else echo -e "\nNumber of threads is missing !!!"; echo "use -h for help"; exit;fi
if [[ $parameter_b ]]; then : ; else echo "Binning will be run."; fi
if [[ $parameter_H ]]; then : ; else echo "Humann3 will be run to profile function and taxonomy."; fi
if [[ $parameter_a != "metaspades" || $parameter_a != "megahit"]]; then echo "Wrong typing. Choose "metaspades" or "megahit"."; exit; fi

################################################################################################################
############################################## load functions ##################################################
export PATH=$PATH:$sfw/minimap2-2.24_x64-linux # we need to map reads to assemblies, minimap2 is 3 times faster than bowtie2
export PATH=$PATH:$sfw/SPAdes-3.15.4-Linux/bin
export PATH=$PATH:$sfw/MEGAHIT-1.2.9-Linux-x86_64-static/bin/
export PATH=$PATH:$sfw/salmon-1.9.0_linux_x86_64/bin  # cannot use conda, since conda only installs very low version 0.14.2, donwload binary from the salmon webpage. change salmon to salmon_tmp, use salmon_tmp as the current version of salmon
export PATH=$PATH:$sfw/kraken-tools
export PATH=$PATH:$sfw/kraken2
export PATH=$PATH:$sfw/Bracken-2.7


################################################################################################################
############################################## step 1, checking files ##########################################
echo -e "step 1, checking files begins .......... \n"

############################################## check file numbers
file_num=$(ls $parameter_i/*1.fastq|sort|wc -l)

####################################### divide threads based on files numbers
# explanation: when file_number >= parameter_t, we let parallel $parameter_t jobs using 1 thread in each parallel; when file_number < parameter_t, we parallel all files with using ceiling(parameter_t/file_number) threads in each parallel.  
# this strategy aims to save runnting time, and is based on the observation that some tools although are given multiple threads, most of the time still running with 1 thread. 
if [[ $file_num -ge $parameter_t ]];then $job=$parameter_t; $job_threads=1; else $job=$file_num; $job_threads=$(( ($parameter_t+$file_num-1)/$file_num ));fi
# so we use $job and $job_threads in the parallel 

############################################## check file naming, reads naming 
## check if the unzipping & renmaing file & renaming reads has been done, if so, skip it 
check_gzip=$(ls $parameter_i|grep 'fastq$') # used for 
check_file=$(ls $parameter_i|grep '_1.fastq$')
check_reads=$(ls $parameter_i|head -n1|while read -r f1;do bioawk -cfastx '{print $name}' $parameter_i/$f1|head -n1|grep '/1$';done)

if [[ -z $check_gzip || -z $check_file || -z $check_reads ]]; then 

  ## check if it is gzipped 
  if [[ -z $check_gzip ]]; then
  parallel -j $job 'var={1}; pigz -p {2} -dc {1} > ${var%.gz*}' :::: <(ls $parameter_i/*) :::: <(echo $job_threads)
  fi

  ## change file names to something_1.fastq or something_1.fastq.gz.
  if [[ -z $check_file ]]; then 
  ls $parameter_i|sort|while read -r line;do
  if (echo $line|grep -q -E ".*_R1_.*fastq|.*_R1_.*fq|.*1.fq"); then
  new_name=$(echo $line|sed -e 's/_R1_.*fastq/_1.fastq/' -e 's/_R1_.*fq/_1.fastq/' -e 's/1.fq/1.fastq/' -e 's/R1.fq/1.fastq/')
  mv $parameter_i/$line $parameter_i/$new_name
  elif (echo $line|grep -q -E ".*_R2_.*fastq|.*_R2_.*fq|.*2.fq"); then
  new_name=$(echo $line|sed -e 's/_R2_.*fastq/_2.fastq/' -e 's/_R2_.*fq/_2.fastq/' -e 's/2.fq/2.fastq/' -e 's/R2.fq/2.fastq/')
  mv $parameter_i/$line $parameter_i/$new_name
  fi
  done

  ## change read names to @something/1 @something/2
  if [[ -z $check_reads ]]; then
  for i in $parameter_i/*1.fastq;do var=$(head -n1 $i); if [[ $var != *"/1" ]]; then sed '/^@/s/ /_/g' $i|bioawk -cfastx '{print "@"$name"/1""\n"$seq"\n""+""\n"$qual}' > $parameter_i/tmp; mv -f $parameter_i/tmp $i; fi; done # the name has to be something/1, something/2, for paired-end reads something has to be same. 
  for i in $parameter_i/*2.fastq;do var=$(head -n1 $i); if [[ $var != *"/1" ]]; then sed '/^@/s/ /_/g' $i|bioawk -cfastx '{print "@"$name"/2""\n"$seq"\n""+""\n"$qual}' > $parameter_i/tmp; mv -f $parameter_i/tmp $i; fi; done # the name has to be something/1, something/2, for paired-end reads something has to be same. 
  fi
else
  :
fi 

### check if forward and reverse files are matched
ls $parameter_i/*1.fastq|sort > $parameter_o/fastq.1.path
ls $parameter_i/*2.fastq|sort > $parameter_o/fastq.2.path
cat $parameter_o/fastq.1.path |rev|cut -d/ -f1|rev|sed 's/_1.fq//' > $parameter_o/sample.name
if [ $(wc -l $parameter_o/fastq.1.path|awk -F" "  '{print $1}') != $(wc -l $parameter_o/fastq.2.path|awk -F" "  '{print $1}') ];then echo "You have different number of forward and reverse files"; exit; fi # test if numbers are same 
awk -F"1.fastq" '{print $1}' $parameter_o/fastq.1.path > $parameter_o/f1
awk -F"2.fastq" '{print $1}' $parameter_o/fastq.2.path > $parameter_o/f2
if ! cmp -s $parameter_o/f1 $parameter_o/f2; then echo “warning !!! you have different forward and reverse names”; exit;fi 
rm -f $parameter_o/f1 $parameter_o/f2 $parameter_o/fastq.*.path $parameter_o/sample.name

echo -e "step 1, checking files finished .......... \n"

################################################################################################################
############################################## step 2, quality control #########################################
echo -e "step 2, quality control begins .......... \n"

##################################### identify sequencing adapters
for i in $sfw/Trimmomatic-0.39/adapters/*;do file=$(ls $parameter_i/*|head -n1); val=$(sed '/^>/d' $i|grep -c -f - $file); echo -e $i"\t"$val >> $parameter_o/adapt.check.out;done
adapt=$(sort -n -k2 $parameter_o/adapt.check.out |tail -n1|cut -f1 2>&1)
echo -e "Sequencing primer used is $adapt\n"
cat $parameter_o/adapt.check.out 
rm -f $parameter_o/adapt.check.out 

##################################### clean reads with kneaddata
source /home/yuanzh/anaconda3/etc/profile.d/conda.sh
conda activate kneaddata

mkdir -p $parameter_o/clean_fastq 

kneaddata --version
time parallel -j $job --xapply 'kneaddata -i1 {1} -i2 {2} -db {5}/kneaddata -o {4}/clean_fastq --sequencer-source {8} --bypass-trf\
--bowtie2 {10}/bowtie2-2.4.5-linux-x86_64 --bowtie2-options "--very-sensitive --dovetail --threads {6}" \
--trimmomatic {10}/Trimmomatic-0.39 \
--trimmomatic-options "ILLUMINACLIP:{11}:2:30:10:2:{9} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50" \
--remove-intermediate-output --reorder --threads {6} --max-memory {7} --output-prefix {3} \
--run-trim-repetitive --fastqc {10}/Fastqc' \
:::: <(ls $parameter_i/*1.fastq|sort) \
:::: <(ls $parameter_i/*2.fastq|sort) \
:::: <(ls $parameter_i/*1.fastq|sort|sed -e 's/\/.*\///' -e 's/_1.fastq//') \
:::: <(echo $parameter_o) \
:::: <(echo $db) \
:::: <(echo $job_threads) \
:::: <(echo $((1024*$mem))m) \
:::: <(echo ${adapt##*/}|cut -d- -f1) \
:::: <(echo $keepbothreads) \
:::: <(echo $sfw) \
:::: <(echo $adapt) #adapt has to be full path 

mkdir -p $parameter_o/multiqc_report
multiqc -d $parameter_o/clean_fastq/ -o $parameter_o/multiqc_report
perc=$(cat $parameter_o/multiqc_report/multiqc_data/multiqc_trimmomatic.txt |cut -f4|sed '1d'|awk '{sum=sum+$0}END{print sum/NR}')
echo $perc
perc=$(ceil $perc) # here I check if surviving reads are good enough > 80%
if (($perc<75)); then echo 'warning !!!, survived reads are too few, consider to use "keepbothreads" ';fi #if forward_only_surviving_pct is high, means our shotgun insert size is smaller than 150bp, so we have many reads forward and reverse are 100% same, so by default trimmomatic only keep the forward, check the keepbothreads

conda deactivate

echo -e "step 2, checking files finished .......... \n"

################################################################################################################
##################################### step 3, assembling #######################################################
echo -e "step 3, assembling begins .......... \n"

mkdir -p $parameter_o/assembling/

if [[ $parameter_a == "metaspades"]]; then 
  time parallel -j $job --xapply 'spades.py -o {4}/assembling/{3} --meta -1 {1} -2 {2} --disable-gzip-output -t {5} -m {6}' \
  :::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort) \
  :::: <(ls $parameter_o/clean_fastq/*paired_2.fastq|sort) \
  :::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort|sed 's/^\/.*clean_fastq\///'|sed 's/_paired_1.fastq//') \
  :::: <(echo $parameter_o) \
  :::: <(echo $job_threads) \
  :::: <(echo $mem)  # most of the time it is in multithreads running. Scaffold is recommended as the resulting sequences. 
  
  # add prefix to output files 
  ls $parameter_o/assembling/$assembling_mode|sort|while read -r line;do mv $parameter_o/assembling/$assembling_mode/${line}/contigs.fasta $parameter_o/assembling/$assembling_mode/${line}/${line}.contigs.fasta;done
  ls $parameter_o/assembling/$assembling_mode|sort|while read -r line;do mv $parameter_o/assembling/$assembling_mode/${line}/scaffolds.fasta $parameter_o/assembling/$assembling_mode/${line}/${line}.scaffolds.fasta;done
  ls $parameter_o/assembling/$assembling_mode|sort|while read -r line;do mv $parameter_o/assembling/$assembling_mode/${line}/spades.log $parameter_o/assembling/$assembling_mode/${line}/${line}.spades.log;done

else [[ $parameter_a == "megahit" ]]; then 
  time parallel -j $job --xapply 'megahit -1 {1} -2 {2} --out-prefix {3} -o {4}/assembling/{3} -t {5} --memory {6}' \
  :::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort) \
  :::: <(ls $parameter_o/clean_fastq/*paired_2.fastq|sort) \
  :::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort|sed 's/^\/.*clean_fastq\///'|sed 's/_paired_1.fastq//') \
  :::: <(echo $parameter_o) \
  :::: <(echo $job_threads) \
  :::: <(echo $mem)
fi

echo -e "step 3, assembling finished ...........\n"

################################################################################################################
###################################### step 4, binning with metadecoder (optional) #############################
if [[ $parameter_b ]]; then 

  echo -e "step 4, binning begins..........\n"

  ## use the pip in the ./local to install metadecoder, pip version 21.2.4 works, not the lastest 22.2.2
  mkdir -p $parameter_o/binning

  parallel -j $job --xapply 'minimap2 -ax sr -t {6}  {1} {2} {3} > {5}/binning/{4}.alignment.sam' \
  :::: <(find $parameter_o/assembling/$assembling_mode -name *.scaffolds.fasta|sort) \
  :::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort) \
  :::: <(ls $parameter_o/clean_fastq/*paired_2.fastq|sort) \
  :::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort|sed -e 's/\/.*\///' -e 's/_paired_1.fastq//') \
  :::: <(echo $parameter_o) \
  :::: <(echo $job_threads) # minimal length 2000 by default. 

  parallel -j $job --xapply 'metadecoder coverage -s {1}  -o {3}/binning/{2}.metadecoder.coverage --threads {4}' \
  :::: <(ls $parameter_o/binning/*|sort) \
  :::: <(ls $parameter_o/binning|sort|sed 's/.alignment.sam//') \
  :::: <(echo $parameter_o) \
  :::: <(echo $job_threads)

  parallel -j $job --xapply 'metadecoder seed --threads {4} -f {1} -o {3}/binning/{2}.metadecoder.seed ' \
  :::: <(find $parameter_o/assembling/$assembling_mode -name *.scaffolds.fasta|sort) \
  :::: <(find $parameter_o/assembling/$assembling_mode -name *.scaffolds.fasta|sort|sed -e 's/\/.*\///' -e 's/.scaffolds.fasta//') \
  :::: <(echo $parameter_o) \
  :::: <(echo $job_threads)

  mkdir -p $parameter_o/binning/bins

  parallel -j $job --xapply 'metadecoder cluster -f {1} -c {2} -s {3} -o {5}/binning/bins/{4}' \
  :::: <(find $parameter_o/assembling/$assembling_mode -name *.scaffolds.fasta|sort) \
  :::: <(ls $parameter_o/binning/*coverage|sort) \
  :::: <(ls $parameter_o/binning/*seed|sort) \
  :::: <( ls $parameter_o/binning/*coverage|sort|sed -e 's/\/.*\///' -e 's/.metadecoder.coverage//') \
  :::: <(echo $parameter_o) \
  :::: <(echo $job_threads) # -o is a file, not a directory, the default threads is ~25, I don’t know how to specify

  rm -f *metadecoder.kmers *metadecoder.dpgmm # it is located where the script is. 

  conda deactivate

  ###################################### binning quality check with checkm 
  conda activate checkm 
  mkdir -p $parameter_o/binning/checkm

  #checkm data setRoot $db/checkm # specify where the database is. 
  export CHECKM_DATA_PATH=$db/checkm # or use this 

  # or do checkm in only one directory,  it is same result as dividing bins 
  # input for checkm has to be a directory, not a file 
  time checkm lineage_wf $parameter_o/binning/bins $parameter_o/binning/checkm -x fasta -t $parameter_t --pplacer_threads $parameter_t -f $parameter_o/binning/checkm/checkm.out -q # pplacer is done bin by bin so very slow if not splitting bins into directories

  ###################################### bins taxonomy annotation
  # to be continued. !!!!!!!!!!


  conda deactivate

  echo -e "step 4, binning finished..........\n"
fi

################################################################################################################
##################################### step 5, taxonomic profiling ##############################################
echo -e "step 5, taxonomic profiling begins .......... \n"

mkdir -p $parameter_o/taxonomy/kraken2

kraken2 --version
time parallel -j $job --xapply 'kraken2 -db {5}/kraken2 --threads {6} --output {4}/taxonomy/kraken2/{3}.out \
--report {4}/taxonomy/kraken2/{3}.report --paired {1} {2}' \
:::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort) \
:::: <(ls $parameter_o/clean_fastq/*paired_2.fastq|sort) \
:::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort|sed 's/^\/.*clean_fastq\///'|sed 's/_paired_1.fastq//') \
:::: <(echo $parameter_o) \
:::: <(echo $db) \
:::: <(echo $job_threads) 

mkdir -p $parameter_o/taxonomy/bracken

time parallel -j $job --xapply 'bracken -d {5}/kraken2  -i {1} -o {4}/taxonomy/bracken/{2} -w {4}/taxonomy/bracken/{3} -r 150 -l S -t {6}' \
:::: <(ls $parameter_o/taxonomy/kraken2/*report|sort) \
:::: <(ls $parameter_o/taxonomy/kraken2/*out|sort|sed  -e 's/\/.*kraken\///' -e 's/\.out$//'|awk '{print $0".S.bracken"}') \
:::: <(ls $parameter_o/taxonomy/kraken2/*out|sort|sed  -e 's/\/.*kraken\///' -e 's/\.out$//'|awk '{print $0".S.bracken.report"}') \
:::: <(echo $parameter_o) \
:::: <(echo $db) \
:::: <(echo $job_threads) # Species level

time parallel -j $job --xapply 'bracken -d {5}/kraken2  -i {1} -o {4}/taxonomy/bracken/{2} -w {4}/taxonomy/bracken/{3} -r 150 -l G -t {6}' \
:::: <(ls $parameter_o/taxonomy/kraken2/*report|sort) \
:::: <(ls $parameter_o/taxonomy/kraken2/*out|sort|sed  -e 's/\/.*kraken\///' -e 's/\.out$//'|awk '{print $0".G.bracken"}') \
:::: <(ls $parameter_o/taxonomy/kraken2/*out|sort|sed  -e 's/\/.*kraken\///' -e 's/\.out$//'|awk '{print $0".G.bracken.report"}') \
:::: <(echo $parameter_o) \
:::: <(echo $db) \
:::: <(echo $job_threads) # Genus level

#use kraken tools to make metaphlan like structure abundance table mpa.txt, based on S.bracken
parallel -j30 --xapply 'kreport2mpa.py -r {1} --display-header --percentages  --no-intermediate-ranks -o {3}/taxonomy/bracken/{2}.mpa.txt' \
:::: <(ls $parameter_o/taxonomy/bracken/*S.bracken.report|sort) \
:::: <(ls $parameter_o/taxonomy/bracken/*S.bracken.report|sort|sed -e 's/\/.*\///' -e 's/.S.bracken.report//' ) \
:::: <(echo $parameter_o)

# merge mpa files
combine_mpa.py -i $parameter_o/taxonomy/bracken/*mpa.txt  -o $parameter_o/taxonomy/bracken/merged.bracken.abundance.txt

echo -e "taxonomic profiling finished…………………….\n"

################################################################################################################
##################################### step 6, functional profiling ##############################################
echo -e "step 6, functional profiling begins .......... \n"

mkdir -p $parameter_o/function/prokka

## prokka to predict and annotate genes 
time parallel -j $job --xapply 'prokka --outdir {2}/function/prokka/{4} --cpus {3} --prefix {4} --force  {1} > {2}/function/prokka/{4}/{4}.log 2>&1' \
:::: <(find $parameter_o/assembling -name *.scaffolds.fasta|sort) \
:::: <(echo $parameter_o) \
:::: <(echo $job_threads) \
:::: <(find $parameter_o/assembling -name *.scaffolds.fasta|sort|sed -e 's/^\/.*\///' -e 's/.scaffolds.fasta//') 

## salmon to quantify genes 
salmon --version

## salmon index
# concat genes predicted from prokka  
find $parameter_o -name *ffn|while read -r line;do var1=${line##*/};bioawk -v var2=$var1 -cfastx '{print ">"var2"__"$name"\n"$seq}' $line >> $parameter_o/function/prokka/concat/concat.prokka.ffn;done # make sure reads name are unique, use “__” double underscore to separate sample name and reads name 
# concat scafolds
find $parameter_o/assembling -name '*scaffolds.fasta'|sort|while read -r line;do var1=${line##*/}; var1=${var1%.scaffolds.fasta}; bioawk -v var2="$var1" -cfastx '{print ">"var2"__"$name"\n"$seq}' $line >> $out/assembling/concat/concat.scaffolds.fasta;done

mkdir -p $parameter_o/function/salmon
time salmon index -t $parameter_o/function/prokka/concat/concat.prokka.ffn -i $parameter_o/function/salmon/salmon_index -p $parameter_t # there are 3 ways to map, CDS only, partial decoy or full genome decoy. 

## salmon quant 
time parallel -j $job --xapply 'salmon quant -i {4}/function/salmon/salmon_index/ -l A -1 {1} -2 {2} -p {5} --meta -o {4}/function/salmon/salmon_quant/{3}' \
:::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort) \
:::: <(ls $parameter_o/clean_fastq/*paired_2.fastq|sort) \
:::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort|sed 's/^\/.*clean_fastq\///'|sed 's/_paired_1.fastq//') \
:::: <(echo $parameter_o) \
:::: <(echo $job_threads)

## merge samples
salmon quantmerge --quants $(find $parameter_o/function -name 'quant.sf'|sed 's/\/quant.sf//') --column TPM  -o $parameter_o/function/salmon/salmon_quant/merged.salmon.prokka.ffn.quant # merge genes from different samples 

echo -e "step 6, functional profiling finished .......... \n"

################################################################################################################
##################################### step 7, functional profiling with humann (optional) ######################
if [[ $parameter_H ]]; then 
  echo -e "step 7, functional profiling with humann3 begins .......... \n"
  
  conda activate humann3_mpa4
  
  mkdir -p $parameter_o/clean_fastq/concat # we should use clean reads

  ## concat 1.fastq 2.fastq for humann3
  parallel -j $cpu --xapply 'cat {1} {2} > {4}/clean_fastq/concat/{3}' \
  :::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort) \
  :::: <(ls $parameter_o/clean_fastq/*paired_2.fastq|sort) \
  :::: <(ls $parameter_o/clean_fastq/*paired_1.fastq|sort|sed -e 's/\/.*\///' -e 's/_paired_1.fastq/.fastq/') \
  :::: <(echo $parameter_o) # the best option for paired-end file is to concatenate forward and reverse. 

  #humann_databases --available # to check available databases. 
  #humann_databases --download uniref uniref90_diamond $db/humann3 #uniref90 is default for genes annotation
  #humann_databases --download chocophlan full $db/humann3 #database, build, directory. Chocophlan is the database which metaphlan uses. 
  #humann_config --update database_folders protein $db/humann3/uniref
  #humann_config --update database_folders nucleotide $db/humann3/chocophlan
  #metaphlan --install --bowtie2db $db/metaphlan3 --nproc $cpu # metaphlan4.beta is available, but database has to be compatible with command
  #export DEFAULT_DB_FOLDER=/mnt/d/shaodong/database/metaphlan4/ # this export seems not working, so still use --bowtie2db in the next step 
  #metaphlan file.fastq --input_type fastq -o profiled.metagenome.txt --bowtie2db /mnt/d/shaodong/database/metaphlan4 --nproc 10 


  mkdir -p $parameter_o/function/humann3
  mkdir -p $parameter_o/taxonomy/metaphlan4

  time parallel -j $job --xapply 'humann -i {1} -o {3}/function/humann3 --threads {5} --memory-use maximum --input-format fastq --metaphlan {4} \
   --metaphlan-options "-t rel_ab  --bowtie2db {6}/metaphlan4/ --nproc {5}" ' \
  :::: <(ls $parameter_o/clean_fastq/concat/*|sort) \
  :::: <(ls $parameter_o/clean_fastq/concat|sort|sed 's/.fastq//') \
  :::: <(echo $parameter_o) \
  :::: <(echo $sfw/metaphlan4) \
  :::: <(echo $job_threads)
  :::: <(echo $db)

  mkdir -p $parameter_o/function/humann3/merged/

  ## merge humann3 otuput
  humann_join_tables -i $parameter_o/function/humann3 -o $parameter_o/function/humann3/merged/merged_genefamilies.tsv --file_name genefamilies # merge humann tables
  humann_join_tables -i $parameter_o/function/humann3 -o $parameter_o/function/humann3/merged/merged_pathwayabundance.tsv --file_name pathabundance
  humann_join_tables -i $parameter_o/function/humann3 -o $parameter_o/function/humann3/merged/merged_pathwaycoverage.tsv --file_name pathcoverage

  ## copy and merge metaphlan4 output to $parameter_o/taxonomy/metaphlan4
  find $parameter_o/function/humann3/ -name *metaphlan_bugs_list.tsv|sort|while read -r line;do nam=$(echo $line|sed -e 's/\/.*\///' -e 's/_metaphlan_bugs_list.tsv//');\
  cp $line $parameter_o/taxonomy/metaphlan4/${nam}.metaphlan.profile.tsv;done # cp metaphlan.profile.tsv to taxonomy
  merge_metaphlan_tables.py $parameter_o/taxonomy/metaphlan4/*metaphlan.profile.tsv > $parameter_o/taxonomy/metaphlan4/merged.metaphlan.abundance.tsv

  conda deactivate

  echo -e "step 7, functional profiling with humann3 finished .......... \n"
fi

################################################################################################################
##################################### compress fastq files  ####################################################

## compress raw fastq input files 
pigz -p $parameter_t  $parameter_i/* 

## if use $parameter_gz, then compress clean_fastq files 
if $parameter_gz; then find $parameter_o/clean_fastq -type f -name '*fastq'|pigz -p $parameter_t 

echo -e "The pipeline is finished .......... \n"

