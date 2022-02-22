#user:whhe  
#data:2021-07-20
#=========================input info==========================#
#vcf:input VCF file with phased GT, LRR, and BAF
#pfx: output prefix
#thr: number of threads to use
#crt: file with call rate information (first column sample ID, second column call rate)
#sex: file with computed gender information (first column sample ID, second column gender: 1=male; 2=female)
#xcl:VCF file with additional list of variants to exclude (optional)
#ped:pedigree file to use if parent child duos are present

#==========================run-script==========================#
#$1: input_dir 
#$2: output_dir
#=========================version==============================#
# version_1.0.1 :    advance in cyto site of cnv 
inputdir=$1
outdir=$2
call_rate=$3

bcftools +gtc2vcf -i -g $inputdir/data > $outdir/sample.txt
#bcftools +affy2vcf --cel --chps ./data/

#对idat 进行Gencall算法处理
LANG="en_US.UTF-8" $HOME/bin/iaap-cli/iaap-cli gencall  /dev/work/project/snp_850k/db/CytoSNP-850Kv1-2_NS550_D1.bpm /dev/work/project/snp_850k/db/CytoSNP-850Kv1-2_NS550_D1_ClusterFile.egt $outdir/01.vcf/ --idat-folder $inputdir/data/ --output-gtc --gender-estimate-call-rate-threshold -0.1

#转格式 
bcftools +gtc2vcf --no-version -Ou --bpm /dev/work/project/snp_850k/db/CytoSNP-850Kv1-2_NS550_D1.bpm --csv /dev/work/project/snp_850k/db/CytoSNP-850Kv1-2_NS550_D1_CSV.csv --egt /dev/work/project/snp_850k/db/CytoSNP-850Kv1-2_NS550_D1_ClusterFile.egt --gtcs $outdir/01.vcf/ --fasta-ref /dev/work/db/hg19/hg19.fa --extra $outdir/01.vcf/sample.tsv|bcftools sort -Ou -T ./bcftools-sort.XXXXXX |bcftools norm --no-version -Ob -o $outdir/01.vcf/sample.bcf -c x -f /dev/work/db/hg19/hg19.fa
bcftools index -f $outdir/01.vcf/sample.bcf

#低call-rate sample 处理
mkdir -p $outdir/00.qc
awk -F '\t' 'NR>1{gsub(/.gtc/,"",$1);print $1"\t"$22"\t"$21}' $outdir/01.vcf/sample.tsv |sed '1isample_id\tcomputed_gender\tcall_rate' > $outdir/00.qc/all.sample.info.xls
awk -F "\t" -vs=$call_rate '$3<s{print $1}' $outdir/00.qc/all.sample.info.xls > $outdir/00.qc/to_recheck.sample.xls

echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">'|bcftools annotate --no-version -Ou -a /dev/work/project/snp_850k/db/segdups.bed.gz -c CHROM,FROM,TO,JK -h /dev/stdin $outdir/01.vcf/sample.bcf |bcftools view --no-version -Ou -S ^$outdir/00.qc/to_recheck.sample.xls|bcftools +fill-tags --no-version -Ou -t ^Y,MT,chrY,chrM -- -t ExcHet,F_MISSING |bcftools view --no-version -Ou -G |bcftools annotate --no-version -Ob -o $outdir/01.vcf/sample.filter.bcf -i 'FILTER!="." && FILTER!="PASS" || INFO/JK<.02 || INFO/ExcHet<1e-6 || INFO/F_MISSING>1-.97' -x ^INFO/JK,^INFO/ExcHet,^INFO/F_MISSING 
bcftools index -f $outdir/01.vcf/sample.filter.bcf

#call filter cnv
mkdir -p $outdir/02.cnv
bcftools +mocha -g GRCh37 --input-stats $outdir/00.qc/all.sample.info.xls --no-version --output-type b --output $outdir/02.cnv/filter.bdev.bcf --variants $outdir/01.vcf/sample.filter.bcf --calls $outdir/02.cnv/sample.filter.calls.tsv --stats $outdir/02.cnv/filter.stats.tsv --ucsc-bed $outdir/02.cnv/filter.ucsc.bed --cnp /dev/work/project/snp_850k/db/cnp.grch37.bed --mhc 6:27486711-33448264 --kir 19:54574747-55504099 $outdir/01.vcf/sample.bcf

#call all cnv 
bcftools +mocha -g GRCh37 --input-stats $outdir/00.qc/all.sample.info.xls --no-version --output-type b --output  $outdir/02.cnv/all.bdev.bcf --calls $outdir/02.cnv/sample.all.calls.tsv --stats $outdir/02.cnv/all.stats.tsv --ucsc-bed $outdir/02.cnv/all.ucsc.bed --cnp /dev/work/project/snp_850k/db/cnp.grch37.bed --mhc 6:27486711-33448264 --kir 19:54574747-55504099 $outdir/01.vcf/sample.bcf
bcftools index $outdir/02.cnv/all.bdev.bcf
awk -F '\t' 'NR>1{$6=sprintf("%.2f",$6/1000000);gsub(/Gain/,"dup",$21);gsub(/Loss/,"del",$21);print $3"\t"$4"\t"$5"\t"$6"\t"$1"\tmosaic "$21"\t"$14"\t"$12}' $outdir/02.cnv/sample.all.calls.tsv|bedtools intersect -a - -b /dev/work/project/snp_850k/db/hg19_cytoBand.txt -wa -wb|awk -F '\t' '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8]=a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8]","$12}END{for(i in a){print i"\t"a[i]}}'|sed 's/,//' > $outdir/02.cnv/sample.tmp.vcf
perl /dev/work/project/snp_850k/scr/get_cyto.pl $outdir/02.cnv/sample.tmp.vcf $outdir/02.cnv/sample.tmp1.vcf
#dup del LOH mosaic-rate calculate

cat $outdir/02.cnv/sample.tmp1.vcf |sort -V|bedtools intersect -a - -b /dev/work/project/snp_850k/db/hg19_cytoBand_pq.new.txt  -wao| awk -F '\t' '{if($(NF-2)==-1 || $(NF-1)~/p/ && $(NF-5)~/p/ || $(NF-1)~/q/ && $(NF-6)~/q/){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}else{if($(NF-1)~/p/){print $1"\t"$2"\t"$(NF-2)"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$(NF-1)}else{print $1"\t"$(NF-3)"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$(NF-1)"\t"$10}}}'|bedtools intersect -a - -b /dev/work/project/snp_850k/db/hg19_snp_chr.xls -wao|awk -F '\t' '{if($2<$(NF-2)){$2=$(NF-2)}if($3>$(NF-1)){$3=$(NF-1)};for(i=1;i<11;i++){printf $i"\t"}print""}'|sed 's/\t$//g'|bedtools intersect -a - -b /dev/work/project/snp_850k/db/hg19.gene.bed -wao|awk -F '\t' '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10]=a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10]","$14}END{for(i in a){print i"\t"a[i]}}'|sed 's/,//'|awk -F '\t' '{if($6~/dup/){rt=sprintf("%.2f",$8/(2/3-1/2)*100)}else{rt=sprintf("%.2f",2*$8/(0.5+$8)*100)};if(rt=="0.00"){rt="NA"}if($6~/LOH/){cns=2;print $5"\t"$1"\t"$9"-"$10"\t("$2"-"$3")x"cns"\t"$4"Mb\t"$6"\t"rt"\t"$11}else{cns=sprintf("%d",$7);cne=cns+1;print $5"\t"$1"\t"$9"-"$10"\t("$2"-"$3")x"cns"-"cne"\t"$4"Mb\t"$6"\t"rt"\t"$11}}'|sort -V |sed '1isample\tchr\tcytoBand_region\tcnv\tlength\ttype\tmosaic_rate(%)\tgene' > $outdir/02.cnv/all.sample.cnv.result.xls

#-----call-hmm--advance--#
#cd $outdir/01.vcf/ && sh /dev/work/project/snp_850k/scr/Allelic_shift_pipeline.sh && cd ../ 

#generate path
for sample in `awk 'NR>1{print $1}' $outdir/00.qc/all.sample.info.xls`
do
     mkdir -p $outdir/02.cnv/$sample/cnv_plot
     awk -F '\t' -vs=$sample 'NR==1||$1==s{print $0}' $outdir/02.cnv/all.sample.cnv.result.xls >  $outdir/02.cnv/$sample/$sample.cnv.xls	
     cat $outdir/02.cnv/$sample/$sample.cnv.xls|sed -e 's/(//g' -e 's/)//g'|awk -F '\t' 'NR>1{split($3,a,"-");split($4,b,"-|x");print $1"\t"$2"\t"b[1]"\t"b[2]"\t"b[3]"\t"a[1]"\t"$3"\t"$4"\t"$6"\t"$7"\t"$5"\t"$8}'|sed '1isample\tchr\tstart\tend\tcn\tcyto_start\tcyto_region\tcnv_info\ttype\tmocha_rate\tlength\tgene' >  $outdir/02.cnv/$sample/$sample.cnv.new.xls

#plot sample

     cat /dev/work/project/snp_850k/db/hg19.chr.site.txt | while read line
     do
     chr=`echo $line|awk -F ':' '{print $1}'`	     
     mocha_plot.R --mocha --stats $outdir/02.cnv/all.stats.tsv --vcf $outdir/02.cnv/all.bdev.bcf --png $outdir/02.cnv/$sample/cnv_plot/$sample.chr$chr.png --samples $sample --regions $line --cytoband /dev/work/project/snp_850k/db/cytoBand.txt.gz
     done
done
	
summary_plot.R --stats $outdir/02.cnv/all.stats.tsv --calls $outdir/02.cnv/sample.all.calls.tsv --pdf $outdir/02.cnv/sample.all.calls.pdf
pileup_plot.R --cytoband /dev/work/project/snp_850k/db/cytoBand.txt.gz --stats $outdir/02.cnv/all.stats.tsv --calls $outdir/02.cnv/sample.all.calls.tsv --pdf /$outdir/02.cnv/sample.chr.cn.pdf
sed -n '2,$p' $outdir/02.cnv/all.sample.cnv.result.xls | while read line
do
	sample=`echo $line|awk '{print $1}'`
	site=`echo $line|sed  -e 's/(//g' -e 's/)//g'|awk -F '\t' '{split($4,a,"x");print $2":"a[1]}'`
        png_name=`echo "$sample_$site"|sed 's/:/_/g'`
	mocha_plot.R --mocha --stats $outdir/02.cnv/all.stats.tsv --vcf $outdir/02.cnv/all.bdev.bcf --png $outdir/02.cnv/$sample/cnv_plot/$png_name.png --samples $sample --regions $site --cytoband /dev/work/project/snp_850k/db/cytoBand.txt.gz
done

#all sample filter
awk -F "\t" 'NR==FNR && FNR==1 {for (i=1; i<=NF; i++) f[$i] = i}NR==FNR && FNR>1 && ($(f["call_rate"])<.97 || $(f["baf_auto"])>.03) {xcl[$(f["sample_id"])]++}NR>FNR && FNR==1 {for (i=1; i<=NF; i++) g[$i] = i; print}NR>FNR && FNR>1 {len=$(g["length"]); bdev=$(g["bdev"]); rel_cov=$(g["rel_cov"])}NR>FNR && FNR>1 && !($(g["sample_id"]) in xcl) && $(g["type"])!~"^CNP" &&( $(g["chrom"])~"X" && $(g["computed_gender"])=="M" || bdev<0.1 || $(g["n50_hets"])<2e5 ) &&( $(g["bdev_se"])!="nan" || $(g["lod_baf_phase"])!="nan" && $(g["lod_baf_phase"]) > 10.0 ) &&( rel_cov<2.1 || bdev<0.05 || len>5e5 && bdev<0.1 && rel_cov<2.5 || len>5e6 && bdev<0.15 )' $outdir/02.cnv/all.stats.tsv > $outdir/02.cnv/all.filter.stats.tsv
awk 'NR==FNR {x[$1"_"$3"_"$4"_"$5]++} NR>FNR && ($0~"^track" || $4"_"$1"_"$2"_"$3 in x)' $outdir/02.cnv/all.filter.stats.tsv $outdir/02.cnv/all.ucsc.bed > $outdir/02.cnv/all.ucsc.filter.bed

# collect report result 
mkdir -p $outdir/report_result
cp -r $outdir/02.cnv/*/ $outdir/report_result
cp -r $outdir/02.cnv/sample.all.calls.tsv $outdir/report_result/sample.all.calls.xls
cp -r $outdir/00.qc/*xls $outdir/report_result
tar 'fzc' $outdir/report_result.tar.gz $outdir/report_result
