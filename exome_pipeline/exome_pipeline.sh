## FASTQC

cd /mnt/bdata/shared/SF10711/exome/raw_data
for sample in $(ls | grep Sample.*); do
  for ea in $(ls "$sample" | grep .*fastq.gz); do 
    fastqc "$sample/$ea" --outdir fastqc
  done;
done;

## TRIM

cd /mnt/bdata/shared/SF10711/exome/raw_data
bbduk="/home/shared/programs/bbmap/bbduk.sh"
ref="/home/shared/programs/bbmap/resources/adapters.fa"
for sample in $(ls | grep Sample.*); do
  for ea in $(ls "$sample" *R1* | grep .*fastq); do
    R1="$sample/$ea"; R2=$(echo $R1 | sed "s/R1/R2/")
    $bbduk -Xmx1g in1=$R1 in2=$R2 \
    out1="$sample/$(echo $ea | sed s/.fastq/_trimmed.fastq/)" \
    out2="$sample/$(echo $(echo $ea | sed s/R1/R2/) | sed s/.fastq/_trimmed.fastq/)" \
    ref=$ref t=5 ktrim=r k=23 kmin=11 hdist=1 minlen=60 tpe tbo
  done;
done;

cd /mnt/bdata/shared/SF10711/exome/raw_data
for sample in $(ls | grep Sample.*); do
  cd "/mnt/bdata/shared/SF10711/exome/raw_data/$sample"
  pigz -p 5 $(ls *trimmed*)
done;
  
## FASTQC

cd /mnt/bdata/shared/SF10711/exome/raw_data
for sample in $(ls | grep Sample.*); do
  cd "/mnt/bdata/shared/SF10711/exome/raw_data/$sample"
  for ea in $(ls *trimmed*); do 
    fastqc $ea --outdir ../fastqc_trimmed
  done;
done;

## INDEX

ref="/home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa" 
/home/shared/programs/bwa-mem2/./bwa-mem2 index $ref

## ALIGN

ref="/home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa" 
bwa="/home/shared/programs/bwa-mem2/./bwa-mem2"
cd /mnt/bdata/shared/SF10711/exome/raw_data
for sample in $(ls | grep Sample.*) 
  do cd /mnt/bdata/shared/SF10711/exome/raw_data/$sample
    for ea in $(ls *R1* | grep .*trimmed.*); do 
      R1=$ea; R2=$(echo $R1 | sed "s/R1/R2/")
      $bwa mem -t 7 -T 0 $ref $R1 $R2 |
      samtools view --threads 7 -hb \
      -o ../../bam/$(echo $ea | sed "s/fastq.gz/bam/")
  done;
done;

## Add RG info

cd /mnt/bdata/shared/SF10711/exome/bam
picard="/home/shared/programs/Picard/picard.jar"
for sample in $(ls ../raw_data/ | grep Sample.* | sed "s/Sample_//"); do
  for lane in "L001" "L002"; do 
    for ea in $(ls "$sample"* | grep ".*$lane.*"); do 
      sem -j 5 java -jar $picard AddOrReplaceReadGroups \
      I=$ea O=temp/"$ea" RGLB=$sample RGPL=illumina \
      RGPU=$lane RGSM=$sample
    done;
  done;
done;

## MERGE & SORT

cd /mnt/bdata/shared/SF10711/exome/bam
picard="/home/shared/programs/Picard/picard.jar"
for sample in $(ls ../raw_data | grep Sample.* | sed "s/Sample_//"); do
  bamlist=$(for f in $(ls "$sample"*); do echo -n "I=$f "; done)
  sem -j 5 java -jar $picard MergeSamFiles $bamlist O=merged_sorted/"$sample"_sorted.bam
done

## ALIGN METRICS

cd /mnt/bdata/shared/SF10711/exome/bam/merged_sorted
picard="/home/shared/programs/Picard/picard.jar"
ref="/home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa" 
for ea in *bam; do 
  java -jar $picard CollectAlignmentSummaryMetrics \
  -R $ref -I $ea -O qc/$(echo $ea | sed "s/.bam//")_metrics.txt
done;

## MARK DUPS

cd /mnt/bdata/shared/SF10711/exome/bam/merged_sorted
picard="/home/shared/programs/Picard/picard.jar"
for ea in *bam; do
  sem -j 2 java -jar $picard MarkDuplicates -I $ea \
  -O mark_dups/"$(echo $ea | sed s/.bam//)"_rmdups.bam \
  -M mark_dups/"$(echo $ea | sed s/.bam//)"_dup_metrics.txt \
  --REMOVE_DUPLICATES true
done;

## BASE RECALIBRATION

## Index fasta:

picard="/home/shared/programs/Picard/picard.jar"
ref="/home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa" 
samtools faidx $ref
java -jar $picard CreateSequenceDictionary -R $ref

## Recal:

cd /mnt/bdata/shared/SF10711/exome/bam/merged_sorted/mark_dups
ref="/home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa" 
dbsnp="/home/shared/hg_align_db/dbSNP/Homo_sapiens_assembly38.dbsnp138.vcf"
gatk="/home/shared/programs/gatk-4.2.2.0/./gatk"

for ea in *bam; do
  $gatk BaseRecalibrator -I $ea -R $ref \
  --known-sites $dbsnp -O recal/"$(echo $ea | sed s/.bam//)"_recal.table
done;

for ea in *bam; do
  $gatk ApplyBQSR -R $ref -I $ea \
  --bqsr-recal-file recal/"$(echo $ea | sed s/.bam//)"_recal.table \
  -O recal/"$(echo $ea | sed s/.bam//)"_recal.bam
done;

## QC
## https://pmbio.org/module-03-align/0003/05/01/PostAlignment_QC/

cd /mnt/bdata/shared/SF10711/exome/bam/merged_sorted/mark_dups/recal
picard="/home/shared/programs/Picard/picard.jar"
ref="/home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa" 
for ea in *bam; do
  java -jar $picard CollectAlignmentSummaryMetrics \
  -R $ref -I $ea -O qc/$(echo $ea | sed "s/.bam//")_picard_metrics.txt
done;

for ea in *bam; do 
  java -jar $picard CollectInsertSizeMetrics -I $ea \
  -O qc/$(echo $ea | sed "s/.bam//")_picard_size_metrics.txt \
  -H qc/$(echo $ea | sed "s/.bam//")_picard_size_hist.pdf
done;

cd /mnt/bdata/shared/SF10711/exome/bam/merged_sorted/mark_dups/recal/qc
multiqc .

## Make BED file of exonic regions only:

cd /home/shared/hg_align_db/GRCh38_gencode_primary
grep 'transcript_type "protein_coding"' gencode.v38.primary_assembly.annotation.gtf |
awk '($3=="exon") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' | sort -T . -k1,1 -k2,2n | bedtools merge |
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,0,"."}' > gencode.v38.primary_assembly.annotation_exome.bed

## Format BED file:

cd /mnt/bdata/shared/SF10711/exome/design_files/SeqCapEZ_Exome_v3.0_Design_Annotation_files/
/opt/liftOver SeqCap_EZ_Exome_v3_hg19_capture_targets.bed /home/shared/hg_align_db/liftover/hg19ToHg38.over.chain SeqCap_EZ_Exome_v3_hg38_capture_targets.bed unMapped
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,0,"."}' SeqCap_EZ_Exome_v3_hg38_capture_targets.bed > SeqCap_EZ_Exome_v3_hg38_capture_targets_VEP_format.bed

## Get alignment metrics:

cd /mnt/bdata/shared/SF10711/exome/bam/merged_sorted/mark_dups/recal
target="/mnt/bdata/shared/SF10711/exome/design_files/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg38_capture_targets_VEP_format.bed"
qualimap="/home/shared/programs/qualimap_v2.2.1/./qualimap"
for ea in *bam; do
  $qualimap bamqc \
  --java-mem-size=20G -bam $ea \
  -outdir qc/qualimap/SeqCap \
  -outformat pdf \
  -outfile "$(echo $ea | sed s/.bam//)" \
  -nt 13 -gff $target
done;

## VARSCAN
## https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#somatic-variant-calling-workflow

cd /mnt/bdata/shared/SF10711/exome/bam/merged_sorted/mark_dups/recal
ref="/home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa"
varscan2="/home/shared/programs/VarScan2/VarScan.v2.4.4.jar"
normal="p300SF10711N_sorted_rmdups_recal.bam"
for tumor in SF*bam; do
  samtools mpileup -B -q 1 -Q 20 -f $ref $normal $tumor |
  java -jar $varscan2 somatic \
  --output vcf/"$(echo $tumor | sed s/.bam//)"_somatic \
  --mpileup 1 \
  --min-freq-for-hom 0.95 \
  --min-var-freq 0.05 \
  --p-value 0.1 \
  --output-vcf 1 
done;

cd /mnt/bdata/shared/SF10711/exome/bam/merged_sorted/mark_dups/recal/vcf
varscan2="/home/shared/programs/VarScan2/VarScan.v2.4.4.jar"
for ea in *vcf; do java -jar $varscan2 processSomatic $ea; done;

## VARIANT ANNOTATION

cd /mnt/bdata/shared/SF10711/exome/bam/merged_sorted/mark_dups/recal/vcf
ref="/home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa"
vep="/opt/vep_ensembl/ensembl-vep/vep"
for ea in *Somatic.hc.vcf; do
  vep -i $ea -o vep/"$(echo $ea | sed s/.vcf//)"_VEP.txt \
  --cache --dir_cache "/home/shared/vep_cache/" --assembly GRCh38 \
  --refseq --fasta $ref --hgvs --hgvsg --symbol --pick --tab --fork 1 \
  --offline
done;

