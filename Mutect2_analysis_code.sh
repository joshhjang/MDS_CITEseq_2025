cd ./MUTECT2

module load bbc2/gatk/gatk-4.6.2.0
module load bbc2/samtools/samtools-1.21
module load bbc2/bcftools/bcftools-1.23 

# Reference + resources 
REF="./genomes/hg38/hg38.fa"     # Must match alignment genome build
GERMLINE_RESOURCE="./genomes/hg38/af-only-gnomad.hg38.vcf.gz"   # tumor-only support
COMMON_SITES_VCF="$GERMLINE_RESOURCE"    
BAM_DIR="./bamfiles"
THREADS=16

mkdir ./tmp
mkdir ./preprocess
mkdir ./vcf

SAMPLES=(
118877P2
118877SC
119325P2
119325SC
119837P2
119837P4
119837SC
)
	

for SAMPLE in "${SAMPLES[@]}"; do
	bam="${BAM_DIR}/${SAMPLE}.bam"
	echo $SAMPLE
	sortedbam="$SAMPLE.sorted.bam"
	samtools sort -@ 16 -o "$sortedbam" "$bam"
	samtools index -@ 16 "$sortedbam"

	#need to add read group
	BAM_RG="${SAMPLE}.rg.bam"
	gatk AddOrReplaceReadGroups \
	  -I "$sortedbam" \
	  -O "$BAM_RG" \
	  -RGID "$SAMPLE" \
	  -RGLB "lib1" \
	  -RGPL "ILLUMINA" \
	  -RGPU "unit1" \
	  -RGSM "$SAMPLE" \
	  --CREATE_INDEX true

	# -------------------------
	# 1) (Optional) MarkDuplicates
	# -------------------------
	BAM_MD="${SAMPLE}.md.bam"
	gatk --java-options "-Djava.io.tmpdir=./tmp" MarkDuplicates \
	    -I "$BAM_RG" \
	    -O "$BAM_MD" \
	    -M "./preprocess/${SAMPLE}.markdups.metrics.txt" \
	    --CREATE_INDEX true \
		2>&1 | tee "./preprocess/01_markdups.log"

	# -------------------------
	# 2) SplitNCigarReads (REQUIRED for RNA-seq)
	# -------------------------

	BAM_SPLIT="./preprocess/${SAMPLE}.split.bam"
	if [[ ! -f "$BAM_SPLIT" ]]; then
	  gatk --java-options "-Djava.io.tmpdir=./tmp" SplitNCigarReads \
	    -R "$REF" \
	    -I "$BAM_MD" \
	    -O "$BAM_SPLIT" \
	    --create-output-bam-index true \
	    2>&1 | tee "./preprocess/02_splitncigarreads.log"
	fi
	done
	
	  # -------------------------
	  # 3) Mutect2 (RNA mode) forcing known alleles
	  # -------------------------
	  UNFILTERED_VCF="./vcf/${SAMPLE}.unfiltered.vcf.gz"
	  gatk --java-options "-Xmx16g -Djava.io.tmpdir=./tmp" Mutect2 \
	    -R "$REF" \
	    -I "$BAM_SPLIT" \
	    -tumor "$SAMPLE" \
	    --native-pair-hmm-threads "$THREADS" \
	    --germline-resource "$GERMLINE_RESOURCE" \
	    --af-of-alleles-not-in-resource 0.0000025 \
	    -O "$UNFILTERED_VCF" \
	    2>&1 | tee "./logs/${SAMPLE}.03_Mutect2.log"

	  # -------------------------
	  # 4) FilterMutectCalls (populate FILTER labels)
	  # -------------------------
	  FILTERED_VCF="./vcf/${SAMPLE}.filtered.vcf.gz"
	  gatk --java-options "-Xmx8g -Djava.io.tmpdir=./tmp" FilterMutectCalls \
	    -R "$REF" \
	    -V "$UNFILTERED_VCF" \
	    -O "$FILTERED_VCF" \
	    2>&1 | tee "./logs/${SAMPLE}.04_FilterMutectCalls.log"

	  # -------------------------
	  # Extract VAF table from FILTERED VCF (NO exclusion)
	  #   - FILTER labels populated
	  #   - keep PASS and non-PASS
	  # -------------------------
	  VAF_TSV="./vcf/${SAMPLE}.known_mutations.vaf.tsv"

	  bcftools query \
	    -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%DP\t%AD\t%AF]\n' \
	    "$FILTERED_VCF" \
	    | awk 'BEGIN{OFS="\t"}
	           NR==1{print "CHROM","POS","REF","ALT","FILTER","DP","AD_ref","AD_alt","AF","VAF_from_AD"; next}
	           {
	             dp=$6;
	             split($7,ad,","); ad_ref=ad[1]; ad_alt=ad[2];
	             af=$8;
	             vaf="NA";
	             if (dp!="" && dp!="." && ad_alt!="" && ad_alt!="." && dp>0) vaf=ad_alt/dp;
	             print $1,$2,$3,$4,$5,dp,ad_ref,ad_alt,af,vaf
	           }' \
	    > "$VAF_TSV"

	  echo "Wrote: $VAF_TSV"
done
	
	
	
	
	
	
	