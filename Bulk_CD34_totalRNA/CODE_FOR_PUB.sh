############ BULK CD34 total RNA ############
#trim reads for adapters
trim_galore --fastqc -o ../trimmedreads --paired "${line}_L000_R1_001.fastq.gz" "${line}_L000_R2_001.fastq.gz" > ../trimmedreads/"${line}_trimgalore_summary.txt"

#remove any rRNA reads
cd ../trimmedreads
mkdir ./no_rRNA/"${line}_sortmerna"
sortmerna -ref ../hg38/rRNA_SILVA/smr_v4.3_default_db.fasta -reads "${line}_L000_R1_001_val_1.fq.gz" -reads "${line}_L000_R2_001_val_2.fq.gz" -v -fastx -other ./no_rRNA/"${line}_no_rRNA" -out2 -paired_out -workdir ./no_rRNA/"${line}_sortmerna"
cd ./no_rRNA/
gzip "${line}_no_rRNA_fwd.fq"
gzip "${line}_no_rRNA_rev.fq"

#align reads to hg38 genome
hisat2 -x ./hg38/HISAT2/HISAT2 --rna-strandness RF --dta -1 "${line}_no_rRNA_fwd.fq.gz" -2 "${line}_no_rRNA_rev.fq.gz" -p 24 --summary-file ../../aligned_no_rRNA/HISAT2_stranded/summary/"${line}_HISAT_summary.txt" | samtools view -bS - | samtools sort -T ../../aligned_no_rRNA/HISAT2_stranded/tmp/"${line}.tmp" -O bam -o ../../aligned_no_rRNA/HISAT2_stranded/"${line}.bam"
cd ../../aligned_no_rRNA/HISAT2_stranded/
samtools index "${line}.bam"

#parse HISAT2 alignment QC (needed for TE normalization)
python3 HISAT2_parse_summary.py

##quantify gene expression/counts
mkdir stringtie
stringtie "${line}.bam" -p 24 -l ${line} -o ./stringtie/GTF/"${line}.gtf" --rf -G ../../hg38/gencode.v37.annotation.gtf -A ./stringtie/gene_abundance/"${line}.gene_abundance.txt" -e -b ./stringtie/ballgown/"${line}"
#generate gene count using "prepDE.py3"

#quantify TE expression/counts
featureCounts -s 2 -T 24 -t exon -g transcript_id -p -f -a ../../hg38/GRCh38_GENCODE_rmsk_TE.gtf -o ./featureCount/"${line}_rpms.tsv" "${line}.bam" --primary -Q 20 --ignoreDup 
#normalize TE counts using library size from HISAT2 output
python3 FeatureCount_TE_CPM_calculation.py




############ CITE-seq ############

#Run cellranger to get count matrix (cellranger-6.1.2)

cellranger count --id=MDS_Exp_Brep \
                   --libraries=library_GEX_ADT.csv  \
				   --transcriptome=refdata-gex-GRCh38-2020-A\
                   --feature-ref=CITEseq_220antibody_HTO_FeatureBarcode.csv\ #use CITEseq_220antibody_HTO_FeatureBarcode_NewCocktail.csv for EXP 6,7,8
                   --expect-cells=25000 \
					--chemistry=SC5P-PE \
                   --localcores=24 \
				   --localmem=128 

#Run CellBender on CellRanger output raw_feature_bc_matrix.h5
cellbender remove-background --cuda --input raw_feature_bc_matrix.h5 --output MDS_Exp_CellBender

#Run TRUST4 for TCR/BCR repretoire prediction
run-trust4 --abnormalUnmapFlag -b possorted_genome_bam.bam -o MDS_Exp -t 16 -f hg38_bcrtcr.fa --ref human_IMGT+C.fa --barcode CB --UMI UB
