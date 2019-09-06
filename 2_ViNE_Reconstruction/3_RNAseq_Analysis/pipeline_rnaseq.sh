# Before running this script, I downloaded the Medicago truncatula A17 genome files from https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/ and the Sinorhizobium meliloti Rm2011 genome files from https://iant.toulouse.inra.fr/bacteria/annotation/cgi/rhime2011/rhime2011.cgi

# Download the sequencing data
sh Scripts/downloadData.sh

# Combine technical replicates into a single file
mv rawReads/SRR949227_1.fastq.gz rawReads/zoneI_rep1_F.fastq.gz
mv rawReads/SRR949227_2.fastq.gz rawReads/zoneI_rep1_R.fastq.gz
mv rawReads/SRR949228_1.fastq.gz rawReads/zoneI_rep2_F.fastq.gz
mv rawReads/SRR949228_2.fastq.gz rawReads/zoneI_rep2_R.fastq.gz
mv rawReads/SRR949229_1.fastq.gz rawReads/zoneI_rep3_F.fastq.gz
mv rawReads/SRR949229_2.fastq.gz rawReads/zoneI_rep3_R.fastq.gz
cat rawReads/SRR949230_1.fastq.gz rawReads/SRR949231_1.fastq.gz > rawReads/zoneIId_rep1_F.fastq.gz
cat rawReads/SRR949230_2.fastq.gz rawReads/SRR949231_2.fastq.gz > rawReads/zoneIId_rep1_R.fastq.gz
mv rawReads/SRR949232_1.fastq.gz rawReads/zoneIId_rep2_F.fastq.gz
mv rawReads/SRR949232_2.fastq.gz rawReads/zoneIId_rep2_R.fastq.gz
mv rawReads/SRR949233_1.fastq.gz rawReads/zoneIId_rep3_F.fastq.gz
mv rawReads/SRR949233_2.fastq.gz rawReads/zoneIId_rep3_R.fastq.gz
cat rawReads/SRR949234_1.fastq.gz rawReads/SRR949235_1.fastq.gz rawReads/SRR949236_1.fastq.gz > rawReads/zoneIIp_rep1_F.fastq.gz
cat rawReads/SRR949234_2.fastq.gz rawReads/SRR949235_2.fastq.gz rawReads/SRR949236_2.fastq.gz > rawReads/zoneIIp_rep1_R.fastq.gz
cat rawReads/SRR949237_1.fastq.gz rawReads/SRR949238_1.fastq.gz rawReads/SRR949239_1.fastq.gz > rawReads/zoneIIp_rep2_F.fastq.gz
cat rawReads/SRR949237_2.fastq.gz rawReads/SRR949238_2.fastq.gz rawReads/SRR949239_2.fastq.gz > rawReads/zoneIIp_rep2_R.fastq.gz
cat rawReads/SRR949598_1.fastq.gz rawReads/SRR949599_1.fastq.gz rawReads/SRR949600_1.fastq.gz > rawReads/zoneIIp_rep3_F.fastq.gz
cat rawReads/SRR949598_2.fastq.gz rawReads/SRR949599_2.fastq.gz rawReads/SRR949600_2.fastq.gz > rawReads/zoneIIp_rep3_R.fastq.gz
cat rawReads/SRR949240_1.fastq.gz rawReads/SRR949241_1.fastq.gz rawReads/SRR949242_1.fastq.gz > rawReads/zoneIZ_rep1_F.fastq.gz
cat rawReads/SRR949240_2.fastq.gz rawReads/SRR949241_2.fastq.gz rawReads/SRR949242_2.fastq.gz > rawReads/zoneIZ_rep1_R.fastq.gz
cat rawReads/SRR949243_1.fastq.gz rawReads/SRR949244_1.fastq.gz > rawReads/zoneIZ_rep2_F.fastq.gz
cat rawReads/SRR949243_2.fastq.gz rawReads/SRR949244_2.fastq.gz > rawReads/zoneIZ_rep2_R.fastq.gz
cat rawReads/SRR949245_1.fastq.gz rawReads/SRR949246_1.fastq.gz rawReads/SRR949247_1.fastq.gz > rawReads/zoneIZ_rep3_F.fastq.gz
cat rawReads/SRR949245_2.fastq.gz rawReads/SRR949246_2.fastq.gz rawReads/SRR949247_2.fastq.gz > rawReads/zoneIZ_rep3_R.fastq.gz
cat rawReads/SRR949260_1.fastq.gz rawReads/SRR949261_1.fastq.gz rawReads/SRR949262_1.fastq.gz > rawReads/zoneIII_rep1_F.fastq.gz
cat rawReads/SRR949260_2.fastq.gz rawReads/SRR949261_2.fastq.gz rawReads/SRR949262_2.fastq.gz > rawReads/zoneIII_rep1_R.fastq.gz
cat rawReads/SRR949263_1.fastq.gz rawReads/SRR949264_1.fastq.gz rawReads/SRR949265_1.fastq.gz > rawReads/zoneIII_rep2_F.fastq.gz
cat rawReads/SRR949263_2.fastq.gz rawReads/SRR949264_2.fastq.gz rawReads/SRR949265_2.fastq.gz > rawReads/zoneIII_rep2_R.fastq.gz
cat rawReads/SRR949267_1.fastq.gz rawReads/SRR949268_1.fastq.gz rawReads/SRR949269_1.fastq.gz > rawReads/zoneIII_rep3_F.fastq.gz
cat rawReads/SRR949267_2.fastq.gz rawReads/SRR949268_2.fastq.gz rawReads/SRR949269_2.fastq.gz > rawReads/zoneIII_rep3_R.fastq.gz
rm -r rawReads/SRR*.fastq.gz

# Prepare a combined genome with fasta and GFF files for mapping
mkdir combinedGenome/
head -43 downloadedGenomes/Medicago_genome.gff > combinedGenome/combinedGenome.gff
grep '>' downloadedGenomes/Sinorhizobium_genome.fasta | cut -f1,1 -d' ' | sed 's/>lcl|/##sequence-region   /' > combinedGenome/temp.txt
grep '>' downloadedGenomes/Sinorhizobium_genome.fasta | cut -f3,3 -d' ' | sed 's/len=/1 /' > combinedGenome/temp2.txt
paste -d' ' combinedGenome/temp.txt combinedGenome/temp2.txt >> combinedGenome/combinedGenome.gff
rm combinedGenome/temp*
grep 'ID=mRNA' downloadedGenomes/Medicago_genome.gff | cut -f9,9 | cut -f1,1 -d';' | sed 's/ID=mRNA://' > downloadedGenomes/Medicago_locusTags.txt
grep 'ID=mRNA' downloadedGenomes/Sinorhizobium_genome.gff | cut -f9,9 | cut -f1,1 -d';' | sed 's/ID=mRNA://' > downloadedGenomes/Sinorhizobium_locusTags.txt
grep -f 'downloadedGenomes/Medicago_locusTags.txt' downloadedGenomes/Medicago_genome.gff >> combinedGenome/combinedGenome.gff
grep -f 'downloadedGenomes/Sinorhizobium_locusTags.txt' downloadedGenomes/Sinorhizobium_genome.gff >> combinedGenome/combinedGenome.gff
sed -i 's/ID=gene:/gene_id=/' combinedGenome/combinedGenome.gff
sed -i 's/Parent=mRNA:/gene_id=/' combinedGenome/combinedGenome.gff
sed -i 's/"//g' combinedGenome/combinedGenome.gff
cat downloadedGenomes/Medicago_genome.fasta downloadedGenomes/Sinorhizobium_genome.fasta | sed 's/lcl|//' > combinedGenome/combinedGenome.fasta

# Index genome for bowtie2 mapping
mkdir expressionAnalysis/
mkdir expressionAnalysis/bowtieGenomeIndex/
bowtie2-build combinedGenome/combinedGenome.fasta expressionAnalysis/bowtieGenomeIndex/combinedGenome

# Map all samples to the combined genome
mkdir expressionAnalysis/bowtieOutput/
mkdir expressionAnalysis/bowtieSummaries/
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneI_rep1_F.fastq.gz -2 rawReads/zoneI_rep1_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneI_rep1.sam 2> expressionAnalysis/bowtieSummaries/zoneI_rep1.txt
gzip expressionAnalysis/bowtieOutput/zoneI_rep1.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneI_rep2_F.fastq.gz -2 rawReads/zoneI_rep2_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneI_rep2.sam 2> expressionAnalysis/bowtieSummaries/zoneI_rep2.txt
gzip expressionAnalysis/bowtieOutput/zoneI_rep2.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneI_rep3_F.fastq.gz -2 rawReads/zoneI_rep3_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneI_rep3.sam 2> expressionAnalysis/bowtieSummaries/zoneI_rep3.txt
gzip expressionAnalysis/bowtieOutput/zoneI_rep3.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIId_rep1_F.fastq.gz -2 rawReads/zoneIId_rep1_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIId_rep1.sam 2> expressionAnalysis/bowtieSummaries/zoneIId_rep1.txt
gzip expressionAnalysis/bowtieOutput/zoneIId_rep1.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIId_rep2_F.fastq.gz -2 rawReads/zoneIId_rep2_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIId_rep2.sam 2> expressionAnalysis/bowtieSummaries/zoneIId_rep2.txt
gzip expressionAnalysis/bowtieOutput/zoneIId_rep2.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIId_rep3_F.fastq.gz -2 rawReads/zoneIId_rep3_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIId_rep3.sam 2> expressionAnalysis/bowtieSummaries/zoneIId_rep3.txt
gzip expressionAnalysis/bowtieOutput/zoneIId_rep3.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIIp_rep1_F.fastq.gz -2 rawReads/zoneIIp_rep1_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIIp_rep1.sam 2> expressionAnalysis/bowtieSummaries/zoneIIp_rep1.txt
gzip expressionAnalysis/bowtieOutput/zoneIIp_rep1.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIIp_rep2_F.fastq.gz -2 rawReads/zoneIIp_rep2_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIIp_rep2.sam 2> expressionAnalysis/bowtieSummaries/zoneIIp_rep2.txt
gzip expressionAnalysis/bowtieOutput/zoneIIp_rep2.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIIp_rep3_F.fastq.gz -2 rawReads/zoneIIp_rep3_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIIp_rep3.sam 2> expressionAnalysis/bowtieSummaries/zoneIIp_rep3.txt
gzip expressionAnalysis/bowtieOutput/zoneIIp_rep3.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIZ_rep1_F.fastq.gz -2 rawReads/zoneIZ_rep1_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIZ_rep1.sam 2> expressionAnalysis/bowtieSummaries/zoneIZ_rep1.txt
gzip expressionAnalysis/bowtieOutput/zoneIZ_rep1.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIZ_rep2_F.fastq.gz -2 rawReads/zoneIZ_rep2_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIZ_rep2.sam 2> expressionAnalysis/bowtieSummaries/zoneIZ_rep2.txt
gzip expressionAnalysis/bowtieOutput/zoneIZ_rep2.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIZ_rep3_F.fastq.gz -2 rawReads/zoneIZ_rep3_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIZ_rep3.sam 2> expressionAnalysis/bowtieSummaries/zoneIZ_rep3.txt
gzip expressionAnalysis/bowtieOutput/zoneIZ_rep3.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIII_rep1_F.fastq.gz -2 rawReads/zoneIII_rep1_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIII_rep1.sam 2> expressionAnalysis/bowtieSummaries/zoneIII_rep1.txt
gzip expressionAnalysis/bowtieOutput/zoneIII_rep1.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIII_rep2_F.fastq.gz -2 rawReads/zoneIII_rep2_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIII_rep2.sam 2> expressionAnalysis/bowtieSummaries/zoneIII_rep2.txt
gzip expressionAnalysis/bowtieOutput/zoneIII_rep2.sam &
bowtie2 -p 20 -x expressionAnalysis/bowtieGenomeIndex/combinedGenome -1 rawReads/zoneIII_rep3_F.fastq.gz -2 rawReads/zoneIII_rep3_R.fastq.gz -S expressionAnalysis/bowtieOutput/zoneIII_rep3.sam 2> expressionAnalysis/bowtieSummaries/zoneIII_rep3.txt
gzip expressionAnalysis/bowtieOutput/zoneIII_rep3.sam &
wait

# Sort the SAM files
mkdir expressionAnalysis/sortedAlignmentFiles/
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneI_rep1.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneI_rep1.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneI_rep1.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneI_rep2.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneI_rep2.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneI_rep2.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneI_rep3.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneI_rep3.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneI_rep3.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIId_rep1.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIId_rep1.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIId_rep1.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIId_rep2.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIId_rep2.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIId_rep2.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIId_rep3.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIId_rep3.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIId_rep3.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIIp_rep1.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIIp_rep1.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIIp_rep1.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIIp_rep2.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIIp_rep2.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIIp_rep2.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIIp_rep3.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIIp_rep3.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIIp_rep3.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIZ_rep1.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIZ_rep1.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIZ_rep1.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIZ_rep2.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIZ_rep2.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIZ_rep2.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIZ_rep3.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIZ_rep3.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIZ_rep3.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIII_rep1.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIII_rep1.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIII_rep1.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIII_rep2.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIII_rep2.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIII_rep2.sam &
samtools sort -n -O SAM -@ 20 expressionAnalysis/bowtieOutput/zoneIII_rep3.sam.gz > expressionAnalysis/sortedAlignmentFiles/zoneIII_rep3.sam
gzip expressionAnalysis/sortedAlignmentFiles/zoneIII_rep3.sam &
wait

# Count reads per gene
mkdir expressionAnalysis/htseqCount/
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneI_rep1.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneI_rep1.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneI_rep2.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneI_rep2.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneI_rep3.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneI_rep3.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIId_rep1.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIId_rep1.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIId_rep2.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIId_rep2.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIId_rep3.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIId_rep3.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIIp_rep1.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIIp_rep1.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIIp_rep2.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIIp_rep2.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIIp_rep3.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIIp_rep3.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIZ_rep1.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIZ_rep1.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIZ_rep2.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIZ_rep2.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIZ_rep3.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIZ_rep3.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIII_rep1.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIII_rep1.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIII_rep2.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIII_rep2.txt &
htseq-count expressionAnalysis/sortedAlignmentFiles/zoneIII_rep3.sam.gz combinedGenome/combinedGenome.gff --type=gene > expressionAnalysis/htseqCount/zoneIII_rep3.txt &
wait

# Split count tables into species specific tables
mkdir expressionAnalysis/htseqCountMedicago/
mkdir expressionAnalysis/htseqCountSinorhizobium/
grep 'Mtrun' expressionAnalysis/htseqCount/zoneI_rep1.txt > expressionAnalysis/htseqCountMedicago/zoneI_rep1.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneI_rep2.txt > expressionAnalysis/htseqCountMedicago/zoneI_rep2.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneI_rep3.txt > expressionAnalysis/htseqCountMedicago/zoneI_rep3.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIId_rep1.txt > expressionAnalysis/htseqCountMedicago/zoneIId_rep1.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIId_rep2.txt > expressionAnalysis/htseqCountMedicago/zoneIId_rep2.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIId_rep3.txt > expressionAnalysis/htseqCountMedicago/zoneIId_rep3.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIIp_rep1.txt > expressionAnalysis/htseqCountMedicago/zoneIIp_rep1.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIIp_rep2.txt > expressionAnalysis/htseqCountMedicago/zoneIIp_rep2.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIIp_rep3.txt > expressionAnalysis/htseqCountMedicago/zoneIIp_rep3.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIZ_rep1.txt > expressionAnalysis/htseqCountMedicago/zoneIZ_rep1.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIZ_rep2.txt > expressionAnalysis/htseqCountMedicago/zoneIZ_rep2.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIZ_rep3.txt > expressionAnalysis/htseqCountMedicago/zoneIZ_rep3.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIII_rep1.txt > expressionAnalysis/htseqCountMedicago/zoneIII_rep1.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIII_rep2.txt > expressionAnalysis/htseqCountMedicago/zoneIII_rep2.txt &
grep 'Mtrun' expressionAnalysis/htseqCount/zoneIII_rep3.txt > expressionAnalysis/htseqCountMedicago/zoneIII_rep3.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneI_rep1.txt > expressionAnalysis/htseqCountSinorhizobium/zoneI_rep1.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneI_rep2.txt > expressionAnalysis/htseqCountSinorhizobium/zoneI_rep2.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneI_rep3.txt > expressionAnalysis/htseqCountSinorhizobium/zoneI_rep3.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIId_rep1.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep1.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIId_rep2.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep2.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIId_rep3.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep3.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIIp_rep1.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep1.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIIp_rep2.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep2.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIIp_rep3.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep3.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIZ_rep1.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep1.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIZ_rep2.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep2.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIZ_rep3.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep3.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIII_rep1.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep1.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIII_rep2.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep2.txt &
grep 'SM' expressionAnalysis/htseqCount/zoneIII_rep3.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep3.txt &
wait

# Get gene lengths
mkdir expressionAnalysis/geneLength/
perl Scripts/getSizes.pl combinedGenome/combinedGenome.gff > expressionAnalysis/geneLength/allGeneLengths.txt
grep 'Mtrun' expressionAnalysis/geneLength/allGeneLengths.txt > expressionAnalysis/geneLength/MedicagoGeneLengths.txt &
grep 'SM' expressionAnalysis/geneLength/allGeneLengths.txt > expressionAnalysis/geneLength/SinorhizobiumGeneLengths.txt &
wait

# Sort HT-seq output and gene lengths output
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneI_rep1.txt > expressionAnalysis/htseqCountMedicago/zoneI_rep1_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneI_rep2.txt > expressionAnalysis/htseqCountMedicago/zoneI_rep2_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneI_rep3.txt > expressionAnalysis/htseqCountMedicago/zoneI_rep3_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIId_rep1.txt > expressionAnalysis/htseqCountMedicago/zoneIId_rep1_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIId_rep2.txt > expressionAnalysis/htseqCountMedicago/zoneIId_rep2_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIId_rep3.txt > expressionAnalysis/htseqCountMedicago/zoneIId_rep3_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIIp_rep1.txt > expressionAnalysis/htseqCountMedicago/zoneIIp_rep1_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIIp_rep2.txt > expressionAnalysis/htseqCountMedicago/zoneIIp_rep2_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIIp_rep3.txt > expressionAnalysis/htseqCountMedicago/zoneIIp_rep3_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIZ_rep1.txt > expressionAnalysis/htseqCountMedicago/zoneIZ_rep1_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIZ_rep2.txt > expressionAnalysis/htseqCountMedicago/zoneIZ_rep2_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIZ_rep3.txt > expressionAnalysis/htseqCountMedicago/zoneIZ_rep3_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIII_rep1.txt > expressionAnalysis/htseqCountMedicago/zoneIII_rep1_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIII_rep2.txt > expressionAnalysis/htseqCountMedicago/zoneIII_rep2_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountMedicago/zoneIII_rep3.txt > expressionAnalysis/htseqCountMedicago/zoneIII_rep3_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneI_rep1.txt > expressionAnalysis/htseqCountSinorhizobium/zoneI_rep1_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneI_rep2.txt > expressionAnalysis/htseqCountSinorhizobium/zoneI_rep2_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneI_rep3.txt > expressionAnalysis/htseqCountSinorhizobium/zoneI_rep3_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep1.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep1_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep2.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep2_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep3.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep3_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep1.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep1_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep2.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep2_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep3.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep3_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep1.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep1_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep2.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep2_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep3.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep3_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep1.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep1_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep2.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep2_sorted.txt &
sort -k1,1 expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep3.txt > expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep3_sorted.txt &
sort -k1,1 expressionAnalysis/geneLength/MedicagoGeneLengths.txt > expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt &
sort -k1,1 expressionAnalysis/geneLength/SinorhizobiumGeneLengths.txt > expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt &
wait

# Calculate TPM values
mkdir expressionAnalysis/TPMvaluesMedicago/
mkdir expressionAnalysis/TPMvaluesSinorhizobium/
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneI_rep1_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneI_rep1.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneI_rep2_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneI_rep2.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneI_rep3_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneI_rep3.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIId_rep1_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIId_rep1.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIId_rep2_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIId_rep2.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIId_rep3_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIId_rep3.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIIp_rep1_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIIp_rep1.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIIp_rep2_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIIp_rep2.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIIp_rep3_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIIp_rep3.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIZ_rep1_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIZ_rep1.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIZ_rep2_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIZ_rep2.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIZ_rep3_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIZ_rep3.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIII_rep1_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIII_rep1.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIII_rep2_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIII_rep2.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountMedicago/zoneIII_rep3_sorted.txt expressionAnalysis/geneLength/MedicagoGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesMedicago/zoneIII_rep3.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneI_rep1_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep1.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneI_rep2_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep2.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneI_rep3_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep3.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep1_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep1.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep2_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep2.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIId_rep3_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep3.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep1_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep1.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep2_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep2.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIIp_rep3_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep3.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep1_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep1.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep2_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep2.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIZ_rep3_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep3.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep1_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep1.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep2_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep2.txt
perl Scripts/calculateTPM.pl expressionAnalysis/htseqCountSinorhizobium/zoneIII_rep3_sorted.txt expressionAnalysis/geneLength/SinorhizobiumGeneLengths_sorted.txt > expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep3.txt

# Fix the names of the Sinorhizobium genes in the TPM table
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep1.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep2.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep3.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep1.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep2.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep3.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep1.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep2.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep3.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep1.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep2.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep3.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep1.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep2.txt
sed -i 's/SM/sm/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep3.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep1.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep2.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep3.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep1.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep2.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep3.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep1.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep2.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep3.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep1.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep2.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep3.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep1.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep2.txt
sed -i 's/sm_b/smb/g' expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep3.txt

# Combine TPM values into one file per species
mkdir _finalOutput/
paste expressionAnalysis/TPMvaluesMedicago/zoneI_rep1.txt expressionAnalysis/TPMvaluesMedicago/zoneI_rep2.txt expressionAnalysis/TPMvaluesMedicago/zoneI_rep3.txt expressionAnalysis/TPMvaluesMedicago/zoneIId_rep1.txt expressionAnalysis/TPMvaluesMedicago/zoneIId_rep2.txt expressionAnalysis/TPMvaluesMedicago/zoneIId_rep3.txt expressionAnalysis/TPMvaluesMedicago/zoneIIp_rep1.txt expressionAnalysis/TPMvaluesMedicago/zoneIIp_rep2.txt expressionAnalysis/TPMvaluesMedicago/zoneIIp_rep3.txt expressionAnalysis/TPMvaluesMedicago/zoneIZ_rep1.txt expressionAnalysis/TPMvaluesMedicago/zoneIZ_rep2.txt expressionAnalysis/TPMvaluesMedicago/zoneIZ_rep3.txt expressionAnalysis/TPMvaluesMedicago/zoneIII_rep1.txt expressionAnalysis/TPMvaluesMedicago/zoneIII_rep2.txt expressionAnalysis/TPMvaluesMedicago/zoneIII_rep3.txt | cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 > _finalOutput/Medicago_TPM_values.txt
paste expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep1.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep2.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneI_rep3.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep1.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep2.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIId_rep3.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep1.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep2.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIIp_rep3.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep1.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep2.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIZ_rep3.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep1.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep2.txt expressionAnalysis/TPMvaluesSinorhizobium/zoneIII_rep3.txt | cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 > _finalOutput/Sinorhizobium_TPM_values.txt
echo -e "Gene\tzoneI_1\tzoneI_2\tzoneI_3\tzoneIId_1\tzoneIId_2\tzoneIId_3\tzoneIIp_1\tzoneIIp_2\tzoneIIp_3\tzoneIZ_1\tzoneIZ_2\tzoneIZ_3\tzoneIII_1\tzoneIII_2\tzoneIII_3\n$(cat _finalOutput/Medicago_TPM_values.txt)" > _finalOutput/Medicago_TPM_values.txt
echo -e "Gene\tzoneI_1\tzoneI_2\tzoneI_3\tzoneIId_1\tzoneIId_2\tzoneIId_3\tzoneIIp_1\tzoneIIp_2\tzoneIIp_3\tzoneIZ_1\tzoneIZ_2\tzoneIZ_3\tzoneIII_1\tzoneIII_2\tzoneIII_3\n$(cat _finalOutput/Sinorhizobium_TPM_values.txt)" > _finalOutput/Sinorhizobium_TPM_values.txt

# Calculate average TPM values
perl Scripts/averageTPM.pl _finalOutput/Medicago_TPM_values.txt > _finalOutput/Medicago_TPM_values_average.txt
perl Scripts/averageTPM.pl _finalOutput/Sinorhizobium_TPM_values.txt > _finalOutput/Sinorhizobium_TPM_values_average.txt

# Calculate TPM values across zones
perl Scripts/percentTPM.pl _finalOutput/Medicago_TPM_values_average.txt > _finalOutput/Medicago_TPM_values_percent.txt
perl Scripts/percentTPM.pl _finalOutput/Sinorhizobium_TPM_values_average.txt > _finalOutput/Sinorhizobium_TPM_values_percent.txt


