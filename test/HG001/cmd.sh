#bwa mem /root/tools/git_repo/pbAmpliconAnalysis_HLA/database/ref/HLA.fasta HG001.PE150.R1.fastq HG001.PE150.R2.fastq >HG001.sam
samtools view -b -o HG001.bam HG001.sam
samtools sort -o HG001.sort.bam HG001.bam
samtools index HG001.sort.bam
