#E_van_MNase_1.fq.gz	E_van_MNase_2.fq.gz	E_van_H3K4me3_1.fq.gz	E_van_H3K4me3_2.fq.gz
mkdir trim
mkdir MNase
name="E_van_MNase"
fastp -i ${name}_1.fq.gz -I ${name}_2.fq.gz -o trim/${name}_1.fastp.fq.gz -O trim/${name}_2.fastp.fq.gz --cut_tail --cut_window_size 3 --cut_mean_quality 20 -q 20 -l 20 -w 32 -h trim/"$name"_fastp.html -j trim/"$name"_fastp.json
bowtie2 -x Eup.bowtie2 -1 trim/${name}_1.fastp.fq.gz -2 trim/${name}_2.fastp.fq.gz -S MNase/${name}_mapped.sam -p 32 --no-mixed --no-discordant
grep "@" MNase/${name}_mapped.sam > MNase/${name}_header.sam
grep "AS:" MNase/${name}_mapped.sam | grep -v "XS:"  > MNase/${name}_mapped1.sam
cat MNase/${name}_header.sam MNase/${name}_mapped1.sam > MNase/${name}_uniq_mapped.sam
rm MNase/${name}_mapped.sam
rm MNase/${name}_mapped1.sam
rm MNase/${name}_header.sam
samtools view -@ 32 -b -S MNase/${name}_uniq_mapped.sam -o MNase/${name}_uniq_mapped.bam
rm MNase/${name}_uniq_mapped.sam
sambamba sort -t 32 MNase/${name}_uniq_mapped.bam -o MNase/${name}_sorted.bam && \
sambamba markdup -r -t 10 --tmpdir=./MNase MNase/${name}_sorted.bam MNase/${name}_rmdup.bam
sambamba view -h MNase/${name}_rmdup.bam | awk '/^@/ || ($9 >= 90 && $9 <= 200) || ($9 <= -90 && $9 >= -200)' | samtools view -Sb - > MNase/$name.90-200.bam
samtools view MNase/$name.90-200.bam > rm MNase/$name.90-200.sam
samtools index MNase/$name.90-200.bam

mkdir H3K4me3
name="E_van_H3K4me3"
fastp -i ${name}_1.fq.gz -I ${name}_2.fq.gz -o trim/${name}_1.fastp.fq.gz -O trim/${name}_2.fastp.fq.gz --cut_tail --cut_window_size 3 --cut_mean_quality 20 -q 20 -l 20 -w 32 -h trim/"$b"_fastp.html -j trim/"$b"_fastp.json
bowtie2 -x bowtie2 -1 trim/${name}_1.fastp.fq.gz -2 trim/${name}_2.fastp.fq.gz -S H3K4me3/${name}_mapped.sam -p 32 --no-mixed --no-discordant
grep "@" H3K4me3/${name}_mapped.sam > H3K4me3/${name}_header.sam
grep "AS:" H3K4me3/${name}_mapped.sam | grep -v "XS:"  > H3K4me3/${name}_mapped1.sam
cat H3K4me3/${name}_header.sam H3K4me3/${name}_mapped1.sam > H3K4me3/${name}_uniq_mapped.sam
rm H3K4me3/${name}_mapped.sam
rm H3K4me3/${name}_mapped1.sam
rm H3K4me3/${name}_header.sam
samtools view -@ 32 -b -S H3K4me3/${name}_uniq_mapped.sam -o H3K4me3/${name}_uniq_mapped.bam
rm H3K4me3/${name}_uniq_mapped.sam
sambamba sort -t 32 H3K4me3/${name}_uniq_mapped.bam -o H3K4me3/${name}_sorted.bam && \
sambamba markdup -r -t 10 --tmpdir=./H3K4me3 H3K4me3/${name}_sorted.bam H3K4me3/${name}_rmdup.bam
sambamba view -h H3K4me3/${name}_rmdup.bam | awk '/^@/ || ($9 >= 100 && $9 <= 250) || ($9 <= -100 && $9 >= -250)' | samtools view -Sb - > H3K4me3/$name.100-250.bam
samtools index H3K4me3/$name.100-250.bam
bamCoverage --ignoreDuplicates --outFileFormat bigwig --skipNAs --bam H3K4me3/$name.100-250.bam --outFileName H3K4me3/$name.100-250.bigwig --normalizeUsing RPKM --minMappingQuality 30 --binSize 2 -p 32
