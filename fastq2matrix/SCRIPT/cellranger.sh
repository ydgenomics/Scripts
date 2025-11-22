/opt/cellranger-9.0.1/bin/cellranger mkref \
--genome cotton \
--fasta /data/work/Reference/Gossypium/genome.fa \
--genes /data/work/Reference/Gossypium/genes.gtf \
--localcores 8 \
--localmem 64 \
--output-dir cotton


/opt/cellranger-9.0.1/bin/cellranger count \
