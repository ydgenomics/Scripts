samtools view -u -f 1 -F 12 lib_002.sorted.md.bam > lib_002_map_map.bam

# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 lib_002.sorted.md.bam > lib_002_unmap_map.bam

# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 lib_002.sorted.md.bam > lib_002_map_unmap.bam

# R1 & R2 unmapped
samtools view -u -f 12 -F 256 lib_002.sorted.md.bam > lib_002_unmap_unmap.bam

samtools merge -u lib_002_unmapped.bam lib_002_unmap_map.bam lib_002_map_unmap.bam lib_002_unmap_unmap.bam