echo "Mapping reads with minimap2..."
./minimap2/minimap2 -t32 -cx map-pb CHM13Y.fa Reads/CHM13Y_reads_5%.fa  > Results/minimap2/CHM13Y_L_5%.gaf
echo "Mapping reads with minigraph..."
./minigraph/minigraph -t32 -cx lr CHM13Y.fa Reads/CHM13_reads_5%.fa > Results/minigraph/CHM13Y_L_5%.gaf
echo "Mapping reads with minichain..."
../minichain -t32 -cx lr CHM13Y.fa Reads/CHM13_reads_5%.fa  -y1 > Results/minichain/CHM13Y_L_5%.gaf
echo "Mapping reads to Linear CHM13Y genome."