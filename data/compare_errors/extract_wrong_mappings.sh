cp ../k8-Linux .
cp ../paftools.js .
echo "Extracting wrong mappings for minichain..."
./k8-Linux paftools.js mapeval -Q0 minichain/Linear.gaf | sed -n "/^E/p" -  > minichain/wrong_mappings_Linear.gaf
./k8-Linux paftools.js mapeval minichain/Linear.gaf | awk '{ print $6 }' | tail -n 1 - > minichain/count_mappings_Linear.gaf
./k8-Linux paftools.js mapeval -Q0 minichain/10H.gaf | sed -n "/^E/p" -  > minichain/wrong_mappings_10H.gaf
./k8-Linux paftools.js mapeval minichain/10H.gaf | awk '{ print $6 }' | tail -n 1 - > minichain/count_mappings_10H.gaf
./k8-Linux paftools.js mapeval -Q0 minichain/20H.gaf | sed -n "/^E/p" -  > minichain/wrong_mappings_20H.gaf
./k8-Linux paftools.js mapeval minichain/20H.gaf | awk '{ print $6 }' | tail -n 1 - > minichain/count_mappings_20H.gaf
./k8-Linux paftools.js mapeval -Q0 minichain/40H.gaf | sed -n "/^E/p" -  >minichain/wrong_mappings_40H.gaf
./k8-Linux paftools.js mapeval minichain/40H.gaf | awk '{ print $6 }' | tail -n 1 - > minichain/count_mappings_40H.gaf
./k8-Linux paftools.js mapeval -Q0 minichain/60H.gaf | sed -n "/^E/p" -  > minichain/wrong_mappings_60H.gaf
./k8-Linux paftools.js mapeval minichain/60H.gaf | awk '{ print $6 }' | tail -n 1 - > minichain/count_mappings_60H.gaf
./k8-Linux paftools.js mapeval -Q0 minichain/80H.gaf | sed -n "/^E/p" -  > minichain/wrong_mappings_80H.gaf
./k8-Linux paftools.js mapeval minichain/80H.gaf | awk '{ print $6 }' | tail -n 1 - > minichain/count_mappings_80H.gaf
./k8-Linux paftools.js mapeval -Q0 minichain/95H.gaf | sed -n "/^E/p" -  >minichain/wrong_mappings_95H.gaf
./k8-Linux paftools.js mapeval minichain/95H.gaf | awk '{ print $6 }' | tail -n 1 - > minichain/count_mappings_95H.gaf
echo "Extracted wrong mappings for minichain!"

echo "Extracting wrong mappings for minigraph..."
./k8-Linux paftools.js mapeval -Q0 minigraph/Linear.gaf | sed -n "/^E/p" -  > minigraph/wrong_mappings_Linear.gaf
./k8-Linux paftools.js mapeval minigraph/Linear.gaf | awk '{ print $6 }' | tail -n 1 - > minigraph/count_mappings_Linear.gaf
./k8-Linux paftools.js mapeval -Q0 minigraph/10H.gaf | sed -n "/^E/p" -  > minigraph/wrong_mappings_10H.gaf
./k8-Linux paftools.js mapeval minigraph/10H.gaf | awk '{ print $6 }' | tail -n 1 - > minigraph/count_mappings_10H.gaf
./k8-Linux paftools.js mapeval -Q0 minigraph/20H.gaf | sed -n "/^E/p" -  > minigraph/wrong_mappings_20H.gaf
./k8-Linux paftools.js mapeval minigraph/20H.gaf | awk '{ print $6 }' | tail -n 1 - > minigraph/count_mappings_20H.gaf
./k8-Linux paftools.js mapeval -Q0 minigraph/40H.gaf | sed -n "/^E/p" -  > minigraph/wrong_mappings_40H.gaf
./k8-Linux paftools.js mapeval minigraph/40H.gaf | awk '{ print $6 }' | tail -n 1 - > minigraph/count_mappings_40H.gaf
./k8-Linux paftools.js mapeval -Q0 minigraph/60H.gaf | sed -n "/^E/p" -  > minigraph/wrong_mappings_60H.gaf
./k8-Linux paftools.js mapeval minigraph/60H.gaf | awk '{ print $6 }' | tail -n 1 - > minigraph/count_mappings_60H.gaf
./k8-Linux paftools.js mapeval -Q0 minigraph/80H.gaf | sed -n "/^E/p" -  > minigraph/wrong_mappings_80H.gaf
./k8-Linux paftools.js mapeval minigraph/80H.gaf | awk '{ print $6 }' | tail -n 1 - > minigraph/count_mappings_80H.gaf
./k8-Linux paftools.js mapeval -Q0 minigraph/95H.gaf | sed -n "/^E/p" -  > minigraph/wrong_mappings_95H.gaf
./k8-Linux paftools.js mapeval minigraph/95H.gaf | awk '{ print $6 }' | tail -n 1 - > minigraph/count_mappings_95H.gaf
echo "Extracted wrong mappings for minigraph!"


echo "Extracting wrong mappings for GraphAligner..."
./k8-Linux paftools.js mapeval -Q0 GraphAligner/Linear.gaf | sed -n "/^E/p" -  > GraphAligner/wrong_mappings_Linear.gaf
./k8-Linux paftools.js mapeval GraphAligner/Linear.gaf | awk '{ print $6 }' | tail -n 1 - > GraphAligner/count_mappings_Linear.gaf
./k8-Linux paftools.js mapeval -Q0 GraphAligner/10H.gaf | sed -n "/^E/p" -  > GraphAligner/wrong_mappings_10H.gaf
./k8-Linux paftools.js mapeval GraphAligner/10H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphAligner/count_mappings_10H.gaf
./k8-Linux paftools.js mapeval -Q0 GraphAligner/20H.gaf | sed -n "/^E/p" -  > GraphAligner/wrong_mappings_20H.gaf
./k8-Linux paftools.js mapeval GraphAligner/20H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphAligner/count_mappings_20H.gaf
./k8-Linux paftools.js mapeval -Q0 GraphAligner/40H.gaf | sed -n "/^E/p" -  > GraphAligner/wrong_mappings_40H.gaf
./k8-Linux paftools.js mapeval GraphAligner/40H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphAligner/count_mappings_40H.gaf
./k8-Linux paftools.js mapeval -Q0 GraphAligner/60H.gaf | sed -n "/^E/p" -  > GraphAligner/wrong_mappings_60H.gaf
./k8-Linux paftools.js mapeval GraphAligner/60H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphAligner/count_mappings_60H.gaf
./k8-Linux paftools.js mapeval -Q0 GraphAligner/80H.gaf | sed -n "/^E/p" -  > GraphAligner/wrong_mappings_80H.gaf
./k8-Linux paftools.js mapeval GraphAligner/80H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphAligner/count_mappings_80H.gaf
./k8-Linux paftools.js mapeval -Q0 GraphAligner/95H.gaf | sed -n "/^E/p" -  > GraphAligner/wrong_mappings_95H.gaf
./k8-Linux paftools.js mapeval GraphAligner/95H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphAligner/count_mappings_95H.gaf
echo "Extracted wrong mappings for GraphAligner!"


echo "Extracting wrong mappings for GraphChainer..."
./k8-Linux paftools.js mapeval -Q60 GraphChainer/Linear.gaf | sed -n "/^E/p" -  > GraphChainer/wrong_mappings_Linear.gaf
./k8-Linux paftools.js mapeval GraphChainer/Linear.gaf | awk '{ print $6 }' | tail -n 1 - > GraphChainer/count_mappings_Linear.gaf
./k8-Linux paftools.js mapeval -Q60 GraphChainer/10H.gaf | sed -n "/^E/p" -  > GraphChainer/wrong_mappings_10H.gaf
./k8-Linux paftools.js mapeval GraphChainer/10H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphChainer/count_mappings_10H.gaf
./k8-Linux paftools.js mapeval -Q60 GraphChainer/20H.gaf | sed -n "/^E/p" -  > GraphChainer/wrong_mappings_20H.gaf
./k8-Linux paftools.js mapeval GraphChainer/20H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphChainer/count_mappings_20H.gaf
./k8-Linux paftools.js mapeval -Q60 GraphChainer/40H.gaf | sed -n "/^E/p" -  > GraphChainer/wrong_mappings_40H.gaf
./k8-Linux paftools.js mapeval GraphChainer/40H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphChainer/count_mappings_40H.gaf
./k8-Linux paftools.js mapeval -Q60 GraphChainer/60H.gaf | sed -n "/^E/p" -  > GraphChainer/wrong_mappings_60H.gaf
./k8-Linux paftools.js mapeval GraphChainer/60H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphChainer/count_mappings_60H.gaf
./k8-Linux paftools.js mapeval -Q60 GraphChainer/80H.gaf | sed -n "/^E/p" -  > GraphChainer/wrong_mappings_80H.gaf
./k8-Linux paftools.js mapeval GraphChainer/80H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphChainer/count_mappings_80H.gaf
./k8-Linux paftools.js mapeval -Q60 GraphChainer/95H.gaf | sed -n "/^E/p" -  > GraphChainer/wrong_mappings_95H.gaf
./k8-Linux paftools.js mapeval GraphChainer/95H.gaf | awk '{ print $6 }' | tail -n 1 - > GraphChainer/count_mappings_95H.gaf
echo "Extracted wrong mappings for GraphChainer!"