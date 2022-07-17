conda activate minichain_data
echo "Evaluating results for graphs..."
./k8 paftools.js mapeval Results/GraphChainer/CHM13Y_95H_5%_stable.gaf > CHM13Y_95H/GraphChainer_95H.eval
./k8 paftools.js mapeval Results/GraphChainer/CHM13Y_10H_5%_stable.gaf > CHM13Y_10H/GraphChainer_10H.eval

./k8 paftools.js mapeval Results/GraphAligner/CHM13Y_95H_5%_stable.gaf > CHM13Y_95H/GraphAligner_95H.eval
./k8 paftools.js mapeval Results/GraphAligner/CHM13Y_10H_5%_stable.gaf > CHM13Y_10H/GraphAligner_10H.eval

./k8 paftools.js mapeval Results/minigraph/CHM13Y_95H_5%.gaf > CHM13Y_95H/minigraph_95H.eval
./k8 paftools.js mapeval Results/minigraph/CHM13Y_10H_5%.gaf > CHM13Y_10H/minigraph_10H.eval

./k8 paftools.js mapeval Results/minichain/CHM13Y_95H_5%.gaf > CHM13Y_95H/minichain_95H.eval
./k8 paftools.js mapeval Results/minichain/CHM13Y_10H_5%.gaf > CHM13Y_10H/minichain_10H.eval

echo "Evaluating results for linear genome..."
./k8 paftools.js mapeval Results/minigraph/CHM13Y_L_5%.gaf > CHM13Y/minigraph.eval
./k8 paftools.js mapeval Results/minichain/CHM13Y_L_5%.gaf > CHM13Y/minichain.eval
./k8 paftools.js mapeval Results/minimap2/CHM13Y_L_5%.gaf > CHM13Y/minimap2.eval

gnuplot CHM13Y_95H/roc.gp
convert -density 1200 CHM13Y_95H/roc-color.eps -resize 1920x1080 95H.jpg
gnuplot CHM13Y_10H/roc.gp
convert -density 1200 CHM13Y_10H/roc-color.eps -resize 1920x1080 10H.jpg
gnuplot CHM13Y/roc.gp
convert -density 1200 CHM13Y/roc-color.eps -resize 1920x1080 Linear.jpg