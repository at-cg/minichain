set t po eps enh co so "Helvetica,26"

set style line 1 lt 1 pt 1 lc rgb "#e41a1c" lw 2;
set style line 2 lt 1 pt 2 lc rgb "#377eb8" lw 2;
set style line 3 lt 1 pt 3 lc rgb "#4daf4a" lw 2;
set style line 4 lt 1 pt 4 lc rgb "#984ea3" lw 2;
set style line 5 lt 1 pt 6 lc rgb "#ff7f00" lw 2;
set style line 6 lt 1 pt 8 lc rgb "#f781bf" lw 2;

set out "95H.eps"

set pointsize 2.0

set key font ",32"
set key at graph 1.01, 0.025

set xlab "{/*1.5 Incorrectly aligned reads}"
set ylab "{/*1.5 Aligned reads}" off -0.5
set ytics 2
#set xtics 10
set yran [90:100.05]

set xtics font ", 35"
set ytics font ", 35"
set rmargin 3.5

set xrange [:10]
set log x
#set format x "10^{%L}"
set format x "%.01f%%"
set format y "%.0f%%"
set mxtics 10
set key bot right
plot "<./eval2roc.pl minichain_95H.eval" u 2:3 t "Minichain" w lp ls 1, \
     "<./eval2roc.pl minigraph_95H.eval" u 2:3 t "Minigraph" w lp ls 2, \
     "<./eval2roc.pl GraphAligner_95H.eval" u 2:3 t "GraphAligner" w lp ls 3, \
     "<./eval2roc.pl GraphChainer_95H.eval" u 2:3 t "GraphChainer" w lp ls 4,
unset label
