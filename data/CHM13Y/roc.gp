set t po eps enh co so "Helvetica,26"

set style line 1 lt 1 pt 1 lc rgb "#e41a1c" lw 2;
set style line 2 lt 1 pt 2 lc rgb "#377eb8" lw 2;
set style line 3 lt 1 pt 3 lc rgb "#4daf4a" lw 2;
set style line 4 lt 1 pt 4 lc rgb "#984ea3" lw 2;
set style line 5 lt 1 pt 6 lc rgb "#ff7f00" lw 2;
set style line 6 lt 1 pt 8 lc rgb "#f781bf" lw 2;

set out "roc-color.eps"

set pointsize 2.0

set xlab "Error rate of mapped PacBio reads"
set ylab "Fraction of mapped reads" off +1.8
set ytics 0.02
set yran [0.9:1.02]

set log x
set format x "10^{%L}"
set key bot right
plot "<./eval2roc.pl minichain.eval" u 2:3 t "minichain" w lp ls 1, \
     "<./eval2roc.pl minigraph.eval" u 2:3 t "minigraph" w lp ls 2, \
     "<./eval2roc.pl GraphAligner.eval" u 2:3 t "GraphAligner" w lp ls 3, \
     "<./eval2roc.pl GraphChainer.eval" u 2:3 t "GraphChainer" w lp ls 4

unset label
