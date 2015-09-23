
file="../../output/bilayer_stack_evaporate/timeseries"
file2="../../output/bilayer_stack_evaporate_relax/timeseries"

set logscale x
set format x "10^{%L}"
set xr [1e3:*]
set xlabel "MC cycles"
set ylabel "E/monomer (kcal/mol)"

plot file u (column(1)):(column(2)) with lines lw 8 lt 1 title "backbone", \
file u (column(1)):(column(3)) with lines lw 8 lt 2 title "sidechain", \
file u (column(1)):(column(4)+column(8)) with lines lw 8 lt 3 title "phen, same", \
file u (column(1)):(column(6)+column(7)+column(10)+column(11)) with lines lw 8 lt 4 title "charged, same", \
file u (column(1)):(column(12)) with lines lw 8 lt 5 title "phen, cross", \
file u (column(1)):(column(14)+column(15)) with lines lw 8 lt 6 title "charged, cross", \
file u (column(1)):(column(2)+column(3)+column(4)+column(5)+column(6)+column(7)+column(8)+column(9)+column(10)+column(11)+column(12)+column(13)+column(14)+column(15)) with lines lw 8 lt 7 title "total", \
file2 u (column(1) + 100000):(column(2)) with lines lw 8 lt 1 title "", \
file2 u (column(1) + 100000):(column(3)) with lines lw 8 lt 2 title "", \
file2 u (column(1) + 100000):(column(4)+column(8)) with lines lw 8 lt 3 title "", \
file2 u (column(1) + 100000):(column(6)+column(7)+column(10)+column(11)) with lines lw 8 lt 4 title "", \
file2 u (column(1) + 100000):(column(12)) with lines lw 8 lt 5 title "", \
file2 u (column(1) + 100000):(column(14)+column(15)) with lines lw 8 lt 6 title "", \
file2 u (column(1) + 100000):(column(2)+column(3)+column(4)+column(5)+column(6)+column(7)+column(8)+column(9)+column(10)+column(11)+column(12)+column(13)+column(14)+column(15)) with lines lw 8 lt 7 title ""

reset

