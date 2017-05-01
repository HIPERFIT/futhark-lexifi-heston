(cat style.gnuplot
 echo 'set output "speedup.svg"'
 echo 'plot "plot_data" with linespoints ls 1 title "Futhark on K40"'
) | gnuplot
