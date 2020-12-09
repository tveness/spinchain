set terminal png
set output "hist_mx.png"

binwidth = 0.005
set boxwidth binwidth
bin(x,width)=width*floor(x/width) + width/2.0 

set xlabel "M_x"
set ylabel "Count"

plot "hist_mc.dat" u (bin($2,binwidth)):(1.0) smooth freq with boxes title "Monte Carlo", "hist_dyn.dat" u (bin($2,binwidth)):(1.0) smooth freq with boxes title "Dynamics"

binwidth = 0.01
set boxwidth binwidth
set xlabel "M_y"
set ylabel "Count"

set output "hist_my.png"
plot "hist_mc.dat" u (bin($3,binwidth)):(1.0) smooth freq with boxes title "Monte Carlo", "hist_dyn.dat" u (bin($3,binwidth)):(1.0) smooth freq with boxes title "Dynamics"


set xlabel "M_z"
set ylabel "Count"

set output "hist_mz.png"
plot "hist_mc.dat" u (bin($4,binwidth)):(1.0) smooth freq with boxes title "Monte Carlo", "hist_dyn.dat" u (bin($4,binwidth)):(1.0) smooth freq with boxes title "Dynamics"

set xlabel "e"
set ylabel "Count"

set output "hist_e.png"
plot "hist_mc.dat" u (bin($5,binwidth)):(1.0) smooth freq with boxes title "Monte Carlo", "hist_dyn.dat" u (bin($5,binwidth)):(1.0) smooth freq with boxes title "Dynamics"
