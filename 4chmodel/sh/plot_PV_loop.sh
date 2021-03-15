# gnuplot -e "heart='82'" plot_PV_loop.sh 

set terminal png background rgb 'white' fontscale 1.8 size 900,700 font 'Helvetica'

# change a color of border.
set border lw 3 lc rgb "black"

# change text colors of  tics
set xtics textcolor rgb "black"
set ytics textcolor rgb "black"

# change text colors of labels
set xlabel "X" textcolor rgb "black"
set ylabel "Y" textcolor rgb "black"

# change a text color of key
set key textcolor rgb "black"

# change a text color of key
set title textcolor rgb "black"
set key left top

set xtics 30,30,270
set ytics 0,60,180

set xrange[30:280]
set yrange[-20:180]

set title 'PV loop'
set xlabel 'Volume (mL)'
set ylabel 'Pressure (mmHg)'
set key on
set grid

n=0
do for [ii=1:800:20] {
    n=n+1
    set output sprintf('/media/crg17/Seagate Expansion Drive/h_case'.heart.'/simulations/cohort/'.heart.'HC_wk350_tpeak100/pvloops/pvloops%03.0f.png',n)
    plot '/media/crg17/Seagate Expansion Drive/h_case'.heart.'/simulations/cohort/'.heart.'HC_wk350_tpeak100/cav.LV.csv' using 3:2 every ::1::ii title 'Left ventricle' with lines lw 3 lt rgb "blue", \
   		 '/media/crg17/Seagate Expansion Drive/h_case'.heart.'/simulations/cohort/'.heart.'HC_wk350_tpeak100/cav.RV.csv' using 3:2 every ::1::ii title 'Right ventricle' with lines lw 3 lc rgb "red" 
}
