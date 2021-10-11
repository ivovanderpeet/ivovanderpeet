reset
## please specify here where your data file is.
data = "./output.dat"

## plot temperature profile
set terminal png font arial 12 size 640, 480
set output "temperature.png"

set xlabel "X axis"
set ylabel "Y axis"
set zlabel "T [K]"
set key at screen 1.0, 0.9
set contour
##set cntrparam level 20
set cntrparam level incr 270,20,370
## set cntrparam level disc 300, 350, 370

splot data u 1:2:6 w l

set output

set output "temperature_color.png"
set pm3d at b
replot
set output

reset
## plot velocity vector map
set terminal png font arial 12 size 640, 480
set output "velocity_map.png"

set xlabel "X axis"
set ylabel "Y axis"

plot data u 1:2:($3*2):($4*2) every 1:2 w vec head filled

set output

reset
## plot velocity profile at output
set terminal png font arial 12 size 640, 480
set output "velocity_out.png"

set xlabel "Y axis"
set ylabel "u [m/s]"
set y2label "T [K]"
set ytics nomirror
set y2tics nomirror

plot data u 2:3 every :::48::48 w l lw 2.5 t"velocity",\
     data u 2:6 every :::48::48 w l lw 2.5 axes x1y2 t "temperature"

set output






