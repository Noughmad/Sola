set terminal eps color solid
set pm3d map
set palette model RGB
set palette model RGB defined (-1 "blue", 0 0.75 0.85 0, 1 "red")

unset tics
unset border

set xzeroaxis
HalfWidth=264

set output "g_snap_gauss.eps"
splot "./Data/cross_1.dat" using 1:($2-HalfWidth):3 matrix

set output "g_snap_split.eps"
splot "./Data/cross_500.dat" using 1:($2-HalfWidth):3 matrix

set output "g_snap_inter.eps"
splot "./Data/cross_2300.dat" using 1:($2-HalfWidth):3 matrix