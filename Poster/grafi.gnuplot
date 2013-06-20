set terminal png size 2400, 1200
set pm3d map interpolate 10,1
set palette model RGB
# set palette model RGB defined (-1 "blue", 0 0.75 0.85 0, 1 "red")
set palette model RGB defined (-1 "blue", 0 0.8 0.8 0.8, 1 "red")

unset tics
unset border

# set xzeroaxis
HalfWidth=136

set output "g_snap_gauss.png"
splot "./Data/cross_0.dat" using 1:($2-HalfWidth):3 matrix

set output "g_snap_split.png"
splot "./Data/cross_500.dat" using 1:($2-HalfWidth):3 matrix

set output "g_snap_inter.png"
splot "./Data/cross_1400.dat" using 1:($2-HalfWidth):3 matrix