EXE="`pwd`/vrvica/build/vrvica"

for kot in 0,9 1,0 1,1 1,2
do
    mkdir "g_animacija_${kot}"
    cd "g_animacija_${kot}"
    "${EXE}" 100 ${kot} 1000
    ffmpeg -r 25 -i g_frame_%04d.jpg g_animacija_${kot}.avi
    cp g_animacija_${kot}.avi ..
    cd ..
done
