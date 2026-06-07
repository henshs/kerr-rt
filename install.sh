# Radiation from Keplerian disk simaluations.

echo "**************************************************************"
echo "Preparing initial setup"
echo "**************************************************************"

mkdir simulations
mkdir simulations/bin
mkdir simulations/init
mkdir simulations/data
mkdir simulations/data/sim_data
mkdir simulations/data/image_data
mkdir simulations/png

echo "**************************************************************"
echo "keplerian_disk_firewall - collecting rays from the whole disk"
echo "**************************************************************"

cd kdfirewall
make
cp bin/kdfirewall ../simulations/bin
cd ..

echo "**************************************************************"
echo "keplerian_disk_image -  keplerian disk image builder"
echo "**************************************************************"

cd kdimage
make
cp bin/kdimage ../simulations/bin
cd ..

echo "**************************************************************"
echo "kerr rms builder"
echo "**************************************************************"

cd kerr_rms
make
cp bin/kerrrms ../simulations/bin
cd ..


echo "**************************************************************"
echo "created a sample run script"
echo "**************************************************************"

sample_run="simulations/run.sh"

cat << 'EOF' > "$sample_run"
#!/bin/bash
a=0.1
tho=85.0

uo=1e-3
rmax=20
alphamin=-25.0
alphamax=25.0
betamin=-25.0
betamax=25.0

w=1200
h=1200

font=$((w / 60))

num=1.000001

mkdir -p init ./data/sim_data ./data/image_data png

rms=$(./bin/kerrrms $a 1)
rmin=$(echo "$num * $rms" | bc -l)

initdir="./init/"
finit="${initdir}/all_photons_init_bh_a${a}_tho${tho}.tab"

echo "a $a" > "$finit"
echo "obs $uo $tho" >> "$finit"
echo "sky $alphamin $alphamax $betamin $betamax" >> "$finit"
echo "skysize $w $h" >> "$finit"
echo "kdisk $rms $rmax" >> "$finit"
echo "opacity $opacity" >> "$finit"

datadir="./data/"
fdata="${datadir}sim_data/data_a_${a}_tho_${tho}.dat"
fimagedata="${datadir}image_data/image_a_${a}_tho_${tho}.dat"
pngdir="./png/"
pngfile="${pngdir}disk_a_${a}_tho_${tho}.png"
    
echo "Generating Disk Data"
echo "$w $h $alphamin $alphamax $betamin $betamax" > "$fdata"
./bin/kdfirewall "$finit" >> "$fdata"
            
echo "Generating Disk Image"
./bin/kdimage "$fdata" > "$fimagedata"

# GNUPLOT SCRIPT
plot="plot_script.plt"

cat << GNUPLOT_EOF > "$plot"
#!/usr/bin/gnuplot
set term pngcairo size 800,800 font "Times,21"
set output "$pngfile"
unset key
set style line 1 lc rgb 'gray' lt 1 lw 1
set size square
set palette rgb 21,22,23
set xrange [$alphamin:$alphamax]
set yrange [$betamin:$betamax]
set xlabel "{/Symbol a}"
set ylabel "{/Symbol b}" rotate by 0
plot "$fimagedata" u 1:2:3 w image
GNUPLOT_EOF

gnuplot "$plot"
rm "$plot"
echo "Disk Image generated"
EOF

chmod +x "$sample_run"



