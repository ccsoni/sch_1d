#!/bin/bash

# run_model model_name N hbar
function run_model() {
    ./sch_1d -T 3.14 -t 0.01 -m $1 -N $2 -d $3
    python3 wavefunc_anim.py "$1/$1_???.dat" 3
    ./convert_png_to_animation_gif.sh "$1/$1_???.png" $1/RePsi.gif
    python3 wavefunc_anim.py "$1/$1_???.dat" 4
    ./convert_png_to_animation_gif.sh "$1/$1_???.png" $1/ImPsi.gif
    python3 wavefunc_anim.py "$1/$1_???.dat" 1
    ./convert_png_to_animation_gif.sh "$1/$1_???.png" Prob.gif
    python3 DF_anim.py "$1/$1_DF_???.dat" 1
    ./convert_png_to_animation_gif.sh "$1/$1_DF_???.png" DF.gif
}

#run_model harmonic_N128_d0.05_5pnt 128 0.05
#run_model harmonic_N256_d0.05_5pnt 256 0.05
run_model harmonic_N256_d0.01_5pnt 256 0.01
