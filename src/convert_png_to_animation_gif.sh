#!/bin/bash

fname=$1

base="${fname%.*}"
prefix="${base//\?}"

#magick -delay 30 -loop 0  $1 ${prefix}.gif
magick -delay 30 -loop 0  $1 $2
