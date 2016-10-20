#!/bin/bash

mkdir roborobo3
cd roborobo3

cat ~/.ssh/id_rsa.pub
#add key to github

git clone git@github.com:inakiFernandez/roborobo3.git .

./makefile-manager -i TemplateWander
./makefile-manager -i TemplateBoids
./makefile-manager -i TemplateRandomwalk
./makefile-manager -i TemplateMedea
./makefile-manager -a Original

#install SDL2
#install boost

#FIXED SDL_Init
#install export
#xset -display ${DISPLAY} dpms force on
#export DISPLAY=:1.0
#sudo xinit -- $DISPLAY

make clean

make