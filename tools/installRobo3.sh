#!/bin/bash

mkdir roborobo3
cd roborobo3

cat ~/.ssh/id_rsa.pub
#By hand: add key to github

git clone git@github.com:inakiFernandez/roborobo3.git .

./makefile-manager -i TemplateWander
./makefile-manager -i TemplateBoids
./makefile-manager -i TemplateRandomwalk
./makefile-manager -i TemplateMedea
./makefile-manager -a Original

#install SDL2
#install boost

#FIXED SDL_Init by hand. This is already fixed on code
#install export
#xset -display ${DISPLAY} dpms force on
#export DISPLAY=:1.0
#sudo xinit -- $DISPLAY

make clean

make

#######################
#For data analysis(plots, videos, etc)
#install avconv or ffmpeg
#install realpath
#install virtualenv
#install python3 (.4 at least)
#install imagick

#Set up virtualenv
cd ..
virtualenv --system-site-packages -p python3 virtualenv-roborobo

$currentdir/virtualenv-roborobo/bin/pip3 install scipy matplotlib numpy brewer2mpl
$currentdir/virtualenv-roborobo/bin/pip3 install --upgrade scipy matplotlib numpy brewer2mpl

currentDir=`realpath .`
#if tcsh: set up PATH variables to use virtualenv by default in my session
echo "setenv PYTHONPATH  $currentDir/roborobo3/tools:" >> ~/.tcshrc
echo "setenv PYTHONPATH  $currentDir/virtualenv-roborobo:$PYTHONPATH" >> ~/.tcshrc
echo "setenv PATH $currentDir/virtualenv-roborobo/bin:$PATH" >> ~/.tcshrc
echo "alias python3 $currentDir/virtualenv-roborobo/bin/python3" >> ~/.tcshrc

source ~/.tcshrc

