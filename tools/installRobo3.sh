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

make clean

make