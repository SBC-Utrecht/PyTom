#!/bin/bash

git clone git://git.code.sf.net/p/pytom/gitrepo ./pytom

rm -rf pytom/.git
rm -rf pytom/*/.git
rm -rf pytom/*/*/.git
rm -rf pytom/*/*/*/.git
rm -rf pytom/*/*/*/*/.git
rm -rf pytom/*/*/*/*/*/.git
rm -rf pytom/*/*/*/*/*/*/.git
rm -rf pytom/*/*/*/*/*/*/*/.git
rm -rf pytom/*/*/*/*/*/*/*/*/.git
rm -rf pytom/*/*/*/*/*/*/*/*/*/.git
rm -rf pytom/*/*/*/*/*/*/*/*/*/*/.git
cd pytom/pytomc
make clean
cd -

zip -r pytom_$1.zip ./pytom

rm -rf ./pytom

