#!/usr/bin/env bash


export PYTHONPATH='../..':'../pytomc/swigModules/':$PYTHONPATH

epydoc  --graph all --pstat=./pstat -v --config=./epydoc.conf --inheritance listed
