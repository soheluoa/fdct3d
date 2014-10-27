#!/bin/sh
# Script to call Seismic Unix Plotting
# fill=2 for grey fill of troughs

#For complete data
pscube n1=500 n2=36 n3=100 f1num=0.00 f2num=0.0 f3num=0.0 f1=0.00 f2=0.0 f3=0.01 d1=0.0021 d2=1.0 d3=1.0 d1num=0.2 d2num=10 d3num=20 label1='t (s)' label2='Crossline (m)' label3='Inline (m)' <reconData1.bin> reconData1.ps

pscube n1=461 n2=36 n3=100 f1num=0.09 f2num=0.0 f3num=0.0 f1=0.085 f2=0.0 f3=0.01 d1=0.0021 d2=1.0 d3=1.0 d1num=0.2 d2num=10 d3num=20 label1='t (s)' label2='Crossline (m)' label3='Inline (m)' <reconData_slice1.bin> reconData_slice1.ps

pscube n1=421 n2=36 n3=100 f1num=0.17 f2num=0.0 f3num=0.0 f1=0.168 f2=0.0 f3=0.01 d1=0.0021 d2=1.0 d3=1.0 d1num=0.2 d2num=10 d3num=20 label1='t (s)' label2='Crossline (m)' label3='Inline (m)' <reconData_slice2.bin> reconData_slice2.ps
