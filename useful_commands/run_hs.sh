#!/bin/bash

#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

GGSPenny -nt 8 -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/run_pro_hs.mac -ro GGSOutput_pro_hs.root
GGSPenny -nt 8 -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/run_ele_hs.mac -ro GGSOutput_ele_hs.root

#./exe/Digitization GGSOutput_pro_hs.root DigitOut_pro_hs.root
#./exe/Digitization GGSOutput_ele_hs.root DigitOut_ele_hs.root

#./exe/DataAnalysis_calo DigitOut_pro_hs.root AnaOut_pro_hs.root 1000 600 800
#./exe/DataAnalysis_calo DigitOut_ele_hs.root AnaOut_ele_hs.root 1000 600 800
