
echo make && \
echo && \

make && \
echo && \

echo electrons: && \
echo && \

echo 10 GeV && \
echo && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_10.mac -ro simulations/e/10/test.root && \
#exe/Digitization simulations/e/10/test.root simulations/e/10/digit.root && \
#exe/DataAnalysis simulations/e/10/digit.root simulations/e/10/histos.root && \

echo && \
echo 100 GeV && \
echo && \

GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_100.mac -ro simulations/e/100/test.root && \
exe/Digitization simulations/e/100/test.root simulations/e/100/digit.root && \
exe/DataAnalysis simulations/e/100/digit.root simulations/e/100/histos.root && \

echo && \
echo 1 TeV && \
echo && \

GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_1000.mac -ro simulations/e/1000/test.root && \
exe/Digitization simulations/e/1000/test.root simulations/e/1000/digit.root && \
exe/DataAnalysis simulations/e/1000/digit.root simulations/e/1000/histos.root && \

echo && \
echo 10 TeV && \
echo && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_10000.mac -ro simulations/e/10000/test.root && \
#exe/Digitization simulations/e/10000/test.root simulations/e/10000/digit.root && \
#exe/DataAnalysis simulations/e/10000/digit.root simulations/e/10000/histos.root && \


echo && \
echo protons: && \
echo && \

echo 10 GeV && \
echo && \

GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_10.mac -ro simulations/p/10/test.root && \
exe/Digitization simulations/p/10/test.root simulations/p/10/digit.root && \
exe/DataAnalysis simulations/p/10/digit.root simulations/p/10/histos.root && \

echo && \
echo 100 GeV && \
echo && \

GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_100.mac -ro simulations/p/100/test.root && \
exe/Digitization simulations/p/100/test.root simulations/p/100/digit.root && \
exe/DataAnalysis simulations/p/100/digit.root simulations/p/100/histos.root && \

echo && \
echo 1 TeV && \
echo && \

GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_1000.mac -ro simulations/p/1000/test.root && \
exe/Digitization simulations/p/1000/test.root simulations/p/1000/digit.root && \
exe/DataAnalysis simulations/p/1000/digit.root simulations/p/1000/histos.root && \

echo && \
echo 10 TeV && \
echo && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_10000.mac -ro simulations/p/10000/test.root && \
#exe/Digitization simulations/p/10000/test.root simulations/p/10000/digit.root && \
#exe/DataAnalysis simulations/p/10000/digit.root simulations/p/10000/histos.root && \


echo

