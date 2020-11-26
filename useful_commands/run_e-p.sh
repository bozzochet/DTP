
make && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_10.mac -ro simulations/e/10/test.root && \
#exe/Digitization simulations/e/10/test.root simulations/e/10/digit.root && \
exe/DataAnalysis simulations/e/10/digit.root simulations/e/10/histos.root 10 10 0.5 && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/p_10.mac -ro simulations/p/10/test.root && \
#exe/Digitization simulations/p/10/test.root simulations/p/10/digit.root && \
exe/DataAnalysis simulations/p/10/digit.root simulations/p/10/histos.root 10 1 0.5 && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_100.mac -ro ext/e/100/test.root && \
#exe/Digitization ext/e/100/test.root ext/e/100/digit.root && \
exe/DataAnalysis simulations/e/100/digit.root simulations/e/100/histos.root 100 100 0.5 && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/p_100.mac -ro ext/p/100/test.root && \
#exe/Digitization ext/p/100/test.root ext/p/100/digit.root && \
exe/DataAnalysis simulations/p/100/digit.root simulations/p/100/histos.root 100 10 0.5 && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_1000.mac -ro ext/e/1000/test.root && \
#exe/Digitization ext/e/1000/test.root ext/e/1000/digit.root && \
#exe/DataAnalysis ext/e/1000/digit.root ext/e/1000/histos.root && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/p_1000.mac -ro ext/p/1000/test.root && \
#exe/Digitization ext/p/1000/test.root ext/p/1000/digit.root && \
#exe/DataAnalysis ext/p/1000/digit.root ext/p/1000/histos.root && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/e_10000.mac -ro ext/e/10000/test.root && \
#exe/Digitization ext/e/10000/test.root ext/e/10000/digit.root && \
#exe/DataAnalysis ext/e/10000/digit.root ext/e/10000/histos.root && \

#GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/p_10000.mac -ro ext/p/10000/test.root && \
#exe/Digitization ext/p/10000/test.root ext/p/10000/digit.root && \
#exe/DataAnalysis ext/p/10000/digit.root ext/p/10000/histos.root && \

echo
