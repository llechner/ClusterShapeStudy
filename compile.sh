g++ -std=c++17 $CMSSW_BASE/src/ClusterShape/src/*.cc -o $CMSSW_BASE/src/ClusterShape/clusterShape `root-config --cflags --glibs` -lstdc++fs -I$CMSSW_BASE/src/

