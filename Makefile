CXX=g++
CXXFLAGS=-O3 -g -Wall -fPIC `root-config --cflags` -I${gpfs}/ZTrackAnalyzer/include
LDFLAGS=`root-config --glibs --ldflags` -L${gpfs}/ZTrackAnalyzer/lib

libs = ZTrackUtilities TreeVariables Trigger Cuts
algs = TreeMaker MinbiasTreeMaker TruthTreeMaker TrackingEfficiency TrackingPurity TagAndProbe FCalCalibration BkgEstimator

bins = run

libs : $(libs)
algs : $(algs)
bins : $(bins)
all : libs algs bins

ZTrackUtilities : src/ZTrackUtilities.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o ${gpfs}/ZTrackAnalyzer/lib/lib$@.so ${gpfs}/ZTrackAnalyzer/src/$@.cxx

TreeVariables : src/TreeVariables.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o ${gpfs}/ZTrackAnalyzer/lib/lib$@.so ${gpfs}/ZTrackAnalyzer/src/$@.cxx

Trigger : src/Trigger.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o ${gpfs}/ZTrackAnalyzer/lib/lib$@.so ${gpfs}/ZTrackAnalyzer/src/$@.cxx

Cuts : src/Cuts.cxx TreeVariables
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lTreeVariables -shared -o ${gpfs}/ZTrackAnalyzer/lib/lib$@.so ${gpfs}/ZTrackAnalyzer/src/$@.cxx

TrackingEfficiency : src/TrackingEfficiency.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o ${gpfs}/ZTrackAnalyzer/lib/lib$@.so ${gpfs}/ZTrackAnalyzer/src/$@.cxx

TrackingPurity : src/TrackingPurity.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o ${gpfs}/ZTrackAnalyzer/lib/lib$@.so ${gpfs}/ZTrackAnalyzer/src/$@.cxx

TagAndProbe : src/TagAndProbe.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o ${gpfs}/ZTrackAnalyzer/lib/lib$@.so ${gpfs}/ZTrackAnalyzer/src/$@.cxx

# Main needs to be compiled into binary
run : $(libs:%=lib/lib%.so) src/run.cxx
	$(CXX) $(CXXFLAGS) src/run.cxx $(LDFLAGS) $(libs:%=-l%) $(algs:%=-l%) -o bin/run

# Algorithms need to be compiled into object files and then migrated to a shared library
% : obj/%.o
	$(CXX) -shared $(LDFLAGS) $(libs:%=-l%) -o lib/lib$@.so $<

obj/%.o : src/%.cxx
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean :
	rm -rf ./lib/*.so*
	rm -rf ./bin/*
	rm -rf ./obj/*
