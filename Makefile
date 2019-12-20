CXX=g++
CXXFLAGS=-g -Wall -fPIC `root-config --cflags` -I${gpfs}/ZTrackAnalyzer/include
LDFLAGS=`root-config --glibs` -L${gpfs}/ZTrackAnalyzer/lib

libs = ZTrackUtilities TreeVariables Trigger Cuts
algs = TreeMaker MinbiasTreeMaker TruthTreeMaker TrackingEfficiency TrackingPurity TagAndProbe FCalCalibration

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
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(libs:%=-l%) $(algs:%=-l%) -o bin/run src/run.cxx

# Algorithms need to be compiled into object files and then migrated to a shared library
% : obj/%.o
	$(CXX) -shared $(LDFLAGS) $(libs:%=-l%) -o lib/lib$@.so $<

obj/%.o : src/%.cxx
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean :
	rm -rf ./lib/*.so*
	rm -rf ./bin/*
	rm -rf ./obj/*
