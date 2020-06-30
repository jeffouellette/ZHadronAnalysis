CXX=g++
CXXFLAGS=-O3 -g -Wall -fPIC `root-config --cflags` -Iinclude -I${ROOT_UTILS_PATH}/include -I${ATLAS_PATH}/include
LDFLAGS=`root-config --glibs --ldflags` -Llib -L${ROOT_UTILS_PATH}/lib -L${ATLAS_PATH}/lib -lUtilities -lAtlasUtils -lAtlasStyle

libs = LocalUtilities TreeVariables Trigger Cuts
algs = TreeMaker MinBiasTreeMaker TruthSelectedTreeMaker TruthTreeMaker TrackingMomentumResolution TrackingEfficiency TrackingPurity TagAndProbe FCalCalibration BkgEstimator EventSkimTreeMaker

bins = run

libs : $(libs)
algs : $(algs)
bins : $(bins)
all : directories libs algs bins

directories : 
	mkdir -p lib bin obj

LocalUtilities : src/LocalUtilities.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o lib/lib$@.so src/$@.cxx

TreeVariables : src/TreeVariables.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o lib/lib$@.so src/$@.cxx

Trigger : src/Trigger.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o lib/lib$@.so src/$@.cxx

Cuts : src/Cuts.cxx TreeVariables
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lTreeVariables -shared -o lib/lib$@.so src/$@.cxx

#TrackingMomentumResolution : src/TrackingMomentumResolution.cxx
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o lib/lib$@.so src/$@.cxx
#
#TrackingEfficiency : src/TrackingEfficiency.cxx
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o lib/lib$@.so src/$@.cxx
#
#TrackingPurity : src/TrackingPurity.cxx
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o lib/lib$@.so src/$@.cxx
#
#TagAndProbe : src/TagAndProbe.cxx
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o lib/lib$@.so src/$@.cxx

# Main needs to be compiled into binary
run : $(libs:%=lib/lib%.so) src/run.cxx
	$(CXX) $(CXXFLAGS) src/run.cxx $(LDFLAGS) $(libs:%=-l%) $(algs:%=-l%) -o bin/run.exe

# Algorithms need to be compiled into object files and then migrated to a shared library
% : obj/%.o
	$(CXX) -shared $(LDFLAGS) $(libs:%=-l%) -o lib/lib$@.so $<

obj/%.o : src/%.cxx
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean :
	rm -rf ./lib/*.so*
	rm -rf ./bin/*
	rm -rf ./obj/*
