#CXXFLAGS=-O3 -g -Wall `root-config --cflags` -I${ATLAS_PATH}/include -I${HOME}/RootUtils -I${ATLAS_PATH}/ZTrackAnalysis/Macros
CXXFLAGS=-Ofast -g -Wall `root-config --cflags` -I${ROOT_UTILS_PATH}/include -I${ZTRACKS_PATH}/Macros
#LDFLAGS=`root-config --glibs` -L${ATLAS_PATH}/lib -L${HOME}/RootUtils -lGlobalParams -lUtilities -lAtlasUtils -lAtlasStyle
LDFLAGS=`root-config --glibs` -L${ROOT_UTILS_PATH}/lib -lUtilities

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

all : ProcessFile DoMixing

ProcessFile : ProcessFile.cxx *.h
	$(CC) -o ${ZTRACKS_PATH}/Macros/ProcessFile ${ZTRACKS_PATH}/Macros/ProcessFile.cxx

DoMixing : DoMixing.cxx
	$(CC) -o ${ZTRACKS_PATH}/Macros/DoMixing ${ZTRACKS_PATH}/Macros/DoMixing.cxx

clean :
	rm ${ZTRACKS_PATH}/Macros/ProcessFile ${ZTRACKS_PATH}/Macros/DoMixing
