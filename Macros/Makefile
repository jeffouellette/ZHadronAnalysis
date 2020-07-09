CXXFLAGS=-O3 -g -Wall `root-config --cflags` -I${ROOT_UTILS_PATH}/include -I${ZHADRONS_PATH}/Macros
LDFLAGS=`root-config --glibs` -L${ROOT_UTILS_PATH}/lib -lUtilities

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

all : ProcessFile DoMixing

ProcessFile : ProcessFile.cxx *.h
	$(CC) -o ${ZHADRONS_PATH}/Macros/ProcessFile.exe ${ZHADRONS_PATH}/Macros/ProcessFile.cxx

DoMixing : DoMixing.cxx
	$(CC) -o ${ZHADRONS_PATH}/Macros/DoMixing.exe ${ZHADRONS_PATH}/Macros/DoMixing.cxx

clean :
	rm ${ZHADRONS_PATH}/Macros/ProcessFile.exe ${ZHADRONS_PATH}/Macros/DoMixing.exe
