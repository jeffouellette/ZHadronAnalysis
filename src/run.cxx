#include "TreeMaker.h"
#include "MinbiasTreeMaker.h"
#include "TruthTreeMaker.h"
#include "TrackingEfficiency.h"
#include "TrackingPurity.h"
#include "TagAndProbe.h"
#include "FCalCalibration.h"
#include "BkgEstimator.h"
#include "Params.h"

#include <string>
#include <iostream>

using namespace std;
using namespace ZTrackAnalyzer;

namespace ZTrackAnalyzer {

bool isMC = false;
bool isPbPb = false;
bool is2015data = false;

bool doElectronPtUpVar = false;
bool doElectronPtDownVar = false;
bool doMuonPtUpVar = false;
bool doMuonPtDownVar = false;
bool doElectronLHMediumVar = false;
bool doMuonTightVar = false;
bool doHITightVar = false;
bool doPionsOnlyVar = false;

float mcFilterEfficiency = 0;
int mcNumberEvents = 0;

CollisionSystem collisionSystem = pp15; // default is pp15

}

int main (int argc, char** argv) {

  int argn                = 1;
  useScratchDisk          = string (argv[argn++]) == "true";
  const string alg        = string (argv[argn++]);
  const char* subdir      = argv[argn++];
  const int dataSet       = atoi (argv[argn++]);
  isMC                    = string (argv[argn++]) == "true";

  string _collSys         = string (argv[argn++]);
  if (_collSys == "PbPb15")
    collisionSystem = PbPb15;
  else if (_collSys == "pp17")
    collisionSystem = pp17;
  else if (_collSys == "PbPb18")
    collisionSystem = PbPb18;
  else {
    cout << "In Run.cxx: Invalid collision system, exiting." << endl;
    return 1;
  }

  isPbPb = (collisionSystem == PbPb15 || collisionSystem == PbPb18);
  is2015data = (collisionSystem == PbPb15);

  doElectronPtUpVar       = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);
  doElectronPtDownVar     = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);
  doMuonPtUpVar           = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);
  doMuonPtDownVar         = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);
  doElectronLHMediumVar   = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);
  doMuonTightVar          = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);
  doHITightVar            = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);
  doPionsOnlyVar          = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);


  mcFilterEfficiency                = (argc > argn && argv[argn] ? atof (argv[argn++]) : 0);
  mcNumberEvents                    = (argc > argn && argv[argn] ? atoi (argv[argn++]) : 0);
  const char* inFileName            = (argc > argn && argv[argn] ? argv[argn++] : "");
  const char* eventWeightsFileName  = (argc > argn && argv[argn] ? argv[argn++] : "");

  cout << "Info: In run.cxx: Configuration set to";
  cout << "\n\tuseScratchDisk = " << useScratchDisk;
  cout << "\n\talg = " << alg;
  cout << "\n\tsubdir = " << subdir;
  cout << "\n\tdataSet = " << dataSet;
  cout << "\n\tisMC = " << isMC;
  cout << "\n\tCollisionSystem = " << collisionSystem;
  cout << "\n\tisPbPb = " << isPbPb;
  cout << "\n\tis2015data = " << is2015data;
  cout << "\n\tinFileName = " << inFileName;
  cout << "\n\tdoElectronPtUpVar = " << doElectronPtUpVar;
  cout << "\n\tdoElectronPtDownVar = " << doElectronPtDownVar;
  cout << "\n\tdoMuonPtUpVar = " << doMuonPtUpVar;
  cout << "\n\tdoMuonPtDownVar = " << doMuonPtDownVar;
  cout << "\n\tdoElectronLHMediumVar = " << doElectronLHMediumVar;
  cout << "\n\tdoMuonTightVar = " << doMuonTightVar;
  cout << "\n\tdoHITightVar = " << doHITightVar;
  cout << "\n\tdoPionsOnlyVar = " << doPionsOnlyVar;
  cout << "\n\tmcFilterEfficiency = " << mcFilterEfficiency;
  cout << "\n\tmcNumberEvents = " << mcNumberEvents;
  cout << "\n\tinFileName = " << inFileName;
  cout << "\n\teventWeightsFileName = " << eventWeightsFileName;
  cout << endl;


  bool success = false;
  if (alg == "TreeMaker") {
    cout << "Info: In run.cxx: Running TreeMaker algorithm..." << endl;
    success = TreeMaker (subdir, dataSet, inFileName);
  }
  else if (alg == "MinbiasTreeMaker") {
    cout << "Info: In run.cxx: Running MinbiasTreeMaker algorithm..." << endl;
    success = MinbiasTreeMaker (subdir, dataSet, inFileName);
  }
  else if (alg == "TruthTreeMaker") {
    cout << "Info: In run.cxx: Running TruthTreeMaker algorithm..." << endl;
    success = TruthTreeMaker (subdir, dataSet, inFileName);
  }
  else if (alg == "TrackingEfficiency") {
    cout << "Info: In run.cxx: Running TrackingEfficiency algorithm..." << endl;
    success = TrackingEfficiency (subdir, dataSet, inFileName, eventWeightsFileName);
  }
  else if (alg == "TrackingPurity") {
    cout << "Info: In run.cxx: Running TrackingPurity algorithm..." << endl;
    success = TrackingPurity (subdir, dataSet, inFileName, eventWeightsFileName);
  }
  else if (alg == "TagAndProbe") {
    cout << "Info: In run.cxx: Running TagAndProbe algorithm..." << endl;
    success = TagAndProbe (subdir, dataSet, inFileName);
  }
  else if (alg == "FCalCalibration") {
    cout << "Info: In run.cxx: Running FCalCalibration algorithm..." << endl;
    success = FCalCalibration (subdir, dataSet, inFileName);
  }
  else if (alg == "BkgEstimator") {
    cout << "Info: In run.cxx: Running BkgEstimator algorithm..." << endl;
    success = BkgEstimator (subdir, dataSet, inFileName);
  }
  else {
    cout << "Error: In run.cxx: Failed to recognize algorithm! Quitting." << endl;
    return 1;
  }


  if (success) {
    cout << "Info: In run.cxx: Finished making trees!" << endl;
    return 0;
  }
  else {
    cout << "Error: In run.cxx: Algorithm failed!" << endl;
    return 1;
  }
}
