#ifndef __Params_h__
#define __Params_h__

#include "ZTrackUtilities.h"

#include <TString.h>

#include <string>
#include <set>
#include <math.h>

using namespace std;

namespace ZTrackAnalyzer { 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global variable declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const double pi = atan (1)*4;

const float electron_pt_cut = 20;
const float muon_pt_cut = 20;
const float trk_pt_cut = 1;
const float jet_pt_cut = 1;

const double electron_mass = 0.000511;
const double muon_mass = 0.105658;

TString workPath = "/atlasgpfs01/usatlas/workarea/jeff/ZTrackAnalyzer/";
TString extWorkPath = "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/";
TString rootPath = extWorkPath + "rootFiles/";
TString dataPath = extWorkPath + "data/";
TString plotPath = workPath + "Plots/"; 

bool useScratchDisk = false;

extern bool isPbPb;
extern bool is2015data;
extern bool isMC;
extern float mcFilterEfficiency;
extern int mcNumberEvents;

// systematics configuration variables
extern bool doElectronPtUpVar;
extern bool doElectronPtDownVar;
extern bool doMuonPtUpVar;
extern bool doMuonPtDownVar;
extern bool doElectronLHMediumVar;
extern bool doMuonTightVar;
extern bool doHITightVar;
extern bool doPionsOnlyVar;

bool doNchWeighting = false;

const int nElectronTrig = 3;
const TString electronTrigNames[nElectronTrig] = {"HLT_e15_lhloose_L1EM12", "HLT_e15_lhloose_ion_L1EM12", "HLT_e15_loose_ion_L1EM12"};
const float electronTrigMinPtCuts[nElectronTrig] = {15, 15, 15};
const float electronTrigMaxPtCuts[nElectronTrig] = {10000, 10000, 10000};
const float electronTrigMinEtaCuts[nElectronTrig] = {-2.47, -2.47, -2.47};
const float electronTrigMaxEtaCuts[nElectronTrig] = {2.47, 2.47, 2.47};

const int nMuonTrig = 1;
const TString muonTrigNames[nMuonTrig] = {"HLT_mu14"};
const float muonTrigMinPtCuts[nMuonTrig] = {14};
const float muonTrigMaxPtCuts[nMuonTrig] = {10000};
const float muonTrigMinEtaCuts[nMuonTrig] = {-2.5};
const float muonTrigMaxEtaCuts[nMuonTrig] = {2.5};

enum CollisionSystem { pp15 = 0, PbPb15 = 1, pPb16 = 2, Pbp16 = 3, XeXe17 = 4, pp17 = 5, PbPb18 = 6 }; // run 2 HI data sets
extern CollisionSystem collisionSystem;


//enum NucleonPair { pp = 0, pn = 1, np = 2, nn = 3 };
//extern NucleonPair inNuclPair;
//
//const char* NucleonPairToString (const NucleonPair p) {
//  switch (p) {
//    case pp:
//      return "pp";
//    case pn:
//      return "pn";
//    case np:
//      return "np";
//    case nn:
//      return "nn";
//  }
//  return "";
//}
//
//NucleonPair StringToNucleonPair (const std::string str) {
//  if (str == std::string ("pp"))
//    return pp;
//  if (str == std::string ("pn"))
//    return pn;
//  if (str == std::string ("np"))
//    return np;
//  if (str == std::string ("nn"))
//    return nn;
//  return pp;
//}

// Centrality classes for mixing events
const float fileCentBins[13] = {
    66.402, // 80%
   296.17,  // 60%
   885.172, // 40%
  1378.92,  // 30%
  1690.47,  // 25%
  2055.77,  // 20%
  2484.75,  // 15%
  2995.94,  // 10%
  3229.67,  //  8%
  3485.57,  //  6%
  3767.,    //  4%
  4083.38,  //  2%
  5000      //  0%,   entry in array is numFileCentBins-1
};
const int numFileCentBins = sizeof (fileCentBins) / sizeof (fileCentBins[0]);

/**
 * Returns the bin corresponding to this sum fcal et bin.
 */
short GetFileCentBin (const float fcal_et) {
  short i = 0;
  while (i < numFileCentBins) {
    if (fcal_et < fileCentBins[i])
      break;
    i++;
  }
  return i;
}


const set<int> group1 = {286711, 286717, 286748, 286767, 286834, 286854, 286908, 286990};
const set<int> group2 = {287038, 287044, 287068, 287222, 287224, 287259, 287270, 287281};
const set<int> group3 = {287321, 287330, 287334, 287378, 287380, 287382, 287560, 287594};
const set<int> group4 = {287632, 287706, 287728, 287827};
const set<int> group5 = {287843, 287866, 287924, 287931};
const set<int> groupA = {365502, 365512, 365573, 365602, 365627, 365678, 365681, 365709};
const set<int> groupB = {365752, 365834};
const set<int> groupC = {365914, 365932};
const set<int> groupD = {366011, 366029, 366092};
const set<int> groupE = {366142, 366268, 366337};
const set<int> groupF = {366383, 366413, 366476};
const set<int> groupG = {366526, 366528, 366627, 366691, 366754, 366805};
const set<int> groupH = {366860, 366878, 366919, 366931, 366994, 367023};
const set<int> groupI = {367099, 367134, 367165};
const set<int> groupJ = {367170, 367233, 367273};
const set<int> groupK = {367318, 367321, 367363, 367364, 367365, 367384};
//const vector<const set<int>*> runGroups = {&group1, &group2, &group3, &group4, &group5, &groupA, &groupB, &groupC, &groupD, &groupE, &groupF, &groupG, &groupH, &groupI, &groupJ, &groupK};

const set <pair <string, const set<int>*>> runGroups = {
  //{"Group1", &group1},
  //{"Group2", &group2},
  //{"Group3", &group3},
  //{"Group4", &group4},
  //{"Group5", &group5},
  {"GroupA", &groupA},
  {"GroupB", &groupB},
  {"GroupC", &groupC},
  {"GroupD", &groupD},
  {"GroupE", &groupE},
  {"GroupF", &groupF},
  {"GroupG", &groupG},
  {"GroupH", &groupH},
  {"GroupI", &groupI},
  {"GroupJ", &groupJ},
  {"GroupK", &groupK}
};


TString GetRunGroupTString (int rn) {
  for (const auto& group : runGroups) {
    if (group.second->find (rn) != group.second->end ())
      return TString (group.first);
  }
  return "";
}


short GetRunGroup (int rn) {
  short rg = 0;
  for (const auto& group : runGroups) {
    if (group.second->find (rn) == group.second->end ())
      rg++;
    else
      return rg;
  }
  return -1;
}

} // end namespace

#endif
