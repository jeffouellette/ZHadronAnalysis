#ifndef __Params_h__
#define __Params_h__

#include "ZTrackUtilities.h"

#include <TString.h>

#include <string>
#include <set>
#include <math.h>

using namespace std;

namespace ZHadronAnalysis { 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global variable declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const double pi = atan (1)*4;

const float electron_pt_cut = 20;
const float muon_pt_cut = 20;
const float trk_pt_cut = 0.8;
const float jet_pt_cut = 1;

const double electron_mass = 0.000511;
const double muon_mass = 0.105658;

TString workPath = "/atlasgpfs01/usatlas/workarea/jeff/ZHadronAnalysis/Analysis/";
TString extWorkPath = "/atlasgpfs01/usatlas/data/jeff/ZHadronAnalysis/";
TString rootPath = extWorkPath + "rootFiles/";
TString dataPath = extWorkPath + "data/";
TString plotPath = workPath + "Plots/"; 

bool useScratchDisk = false;

extern bool isPbPb;
extern bool is2015data;
extern bool isMC;
extern float crossSectionPicoBarns;
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
const float fileCentBins[30] = {
    66.402, // 80%
   148.625, // 70%, 1
   296.17,  // 60%, 2
   533.608, // 50%, 3
   659.269, // 46%, 4
   804.607, // 42%, 5
   971.487, // 38%, 6
  1162.33,  // 34%, 7
  1378.92,  // 30%, 8
  1497.54,  // 28%, 9
  1624.34,  // 26%, 10
  1759.06,  // 24%, 11
  1902.73,  // 22%, 12
  2055.77,  // 20%, 13
  2218.88,  // 18%, 14
  2393.11,  // 16%, 15
  2579.56,  // 14%, 16
  2779.68,  // 12%, 17
  2995.94,  // 10%, 18
  3110.27,  //  9%, 19
  3229.67,  //  8%, 20
  3354.66,  //  7%, 21
  3485.57,  //  6%, 22
  3622.6,   //  5%, 23
  3767.,    //  4%, 24
  3920.41,  //  3%, 25
  4083.38,  //  2%, 26
  4263.72,  //  1%, 27
  5000,      //  0%, 28, entry in array is numFileCentBins-1
  10000     // overflow for Hijing
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


const int numFileIPBins = 79;
double* GenerateIPMixBins () {
  int i = 0;
  double* ipBins = new double[numFileIPBins+1];
  for (i = 0; i < 5; i++)
    ipBins[i] = i*(1.564)/5.;
  for (i = 5; i < 23; i++)
    ipBins[i] = 1.564+(i-5)*(4.953-1.564)/18.;
  for (i = 23; i < 63; i++)
    ipBins[i] = 4.953+(i-23)*(11.083-4.953)/40.;
  for (i = 63; i < 78; i++)
    ipBins[i] = 11.083+(i-63)*(14.032-11.083)/15.;
  ipBins[78] = 14.032;
  return ipBins;
}
const double* fileIPBins = GenerateIPMixBins ();
/**
 * Returns the bin corresponding to this sum fcal et bin.
 */
short GetFileIPBin (const float ip) {
  short i = 0;
  while (i < numFileIPBins) {
    if (ip < fileIPBins[i])
      break;
    i++;
  }
  return i;
}


const double fileNtrkBins[81] = {
   0,
   4.46017,
   8.92033,
   11.9765,
   14.5842,
   17.1919,
   19.7996,
   22.312,
   24.8163,
   27.3207,
   29.8251,
   32.5971,
   35.3891,
   38.1812,
   41.0823,
   44.1872,
   47.2921,
   50.4606,
   54.0624,
   57.6641,
   61.4922,
   65.738,
   69.9839,
   74.5723,
   79.162,
   84.2263,
   89.3966,
   94.7963,
   100.256,
   106.384,
   112.668,
   119.176,
   125.905,
   132.904,
   140.27,
   147.77,
   155.623,
   163.77,
   172.14,
   180.71,
   189.899,
   199.249,
   209.028,
   219.016,
   229.4,
   240.137,
   251.174,
   263.165,
   274.432,
   286.208,
   298.925,
   311.532,
   324.533,
   338.187,
   352.468,
   366.967,
   381.965,
   397.7,
   414.07,
   430.512,
   447.859,
   466.129,
   484.939,
   504.442,
   525.371,
   546.002,
   568.148,
   591.277,
   614.442,
   639.841,
   665.697,
   693.615,
   721.418,
   750.462,
   781.428,
   814.919,
   851.232,
   892.274,
   935.023,
   988.651,
   1400
};
const int numFileNtrkBins = sizeof (fileNtrkBins) / sizeof (fileNtrkBins[0]);
/**
 * Returns the bin corresponding to this sum fcal et bin.
 */
short GetFileNtrkBin (const float ntrk) {
  short i = 0;
  while (i < numFileNtrkBins) {
    if (ntrk < fileNtrkBins[i])
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

const vector <pair <string, const set<int>*>> runGroups = {
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
const int numRunGroups = runGroups.size ();


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
