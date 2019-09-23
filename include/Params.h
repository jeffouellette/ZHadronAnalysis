#ifndef __Params_h__
#define __Params_h__

#include "ZTrackUtilities.h"

#include <TString.h>

#include <string>
#include <math.h>

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

} // end namespace

#endif
