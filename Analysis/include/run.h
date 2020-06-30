#ifndef __run_h__
#define __run_h__

#include "Params.h"

using namespace ZHadronAnalysis;

namespace ZHadronAnalysis {

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

float crossSectionPicoBarns = 0;
float mcFilterEfficiency = 0;
int mcNumberEvents = 0;

CollisionSystem collisionSystem = pp15; // default is pp15

}

#endif
