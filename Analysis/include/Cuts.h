#ifndef __Cuts_h__
#define __Cuts_h__

#include "TreeVariables.h"

namespace ZHadronAnalysis {

bool IsElectronTrack (TreeVariables* t, const int iTrk, const int iE1 = -1, const int iE2 = -1);
bool IsMuonTrack (TreeVariables* t, const int iTrk, const int iM1 = -1, const int iM2 = -1);

} // end namespace

#endif
