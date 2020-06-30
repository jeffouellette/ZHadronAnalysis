#ifndef __TagAndProbe_h__
#define __TagAndProbe_h__

#include "TreeVariables.h"

namespace ZHadronAnalysis {

bool TagAndProbe (const char* directory,
                  const int dataSet,
                  const char* inFileName = "");

bool GetMuonTrackCut (const TreeVariables* t,
                      const int muon,
                      const string cutLevel = "HILoose");

} // end namespace

#endif
