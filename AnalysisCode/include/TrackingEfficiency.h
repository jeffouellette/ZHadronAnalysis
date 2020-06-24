#ifndef __TrackingEfficiency_h__
#define __TrackingEfficiency_h__

namespace ZTrackAnalyzer {

bool TrackingEfficiency (const char* directory,
                         const int dataSet,
                         const char* inFileName = "",
                         const char* eventWeightsFileName = "");

} // end namespace

#endif
