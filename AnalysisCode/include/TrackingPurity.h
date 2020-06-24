#ifndef __TrackingPurity_h__
#define __TrackingPurity_h__

namespace ZTrackAnalyzer {

bool TrackingPurity (const char* directory,
                     const int dataSet,
                     const char* inFileName = "",
                     const char* eventWeightsFileName = "");

} // end namespace

#endif
