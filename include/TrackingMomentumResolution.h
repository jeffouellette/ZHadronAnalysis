#ifndef __TrackingMomentumResolution_h__
#define __TrackingMomentumResolution_h__

namespace ZTrackAnalyzer {

bool TrackingMomentumResolution (const char* directory,
                                 const int dataSet,
                                 const char* inFileName = "",
                                 const char* eventWeightsFileName = "");

} // end namespace

#endif
