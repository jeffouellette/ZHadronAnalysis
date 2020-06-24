#ifndef __FCalCalibration_h__
#define __FCalCalibration_h__

namespace ZTrackAnalyzer {

// Set of percentile-space centrality cuts, intended for plotting # events vs. centrality %
const float percentileCentBins[81] = {
    5000    /*  0% */,  4263.72 /*  1% */,  4083.38 /*  2% */,  3920.41 /*  3% */,    3767.   /*  4% */,
    3622.6  /*  5% */,  3485.57 /*  6% */,  3354.66 /*  7% */,  3229.67 /*  8% */,    3110.27 /*  9% */,
    2995.94 /* 10% */,  2885.71 /* 11% */,  2779.68 /* 12% */,  2677.6  /* 13% */,    2579.56 /* 14% */,
    2484.75 /* 15% */,  2393.11 /* 16% */,  2304.47 /* 17% */,  2218.88 /* 18% */,    2136.17 /* 19% */,
    2055.77 /* 20% */,  1978.02 /* 21% */,  1902.73 /* 22% */,  1829.74 /* 23% */,    1759.06 /* 24% */,
    1690.47 /* 25% */,  1624.34 /* 26% */,  1560.05 /* 27% */,  1497.54 /* 28% */,    1437.29 /* 29% */,
    1378.92 /* 30% */,  1322.12 /* 31% */,  1267.07 /* 32% */,  1213.91 /* 33% */,    1162.33 /* 34% */,
    1112.2  /* 35% */,  1063.71 /* 36% */,  1016.86 /* 37% */,  971.487 /* 38% */,    927.582 /* 39% */,
    885.172 /* 40% */,  844.192 /* 41% */,  804.607 /* 42% */,  766.305 /* 43% */,    729.251 /* 44% */,
    693.606 /* 45% */,  659.269 /* 46% */,  626.047 /* 47% */,  594.07  /* 48% */,    563.263 /* 49% */,
    533.608 /* 50% */,  505.085 /* 51% */,  477.734 /* 52% */,  451.509 /* 53% */,    426.354 /* 54% */,
    402.144 /* 55% */,  378.968 /* 56% */,  356.885 /* 57% */,  335.738 /* 58% */,    315.523 /* 59% */,
    296.17  /* 60% */,  277.742 /* 61% */,  260.219 /* 62% */,  243.588 /* 63% */,    227.751 /* 64% */,
    212.684 /* 65% */,  198.428 /* 66% */,  184.922 /* 67% */,  172.097 /* 68% */,    160.033 /* 69% */,
    148.625 /* 70% */,  137.874 /* 71% */,  127.745 /* 72% */,  118.249 /* 73% */,    109.333 /* 74% */,
    100.928 /* 75% */,  93.071  /* 76% */,  85.729  /* 77% */,  78.834  /* 78% */,    72.411  /* 79% */,
    66.402    /* 80% , entry in array is numPercentileCentBins-1 */
};
const int percentileCentCuts[81] = {
    0,                  1,                  2,                  3,                    4,
    5,                  6,                  7,                  8,                    9,
    10,                 11,                 12,                 13,                   14,
    15,                 16,                 17,                 18,                   19,
    20,                 21,                 22,                 23,                   24,
    25,                 26,                 27,                 28,                   29,
    30,                 31,                 32,                 33,                   34,
    35,                 36,                 37,                 38,                   39,
    40,                 41,                 42,                 43,                   44,
    45,                 46,                 47,                 48,                   49,
    50,                 51,                 52,                 53,                   54,
    55,                 56,                 57,                 58,                   59,
    60,                 61,                 62,                 63,                   64,
    65,                 66,                 67,                 68,                   69,
    70,                 71,                 72,                 73,                   74,
    75,                 76,                 77,                 78,                   79,
    80
};

const short numPercentileCentBins = sizeof (percentileCentBins) / sizeof (percentileCentBins[0]);

// Returns which centrality percentile this fcal_et is in, e.g. GetPercentileCentrality (67) returns 80.
// If collision is more peripheral than 80%, will return 100. If centrality is more central than upper bound, will return -1.
int GetPercentileCentrality (const float fcal_et) {
  short i = 0;
  while (i < numPercentileCentBins) {
    if (fcal_et > percentileCentBins[i])
      break;
    i++;
  }
  if (i == 0)
    return -1;
  if (i == numPercentileCentBins)
    return 100;
  return i;
}


bool FCalCalibration (const char* directory,
                      const int dataSet,
                      const char* inFileName = "");

} // end namespace

#endif
