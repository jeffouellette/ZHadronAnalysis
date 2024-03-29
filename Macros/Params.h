#ifndef __Params_h__
#define __Params_h__

#include <Utilities.h>
#include <ArrayTemplates.h>

#include <set>
#include <map>
#include <iostream>
#include <algorithm>

#include <math.h>

#include <TColor.h>
#include <TPave.h>
#include <TLine.h>
#include <TMarker.h>

using namespace std;

typedef TGraphAsymmErrors TGAE;

TColor* tcolor = new TColor ();

//const Style_t markerStyles[7] = {kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCrossX, kFullCircle, kFullSquare, kFullDiamond};
const Style_t markerStyles[7] = {kFullCircle, kFullSquare, kFullDiamond, kFullCrossX, kFullCircle, kFullSquare, kFullDiamond};
//const Color_t colors[11] =    {kBlack, kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};
const Color_t colors[11] =    {kBlack, kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};
const Color_t fillColors[11] = {kGray, kRed-9, kAzure-9, kGreen-10, kViolet-2, kMagenta-9, kCyan-6, kOrange, kGreen-2, kGray, kPink-4};
const Color_t modelFillColors[8] = {kGray, kRed+2, kBlue-3, kGreen+3, kViolet-6, kMagenta+2, kCyan+3, kOrange+2};

const Color_t finalColors[4]          = {(Color_t) tcolor->GetColor (0, 0, 0), (Color_t) tcolor->GetColor ( 87, 132, 198), (Color_t) tcolor->GetColor (130,  10, 130), (Color_t) tcolor->GetColor (255,  12,  73)};
const Color_t finalFillColors[4]      = {(Color_t) tcolor->GetColor (0, 0, 0), (Color_t) tcolor->GetColor ( 87, 132, 198), (Color_t) tcolor->GetColor (130,  10, 130), (Color_t) tcolor->GetColor (255,  12,  73)};
const Color_t finalModelFillColors[4] = {(Color_t) tcolor->GetColor (0, 0, 0), (Color_t) tcolor->GetColor ( 87, 132, 198), (Color_t) tcolor->GetColor (130,  10, 130), (Color_t) tcolor->GetColor (255,  12,  73)};
//const Color_t finalModelFillColors[4] = {kBlack, kAzure-9, kViolet-9, kRed-9};
const float fillAlpha = 1;

const double pi = M_PI;

double* doubleLogSpace (const double lo, const double hi, const int num) {
  double* arr = Get1DArray <double> (2*num+1);
  double* logarr = logspace (lo, hi, num-1);
  for (int i = 0; i < num; i++) {
    arr[num+i+1] = logarr[i];
    arr[i] = -logarr[num-1-i];
  }
  arr[num] = 0;
  delete logarr;
  return arr;
}

float max_iaa = 3.2;
//float max_icp = 1.4;
float max_rel_sys = 40;

extern int mixingFraction;

const TString rootPath = TString (std::getenv ("ZHADRONS_DATA_PATH")) + "/rootFiles/";
const TString plotPath = TString (std::getenv ("ZHADRONS_PATH")) + "/Plots/";
const TString tablesPath = TString (std::getenv ("ZHADRONS_PATH")) + "/Tables/";

const double etaTrkBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
const int numEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;
const int numFineEtaTrkBins = 40;
const double* fineEtaTrkBins = linspace (-2.5, 2.5, numFineEtaTrkBins);


// Centrality cuts in impact parameter from Glauber
const double ipCentBins[4] = {14.031, 8.582, 4.952, 0};
const int numIPCentBins = sizeof (ipCentBins) / sizeof (ipCentBins[0]);
// extended range: 0%     1%      10%     20%     30%     40%     50%     60%     70%     80%
//                 0      1.564   4.953   7.009   8.582   9.911   11.083  12.136  13.111  14.032
short GetIPCentBin (const float ip) {
  short i = 0;
  while (i < numIPCentBins) {
    if (ip >= ipCentBins[i])
      break;
    i++;
  }
  return i;
}



// Centrality cuts in GeV:  80%  -  30%  -  10%  -   0%
//const double centBins[4] = {63.719, 1368.75, 2989.31, 5000}; // 2015 recommendations, old Glauber MC
const double centBins[4] = {66.402, 1378.92, 2995.94, 5000}; // updated 2015 recommendations, new Glauber MC
const int centCuts[4]    = {80,     30,     10,      0};

const int numCentBins = sizeof (centBins) / sizeof (centBins[0]); // no minus 1 to include pp bin

short GetCentBin (const float fcal_et) {
  short i = 0;
  while (i < numCentBins) {
    if (fcal_et < centBins[i])
      break;
    i++;
  }
  return i;
}

// Centrality cuts for bin-by-bin corrections (<--> "bbbCorr")
const double bbbCorrCentBins[4] = {66.402, 1378.92, 2995.94, 5000}; // updated 2015 recommendations, new Glauber MC
const int bbbCorrCentCuts[4]    = {80,     30,     10,      0};

const int numBBBCorrCentBins = sizeof (bbbCorrCentBins) / sizeof (bbbCorrCentBins[0]); // no minus 1 to include pp bin

short GetBBBCorrCentBin (const float fcal_et) {
  short i = 0;
  while (i < numBBBCorrCentBins) {
    if (fcal_et < bbbCorrCentBins[i])
      break;
    i++;
  }
  return i;
}

// Centrality cuts for track-based corrections (purity, efficiency)
const double trkCorrCentBins[4] = {66.402, 1378.92, 2995.94, 5000}; // updated 2015 recommendations, new Glauber MC
const int trkCorrCentCuts[4]    = {80,     30,     10,      0};

const int numTrkCorrCentBins = sizeof (trkCorrCentBins) / sizeof (trkCorrCentBins[0]); // no minus 1 to include pp bin

short GetTrkCorrCentBin (const float fcal_et) {
  short i = 0;
  while (i < numTrkCorrCentBins) {
    if (fcal_et < trkCorrCentBins[i])
      break;
    i++;
  }
  return i;
}

// Finer centrality cuts for slightly increased precision measurements
const double fineCentBins[10] = {66.402, 148.625, 296.17, 533.608, 885.172, 1378.92, 2055.77, 2995.94, 3622.6, 5000};
const int fineCentCuts[10]    = {80,     70,      60,     50,      40,      30,      20,      10,      5,      0};
const int numFineCentBins = sizeof (fineCentBins) / sizeof (fineCentBins[0]);
short GetFineCentBin (const float fcal_et) {
  short i = 0;
  while (i < numFineCentBins) {
    if (fcal_et < fineCentBins[i])
      break;
    i++;
  }
  return i;
}


// delta-phi bin edges between Z and tracks
const double phiLowBins[3] = {0,      3*pi/4,     15*pi/16};
const double phiHighBins[3] = {pi/2,  15*pi/16,   pi};
const int numPhiBins = sizeof (phiLowBins) / sizeof (phiLowBins[0]);

// Returns which delta-phi bin corresponds to this relative angle
short GetidPhi (const float dphi) {
  short idPhi = 0;
  while (idPhi < numPhiBins) {
    if (phiLowBins[idPhi] < dphi && dphi < phiHighBins[idPhi])
      break;
    else
      idPhi++;
  }
  return idPhi;
}

// Returns how many bins should be used in the dPhi correlation plot for this bin
short GetNdPhiBins (const short iPtch, const short iCent) {
  short nBins = 40;
  if (iPtch > 3) nBins /= 2;
  if (iCent > 0) nBins /= 2;
  return nBins;
}


// Z boson pT bins used in analysis
const double zPtBins[6] = {0, 5, 15, 30, 60, 10000};
const int nPtZBins = sizeof (zPtBins) / sizeof (zPtBins[0]) - 1;
short GetPtZBin (const float zPt) {
  short iPtZ = 0;
  while (iPtZ < nPtZBins) {
    if (zPt < zPtBins[iPtZ+1])
      break;
    else
      iPtZ++;
  }
  return iPtZ;
}


// min/max values on track kinematics
const double trk_min_pt = 1;
const double trk_max_pt = 60;


// code that instantiates bins in track pT.
// essentially just a set-of-a-set of bin edges where the highest bin edge is min (Z pT)
const int maxNPtchBins = 8;
double allPtchBins[9] = {1, 2, 4, 8, 15, 30, 60, 120, 240};
int* nPtchBins = Get1DArray <int> (nPtZBins);

double** init_pTchBins () {
  double** _pTchBins = Get1DArray <double*> (nPtZBins);
  for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    nPtchBins[iPtZ] = 0;
    while (nPtchBins[iPtZ] < maxNPtchBins && allPtchBins[nPtchBins[iPtZ]] < zPtBins[iPtZ])
      nPtchBins[iPtZ]++;

    if (nPtchBins[iPtZ] == 0)
      nPtchBins[iPtZ] = 1;
  }

  for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    _pTchBins[iPtZ] = Get1DArray<double> (nPtchBins[iPtZ]+1);
    for (int iPtch = 0; iPtch <= nPtchBins[iPtZ]; iPtch++) {
      _pTchBins[iPtZ][iPtch] = allPtchBins[iPtch];
    }
  }

  return _pTchBins;
}
double** pTchBins = init_pTchBins (); // iPtZ, iPtch

// code that instantiates bins in track pT / Z pT.
// essentially just a set-of-a-set of bin edges where the lowest bin is min (track pT) / min (Z pT)
const int maxNXhZBins = 8;
double allXhZBins[9] = {1./120., 1./60., 1./30., 1./15., 1./8., 1./4., 1./2., 1., 2.};
int* nXhZBins = Get1DArray <int> (nPtZBins);

double** init_xhZBins () {
  double** _xhZBins = Get1DArray <double*> (nPtZBins);

  for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    nXhZBins[iPtZ] = maxNXhZBins;
    while (nXhZBins[iPtZ] > 0 && zPtBins[iPtZ] > 0 && allXhZBins[maxNXhZBins-nXhZBins[iPtZ]] < 1./zPtBins[iPtZ])
      nXhZBins[iPtZ]--;
    nXhZBins[iPtZ]--;

    //nXhZBins[iPtZ]--; // for rebinning 1st bin

    if (nXhZBins[iPtZ] == 0)
      nXhZBins[iPtZ] = 1;
  }

  for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    _xhZBins[iPtZ] = Get1DArray<double> (nXhZBins[iPtZ]+1);
    for (int iXhZ = 0; iXhZ <= nXhZBins[iPtZ]; iXhZ++) {
      //_xhZBins[iPtZ][iXhZ] = allXhZBins[iXhZ+(maxNXhZBins-nXhZBins[iPtZ]-1 - (iXhZ==0))]; // for rebinning 1st bin
      _xhZBins[iPtZ][iXhZ] = allXhZBins[iXhZ+(maxNXhZBins-nXhZBins[iPtZ]-1)];
    }
  }

  return _xhZBins;
}
double** xhZBins = init_xhZBins ();

// Returns which bin corresponds to this value of zH for this Z boson pT bin number
short GetiZH (const float zH, const int iPtZ) {
  short iZH = 0;
  while (iZH < nXhZBins[iPtZ]) {
    if (xhZBins[iPtZ][iZH+1] < zH)
      iZH++;
    else
      break;
  }
  return iZH;
}


// Prints track pT bins for this Z pT bin
void PrintPtBins (const short iPtZ) {
  for (int i = 0; i <= nPtchBins[iPtZ]; i++) {
    cout << pTchBins[iPtZ][i] << endl;
  }
}

// Prints track pT / Z pT bins for this Z pT bin
void PrintXhZBins (const short iPtZ) {
  for (int i = 0; i <= nXhZBins[iPtZ]; i++) {
    cout << xhZBins[iPtZ][i] << endl;
  }
}

// Prints track pT / Z pT bin edges in a LaTeX friendly format (intended to be placed in a table)
void LatexXhZBins () {
  cout << "\\hline" << endl;
  for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    cout << "\\multicolumn{2}{";
    if (iPtZ == 2)
      cout << "|";
    cout << " c |";
    if (iPtZ != nPtZBins-1)
      cout << "|";
    cout << "}{";
    cout << "\\pt^\\mathrm{Z} = \\SI{" << zPtBins[iPtZ] << "}{\\GeV}} & ";
  }
  cout << "\\\\ \\hline" << endl;
  for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    cout << "$x_\\mathrm{hZ}$ & $\\pt^\\mathrm{ch}$ [\\GeV] ";
    if (iPtZ != nPtZBins-1)
      cout << "& ";
  }
  cout << "\\\\ \\hline" << endl;
  for (int i = 0; i <= maxNXhZBins; i++) {
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (i < nXhZBins[iPtZ])
        cout << xhZBins[iPtZ][i] << " & " << xhZBins[iPtZ][i] * zPtBins[iPtZ] << " ";
      if (iPtZ != nPtZBins-1)
        cout << "& ";
    }
    //if (i != nXhZBins[iPtZ])
    cout << "\\\\ \\hline";
    cout << endl;
  }

}


// Calculates the greatest-common-denominator of a & b, for use in converting decimal to fractional representations
long myGCD (long a, long b) {
  if (a == 0)
    return b;
  else if (b == 0)
    return a;

  if (a < b)
    return myGCD (a, b % a);
  else
    return myGCD (b, a % b);
}


// Returns a string of phi in units of pi, e.g. GetPiString (1.57079633) would return "#pi/2"
const char* GetPiString (double phi) {
  phi = phi / pi;
  while (phi < 0)
    phi += 2;
  while (phi >= 2)
    phi -= 2;

  if (phi == 0)
    return "0";

  int denom, numer;

  double integerFloor = floor (phi);
  double fracPart = phi - integerFloor;

  if (integerFloor == phi) {
    denom = 1;
    numer = phi;
  }
  else {
    const long precision = 1000000000;
    long _gcd = myGCD (round (fracPart * precision), precision);

    denom = (int)(precision / _gcd);
    numer = (int)(round (fracPart * precision) / _gcd);
  }

  if (numer == 1 && denom == 1)
    return "#pi";
  else if (numer == 1)
    return Form ("#pi/%i", denom);
  else if (denom == 1)
    return Form ("%i#pi", numer);
  else
    return Form ("%i#pi/%i", numer, denom);
}


// Set of very fine centrality cuts, intended for use as mixed event categories
const float superFineCentBins[157] = {
    66.402, // 80%
    //0.5*66.402 + 0.5*72.411,
    72.411, // 79%
    //0.5*72.411 + 0.5*78.834,
    78.834, // 78%
    //0.5*78.834 + 0.5*85.729,
    85.729, // 77%
    //0.5*85.729 + 0.5*93.071,
    93.071, // 76%
    //0.5*93.071 + 0.5*100.928,
    100.928, // 75%
    //0.5*100.928 + 0.5*109.333,
    109.333, // 74%
    //0.5*109.333 + 0.5*118.249,
    118.249, // 73%
    //0.5*118.249 + 0.5*127.745,
    127.745, // 72%
    //0.5*127.745 + 0.5*137.874,
    137.874, // 71%
    //0.5*137.874 + 0.5*148.625,
    148.625, // 70%
    //0.5*148.625 + 0.5*160.033,
    160.033, // 69%
    //0.5*160.033 + 0.5*172.097,
    172.097, // 68%
    //0.5*172.097 + 0.5*184.922,
    184.922, // 67%
    //0.5*184.922 + 0.5*198.428,
    198.428, // 66%
    //0.5*198.428 + 0.5*212.684,
    212.684, // 65%
    //0.5*212.684 + 0.5*227.751,
    227.751, // 64%
    //0.5*227.751 + 0.5*243.588,
    243.588, // 63%
    //0.5*243.588 + 0.5*260.219,
    260.219, // 62%
    //0.5*260.219 + 0.5*277.742,
    277.742, // 61%
    //0.5*277.742 + 0.5*296.17,
    296.17,  // 60%
    //0.5*296.17 + 0.5*315.523,
    315.523, // 59%
    //0.5*315.523 + 0.5*335.738,
    335.738, // 58%
    //0.5*335.738 + 0.5*356.885,
    356.885, // 57%
    //0.5*356.885 + 0.5*378.968,
    378.968, // 56%
    //0.5*378.968 + 0.5*402.144,
    402.144, // 55%
    //0.5*402.144 + 0.5*426.354,
    426.354, // 54%
    //0.5*426.354 + 0.5*451.509,
    451.509, // 53%
    //0.5*451.509 + 0.5*477.734,
    477.734, // 52%
    //0.5*477.734 + 0.5*505.085,
    505.085, // 51%
    //0.5*505.085 + 0.5*533.608,
    533.608, // 50%
    0.5*533.608 + 0.5*563.263,
    563.263, // 49%
    0.5*563.263 + 0.5*594.07,
    594.07,  // 48%
    0.5*594.07 + 0.5*626.047,
    626.047, // 47%
    0.5*626.047 + 0.5*659.269,
    659.269, // 46%
    0.5*659.269 + 0.5*693.606,
    693.606, // 45%
    0.5*693.606 + 0.5*729.251,
    729.251, // 44%
    0.5*729.251 + 0.5*766.305,
    766.305, // 43%
    0.5*766.305 + 0.5*804.607,
    804.607, // 42%
    0.5*804.607 + 0.5*844.192,
    844.192, // 41%
    0.5*844.192 + 0.5*885.172,
    885.172, // 40%
    0.5*885.172 + 0.5*927.582,
    927.582, // 39%
    0.5*927.582 + 0.5*971.487,
    971.487, // 38%
    0.5*971.487 + 0.5*1016.86,
    1016.86,  // 37%
    0.5*1016.86 + 0.5*1063.71,
    1063.71,  // 36%
    0.5*1063.71 + 0.5*1112.2,
    1112.2,   // 35%
    0.5*1112.2 + 0.5*1162.33,
    1162.33,  // 34%
    0.5*1162.33 + 0.5*1213.91,
    1213.91,  // 33%
    0.5*1213.91 + 0.5*1267.07,
    1267.07,  // 32%
    0.5*1267.07 + 0.5*1322.12,
    1322.12,  // 31%
    0.5*1322.12 + 0.5*1378.92,
    1378.92,  // 30%
    0.5*1378.92 + 0.5*1437.29,
    1437.29,  // 29%
    0.5*1437.29 + 0.5*1497.54,
    1497.54,  // 28%
    0.5*1497.54 + 0.5*1560.05,
    1560.05,  // 27%
    0.5*1560.05 + 0.5*1624.34,
    1624.34,  // 26%
    0.5*1624.34 + 0.5*1690.47,
    1690.47,  // 25%
    0.5*1690.47 + 0.5*1759.06,
    1759.06,  // 24%
    0.5*1759.06 + 0.5*1829.74,
    1829.74,  // 23%
    0.5*1829.74 + 0.5*1902.73,
    1902.73,  // 22%
    0.5*1902.73 + 0.5*1978.02,
    1978.02,  // 21%
    0.5*1978.02 + 0.5*2055.77,
    2055.77,  // 20%
    0.5*2055.77 + 0.5*2136.17,
    2136.17,  // 19%
    0.5*2136.17 + 0.5*2218.88,
    2218.88,  // 18%
    0.5*2218.88 + 0.5*2304.47,
    2304.47,  // 17%
    0.5*2304.47 + 0.5*2393.11,
    2393.11,  // 16%
    0.5*2393.11 + 0.5*2484.75,
    2484.75,  // 15%
    0.5*2484.75 + 0.5*2579.56,
    2579.56,  // 14%
    0.5*2579.56 + 0.5*2677.6,
    2677.6,   // 13%
    0.5*2677.6 + 0.5*2779.68,
    2779.68,  // 12%
    0.5*2779.68 + 0.5*2885.71,
    2885.71,  // 11%
    0.5*2885.71 + 0.5*2995.94,
    2995.94,  // 10%
    0.75*2995.94+0.25*3110.27, // middle bin
    0.50*2995.94+0.50*3110.27, // middle bin
    0.25*2995.94+0.75*3110.27, // middle bin
    3110.27,  //  9%
    0.75*3110.27+0.25*3229.67, // middle bin
    0.50*3110.27+0.50*3229.67, // middle bin
    0.25*3110.27+0.75*3229.67, // middle bin
    3229.67,  //  8%
    0.75*3229.67+0.25*3354.66, // middle bin
    0.50*3229.67+0.50*3354.66, // middle bin
    0.25*3229.67+0.75*3354.66, // middle bin
    3354.66,  //  7%
    0.75*3354.66+0.25*3485.57, // middle bin
    0.50*3354.66+0.50*3485.57, // middle bin
    0.25*3354.66+0.75*3485.57, // middle bin
    3485.57,  //  6%
    0.75*3485.57+0.25*3622.6,   // middle bin
    0.50*3485.57+0.50*3622.6,   // middle bin
    0.25*3485.57+0.75*3622.6,   // middle bin
    3622.6,   //  5%
    0.75*3622.6+0.25*3767.,     // middle bin
    0.50*3622.6+0.50*3767.,     // middle bin
    0.25*3622.6+0.75*3767.,     // middle bin
    3767.,    //  4%
    0.75*3767.+0.25*3920.41,    // middle bin
    0.50*3767.+0.50*3920.41,    // middle bin
    0.25*3767.+0.75*3920.41,    // middle bin
    3920.41,  //  3%
    0.75*3920.41+0.25*4083.38,  // middle bin
    0.50*3920.41+0.50*4083.38,  // middle bin
    0.25*3920.41+0.75*4083.38,  // middle bin
    4083.38,  //  2%
    0.75*4083.38+0.25*4263.72,  // middle bin
    0.50*4083.38+0.50*4263.72,  // middle bin
    0.25*4083.38+0.75*4263.72,  // middle bin
    4263.72,  //  1%
    0.9*4263.72+0.1*5000,     // middle bin
    0.8*4263.72+0.2*5000,     // middle bin
    0.7*4263.72+0.3*5000,     // middle bin
    0.6*4263.72+0.4*5000,     // middle bin
    0.5*4263.72+0.5*5000,     // middle bin
    0.4*4263.72+0.6*5000,     // middle bin
    0.3*4263.72+0.7*5000,     // middle bin
    0.2*4263.72+0.8*5000,     // middle bin
    0.1*4263.72+0.9*5000,     // middle bin
    5000      //  0%,   entry in array is numSuperFineCentBins-1
};

const short numSuperFineCentBins = sizeof (superFineCentBins) / sizeof (superFineCentBins[0]);

void PrintSuperFineCentBins () {
  cout << "\\multicolumn{8}{|c|}{$\\sumet^\\mathrm{FCal}$ bin edges (TeV)} \\\\ \\hline" << endl;
  int iSFCBin = 0;
  for (; iSFCBin < numSuperFineCentBins; iSFCBin++) {
    if (iSFCBin % 8 == 0) cout << "\t\t\t";
    cout << superFineCentBins[iSFCBin]*1e-3;
    if (iSFCBin % 8 == 7) cout << " \\\\" << endl;
    else cout << " & ";
  }
  while (iSFCBin % 8 != 7) {
    cout << "& ";
    iSFCBin++;
  }
  cout << "\\\\ \\hline" << endl;
}

// Returns which "super fine" centrality bin number corresponds to this sum fcal et.
short GetSuperFineCentBin (const float fcal_et) {
  short i = 0;
  while (i < numSuperFineCentBins) {
    if (fcal_et < superFineCentBins[i])
      break;
    i++;
  }
  return i;
}


// Set of percentile-space centrality cuts, intended for plotting # events vs. centrality %
const float percentileCentBins[81]={//164] = {
    5000,     //  0%
    4263.72,  //  1%
    4083.38,  //  2%
    3920.41,  //  3%
    3767.,    //  4%
    3622.6,   //  5%
    3485.57,  //  6%
    3354.66,  //  7%
    3229.67,  //  8%
    3110.27,  //  9%
    2995.94,  // 10%
    2885.71,  // 11%
    2779.68,  // 12%
    2677.6,   // 13%
    2579.56,  // 14%
    2484.75,  // 15%
    2393.11,  // 16%
    2304.47,  // 17%
    2218.88,  // 18%
    2136.17,  // 19%
    2055.77,  // 20%
    1978.02,  // 21%
    1902.73,  // 22%
    1829.74,  // 23%
    1759.06,  // 24%
    1690.47,  // 25%
    1624.34,  // 26%
    1560.05,  // 27%
    1497.54,  // 28%
    1437.29,  // 29%
    1378.92,  // 30%
    1322.12,  // 31%
    1267.07,  // 32%
    1213.91,  // 33%
    1162.33,  // 34%
    1112.2,   // 35%
    1063.71,  // 36%
    1016.86,  // 37%
    971.487,  // 38%
    927.582,  // 39%
    885.172,  // 40%
    844.192,  // 41%
    804.607,  // 42%
    766.305,  // 43%
    729.251,  // 44%
    693.606,  // 45%
    659.269,  // 46%
    626.047,  // 47%
    594.07,   // 48%
    563.263,  // 49%
    533.608,  // 50%
    505.085,  // 51%
    477.734,  // 52%
    451.509,  // 53%
    426.354,  // 54%
    402.144,  // 55%
    378.968,  // 56%
    356.885,  // 57%
    335.738,  // 58%
    315.523,  // 59%
    296.17,   // 60%
    277.742,  // 61%
    260.219,  // 62%
    243.588,  // 63%
    227.751,  // 64%
    212.684,  // 65%
    198.428,  // 66%
    184.922,  // 67%
    172.097,  // 68%
    160.033,  // 69%
    148.625,  // 70%
    137.874,  // 71%
    127.745,  // 72%
    118.249,  // 73%
    109.333,  // 74%
    100.928,  // 75%
    93.071,   // 76%
    85.729,   // 77%
    78.834,   // 78%
    72.411,   // 79%
    66.402    // 80%, entry in array is numPercentileCentBins-1
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
  return i-1;
}


const short numppCentBins = 150;
const double* ppCentBins = linspace (-100, 200, numppCentBins);
int GetPPCentBin (const float fcal_et) {
  int i = 0;
  while (i < numppCentBins) {
    if (fcal_et < ppCentBins[i+1])
      break;
    i++;
  }
  return i;
}


const short numQ2Bins = 4;
const double* q2Bins = linspace (0, 0.2, numQ2Bins);

short GetQ2Bin (const float q2) {
  short i = 0;
  while (i < numQ2Bins) {
    if (q2 < q2Bins[i+1])
      break;
    i++;
  }
  return i;
}

const short numPsi2Bins = 16;
const double* psi2Bins = linspace (-pi/2., pi/2., numPsi2Bins);

//const double psi2Bins[13] = {
//  -pi/2.,
//  -5.*pi/12.,
//  -pi/3.,
//  -pi/4.,
//  -pi/6.,
//  -pi/12.,
//  0,
//  pi/12.,
//  pi/6.,
//  pi/4.,
//  pi/3.,
//  5.*pi/12.,
//  pi/2.
//};
//const short numPsi2Bins = sizeof (psi2Bins) / sizeof (psi2Bins[0]) - 1;

short GetPsi2Bin (const float psi2) {
  short i = 0;
  while (i < numPsi2Bins) {
    if (psi2 < psi2Bins[i+1])
      break;
    i++;
  }
  return i;
}

const int numVZBins = 15;
const double* vzBins = linspace (-150, 150, numVZBins);
short GetVZBin (const float vz) {
  short i = 0;
  while (i < numVZBins) {
    if (vz < vzBins[i+1])
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

const map <string, const set<int>*> runGroups = {
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


void MakeBoxOutline (const double x, const double y, const Color_t color, const double boxMultiplierX = 1., const double boxMultiplierY = 1.) {
  const double ytsize = 0.07;
  const double xtsize = 0.18;
  const double y1 = y - 0.25*ytsize*boxMultiplierY;
  const double y2 = y + 0.25*ytsize*boxMultiplierY;
  const double x2 = x - 0.15*xtsize;
  const double x1 = x - 0.55*xtsize*boxMultiplierX;

  TLine mline;
  mline.SetLineWidth (1);
  mline.SetLineColor (color);
  //mline.SetLineStyle (lstyle);
  mline.SetLineStyle (0);
  //Double_t y_new = (y1+y2)/2.;
  //mline.DrawLineNDC (x1, y_new, x2, y_new);
  mline.DrawLineNDC (x1, y1, x2, y1);
  mline.DrawLineNDC (x1, y2, x2, y2);
  mline.DrawLineNDC (x1, y1, x1, y2);
  mline.DrawLineNDC (x2, y1, x2, y2);
}


void MakeTheoryBox (const double x, const double y, const Color_t color, const double colorAlpha, const double boxMultiplierX = 1., const double boxMultiplierY = 1.) {
  const double ytsize = 0.07;
  const double xtsize = 0.18;
  const double y1 = y - 0.25*ytsize*boxMultiplierY;
  const double y2 = y + 0.25*ytsize*boxMultiplierY;
  const double x2 = x - 0.15*xtsize;
  const double x1 = x - 0.55*xtsize*boxMultiplierX;
  TPave *mbox = new TPave (x1, y1, x2, y2, 0, "NDC");
  mbox->SetFillColorAlpha (color, colorAlpha);
  mbox->SetFillStyle (1001);
  mbox->Draw ();

  //TLine mline;
  //mline.SetLineWidth (1);
  //mline.SetLineColor (color);
  ////mline.SetLineStyle (lstyle);
  //mline.SetLineStyle (0);
  ////Double_t y_new = (y1+y2)/2.;
  ////mline.DrawLineNDC (x1, y_new, x2, y_new);
  //mline.DrawLineNDC (x1, y1, x2, y1);
  //mline.DrawLineNDC (x1, y2, x2, y2);
  //mline.DrawLineNDC (x1, y1, x1, y2);
  //mline.DrawLineNDC (x2, y1, x2, y2);
  return;
}


void MakeTheoryLine (const double x, const double y, const Color_t color, const double lineMultiplier = 1.) {
  const double ytsize = 0.07;
  const double xtsize = 0.18;
  const double x2 = x - 0.15*xtsize;
  const double x1 = x - 0.55*xtsize*lineMultiplier;
  TLine mline;
  mline.SetLineColor (color);
  mline.SetLineWidth (4);
  mline.DrawLineNDC (x1, y, x2, y);
  return;
}


void MakeDataBox (const double x, const double y, const Color_t color, const double colorAlpha, const Style_t mstyle, const double msize, const double bmx = 1., const double bmy = 1.) {
  MakeTheoryBox (x, y, color, colorAlpha, bmx, bmy);
  MakeBoxOutline (x, y, color, bmx, bmy);

  const double ytsize = 0.07;
  const double xtsize = 0.18;

  const double y1 = y - 0.25*ytsize;
  const double y2 = y + 0.25*ytsize;
  const double x2 = x - 0.15*xtsize;
  const double x1 = x - 0.55*xtsize*bmx;

  TLine* ml = new TLine ();
  ml->SetNDC();
  ml->SetLineColor (color);
  ml->SetLineStyle (1);
  ml->SetLineWidth (2);

  ml->DrawLineNDC (0.9*x1+0.1*x2, 0.5*(y1+y2), 0.1*x1+0.9*x2, 0.5*(y1+y2));
  ml->DrawLineNDC (0.5*(x1+x2), 0.9*y1+0.1*y2, 0.5*(x1+x2), 0.1*y1+0.9*y2);

  TMarker* marker = new TMarker (0.5*(x1+x2), y, 0);
  marker->SetNDC();
  marker->SetMarkerColor (color);
  marker->SetMarkerStyle (mstyle);
  marker->SetMarkerSize (msize);
  marker->Draw ();

  if (IsFullMarker (mstyle)) {
    TMarker* marker2 = new TMarker (0.5*(x1+x2), y, 0);
    marker2->SetNDC();
    marker2->SetMarkerColor (kBlack);
    marker2->SetMarkerStyle (FullToOpenMarker (mstyle));
    marker2->SetMarkerSize (msize);
    marker2->Draw();
  }
  return;
}


#endif
