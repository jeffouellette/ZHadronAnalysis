#ifndef __ZTrackUtilities_h__
#define __ZTrackUtilities_h__


#include <TString.h>
#include <TH1D.h>
#include <TFile.h>

namespace ZTrackAnalyzer {

/**
 * Establishes path variables appropriately.
 */
void SetupDirectories (const TString dataSubDir, const bool addSubDir = true);


/**
 * Returns a copy of the histogram detailing the Zdc cuts.
 */
TH1D* GetZdcCuts ();


/**
 * Returns the appropriate file in the given directory.
 * For MC, inFileName MUST be specified.
 */
TFile* GetFile (const char* directory, const int dataSet, const char* inFileName);


/**
 * Returns an abbreviated, unique identifier for a given dataset.
 */
TString GetIdentifier (const int dataSet, const char* directory, const char* inFileName);


/**
 * Safely deletes pointer T.
 */
template <typename T> inline void SaferDelete (T** t) {
  if (*t) { delete (*t); (*t) = nullptr; }
}



/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
double* linspace (double lo, double hi, int num);


/**
 * Returns a logarithmically spaced array, where the 0th element is lo and the num-th element is hi.
 */
double* logspace (double lo, double hi, int num);


/**
 * Returns the equivalent angle in the range 0 to 2pi.
 */
double InTwoPi (double phi);


/**
 * Returns the difference between two angles in 0 to pi.
 * If sign is true, will return a signed version such that phi2 = phi1 + dphi
 */
double DeltaPhi (double phi1, double phi2, const bool sign = false);


/**
 * Returns dR between two eta, phi coordinates.
 */
double DeltaR (const double eta1, const double eta2, const double phi1, const double phi2 );


/**
 * Returns true iff this eta, phi coordinate lies in the disabled HEC region.
 */
bool InDisabledHEC (const double eta, double phi, const double dr = 0.4);


/**
 * Returns true iff this eta lies within the EMCal.
 */
bool InEMCal (const float eta);


/**
 * Returns true iff this object is within a given radius in the HCal.
 */
bool InHadCal (const float eta, const float R = 0.4);


/**
 * Returns the energy scale correction factor for this electron as derived in MC between pp and Pb+Pb.
 */
float GetZmassSF_MC (const float fcal_et, const float eta);


/**
 * Returns the energy scale correction factor for the electron as derived in Pb+Pb between data and MC.
 */
float GetZmassSF_PbPb (const float fcal_et, const float eta);

} // end namespace

#endif
