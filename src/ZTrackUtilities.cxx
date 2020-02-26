#ifndef __ZTrackUtilities_cxx__
#define __ZTrackUtilities_cxx__

#include "ZTrackUtilities.h"
#include "Params.h"

#include <TSystemDirectory.h>

#include <math.h>
#include <iostream>
#include <sys/stat.h>

using namespace std;

namespace ZTrackAnalyzer {

/**
 * Establishes path variables appropriately.
 */
void SetupDirectories (const TString dataSubDir) {
  rootPath = extWorkPath + "rootFiles/" + dataSubDir + "/";
  plotPath = workPath + "Plots/" + dataSubDir + "/";

  if (doElectronPtUpVar)
    rootPath = rootPath + "Variations/ElectronPtUpVariation/";
  else if (doElectronPtDownVar)
    rootPath = rootPath + "Variations/ElectronPtDownVariation/";
  else if (doMuonPtUpVar)
    rootPath = rootPath + "Variations/MuonPtUpVariation/";
  else if (doMuonPtDownVar)
    rootPath = rootPath + "Variations/MuonPtDownVariation/";
  else if (doElectronLHMediumVar)
    rootPath = rootPath + "Variations/ElectronLHMediumWPVariation/";
  else if (doMuonTightVar)
    rootPath = rootPath + "Variations/MuonTightWPVariation/";
  else if (doHITightVar)
    rootPath = rootPath + "Variations/TrackHITightWPVariation/";
  else if (doPionsOnlyVar)
    rootPath = rootPath + "Variations/TrackEffPionsVariation/";
  else
    rootPath = rootPath + "Nominal/";
}


/**
 * Returns a copy of the histogram detailing the Zdc cuts.
 */
TH1D* GetZdcCuts () {
  if (!isPbPb)
    return nullptr; // no cuts needed outside of Pb+Pb

  TFile* zdcFile = new TFile (Form ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/ZdcAnalysis/HIRun2PileUp_PbPb5p02_%s.root", is2015data ? "2015" : "2018"), "read");
  TH1D* h_zdcCuts = (TH1D*) ((TH1D*) zdcFile->Get (is2015data ? "hCut" : "h_cuts")->Clone ("h_zdcCuts"));
  h_zdcCuts->SetDirectory (0);
  zdcFile->Close ();
  SaferDelete (zdcFile);
  return h_zdcCuts;
}


/**
 * Returns the appropriate file in the given directory.
 * For MC, inFileName MUST be specified.
 * If dataSet == 0, will assume this is a "PhysCont" sample, i.e. an entire collision system in 1 data set.
 */
TFile* GetFile (const char* directory, const int dataSet, const char* inFileName) {
  TFile* file = nullptr;

  // First figure out the file we are looking for
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (!isMC) {
      if (dataSet == 0)
        fileIdentifier = "PhysCont.AOD.";
      else
        fileIdentifier = to_string (dataSet);
    }
    else {
      cout << "Error: In ZTrackUtilities.C: Cannot identify this MC file! Will return null!" << endl;
      return nullptr;
    }
  }
  else
    fileIdentifier = inFileName;

  // Now get the list of files
  const TString dataPathTemp = dataPath + "/" + directory + "/";
  TSystemDirectory dir (dataPathTemp.Data (), dataPathTemp.Data ());
  TList* sysfiles = dir.GetListOfFiles ();
  if (!sysfiles) {
    cout << "Error: In ZTrackUtilities.C: Cannot get list of files! Will return null!" << endl;
    return nullptr;
  }
  TSystemFile* sysfile;
  TString fname;
  TIter next (sysfiles);

  while ((sysfile = (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
      if (fname.Contains (fileIdentifier)) {
        break;
      }
    }
  }

  if (!fname.Contains (fileIdentifier)) {
    cout << "Error: In ZTrackUtilities.cxx: TFile cannot be found for given data set. Will return null!" << endl;
    return nullptr;
  }

  const TString sourceName = dataPathTemp + fname;

  if (useScratchDisk) {
    TString condorScratchDiskPath = TString (getenv ("_CONDOR_SCRATCH_DIR")) + "/jeff/";

    const TString fullName = condorScratchDiskPath + fname;

    struct stat info;

    if (stat (fullName, &info) == 0) { // if the file can be resolved, use it
      cout << "Resolved file on scratch disk: " << fullName << endl;
    }
    else { // otherwise we need to copy to scratch disk
      cout << "Copying file from pnfs: " << sourceName.Data () << endl << "destination: " << fullName.Data () << endl;
      system (Form ("mkdir -p %s", condorScratchDiskPath.Data ()));
      //system (Form ("rm %s", fullName.Data ()));
      if (dataPath.Contains ("pnfs")) // pnfs disk needs special copy command
        system (Form ("xrdcp -f root://dcxrd.usatlas.bnl.gov:1096//%s %s", sourceName.Data (), fullName.Data ()));
      else // otherwise cp command works
        system (Form ("cp %s %s", sourceName.Data (), fullName.Data ()));
      cout << "Finished copying file: " << fullName.Data () << endl;
    }

    file = new TFile (fullName, "read");
  }
  else {
    cout << "Resolved file at " << sourceName << endl;
    file = new TFile (sourceName, "read");
  }
  
  if (!file) {
    cout << "Error: In ZTrackUtilities.C: TFile not obtained for given data set. Will return null!" << endl;
    return nullptr;
  }
  else return file;
} 


/**
 * Returns an abbreviated, unique identifier for a given dataset.
 */
TString GetIdentifier (const int dataSet, const char* directory, const char* inFileName) {
  TString id;
  if (!isMC) {
    if (dataSet == 0)
      id = "data15hi";
    else
      id = to_string (dataSet);
    if (TString (directory).Contains ("pc"))
      id += "_pc";
    else if (TString (directory).Contains ("cc"))
      id += "_cc";
    return id;
  }

  id = (isPbPb ? (is2015data ? "PbPb15" : "PbPb18") : "pp");

  if (TString (inFileName).Contains ("Sherpa"))
    id = id + "_Sherpa";

  if (isPbPb) {
    if (TString (inFileName).Contains ("pp_Zee"))
      id = id + "_pp_Zee";
    else if (TString (inFileName).Contains ("pn_Zee"))
      id = id + "_pn_Zee";
    else if (TString (inFileName).Contains ("np_Zee"))
      id = id + "_np_Zee";
    else if (TString (inFileName).Contains ("nn_Zee"))
      id = id + "_nn_Zee";
    else if (TString (inFileName).Contains ("Zee"))
      id = id + "_pp_Zee";
    else if (TString (inFileName).Contains ("pp_Zmumu"))
      id = id + "_pp_Zmumu";
    else if (TString (inFileName).Contains ("pn_Zmumu"))
      id = id + "_pn_Zmumu";
    else if (TString (inFileName).Contains ("np_Zmumu"))
      id = id + "_np_Zmumu";
    else if (TString (inFileName).Contains ("nn_Zmumu"))
      id = id + "_nn_Zmumu";
    else if (TString (inFileName).Contains ("Zmumu"))
      id = id + "_pp_Zmumu";
    else if (TString (inFileName).Contains ("Ztautau"))
      id = id + "_pp_Ztautau";
    else if (TString (inFileName).Contains ("ttbar"))
      id = id + "_pp_ttbar";
  } 
  else {
    if (TString (inFileName).Contains ("Zee"))
      id = id + "_Zee";
    else if (TString (inFileName).Contains ("Zmumu"))
      id = id + "_Zmumu";
    else if (TString (inFileName).Contains ("Ztautau"))
      id = id + "_Ztautau";
    else if (TString (inFileName).Contains ("ttbar"))
      id = id + "_ttbar";
  }
  

  if (TString (inFileName).Contains ("Hijing")) {
    id = id + "_Hijing";
    //if (collisionSystem == pp15 || collisionSystem == PbPb15)
    //  id += "_15";
    //else if (collisionSystem == pPb16 || collisionSystem == Pbp16)
    //  id += "_16";
    //else if (collisionSystem == XeXe17 || collisionSystem == pp17)
    //  id += "_17";
    //else if (collisionSystem == PbPb18)
    //  id += "_18";
  }

  if (TString (inFileName).Contains ("ptmin25")) {
    id = id + "_ptmin25";
  }
  else if (TString (inFileName).Contains ("ygt175")) {
    id = id + "_ygt175";
  }

  return id;
}


/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
double* linspace (double lo, double hi, int num) {
  double* arr = new double[num+1];
  double delta = ((double)(hi)-(double)(lo))/(double)(num);
  for (int i = 0; i <= num; i++) {
    arr[i] = lo + i * delta;
  }
  return arr;
}


/**
 * Returns a logarithmically spaced array, where the 0th element is lo and the num-th element is hi.
 */
double* logspace (double lo, double hi, int num) {
  double loghi = log2(hi);
  if (lo == 0) {
    double* arr = linspace(log2(hi/(100*num)), loghi, num);
    for (int i = 0; i <= num; i++) {
      arr[i] = pow(2, arr[i]);
    }
    return arr;
  } else {
    double loglo = log2(lo);
    double* arr = linspace(loglo, loghi, num);
    for (int i = 0; i <= num; i++) {
      arr[i] = pow(2, arr[i]);
    }
    return arr;
  }
}


/**
 * Returns the equivalent angle in the range 0 to 2pi.
 */
double InTwoPi (double phi) {
  while (phi < 0 || 2*pi <= phi) {
   if (phi < 0) phi += 2*pi;
   else phi -= 2*pi;
  }
  return phi;
}


/**
 * Returns the difference between two angles in 0 to pi.
 * If sign is true, will return a signed version such that phi2 = phi1 + dphi
 */
double DeltaPhi (double phi1, double phi2, const bool sign) {
  phi1 = InTwoPi(phi1);
  phi2 = InTwoPi(phi2);
  double dphi = abs(phi1 - phi2);
  while (dphi > pi) dphi = abs (dphi - 2*pi);

  if (sign && InTwoPi (phi2 + dphi) == phi1)
     dphi *= -1;

  return dphi;
}


/**
 * Returns dR between two eta, phi coordinates.
 */
double DeltaR (const double eta1, const double eta2, const double phi1, const double phi2 ) {
 const double deta = eta1 - eta2;
 const double dphi = DeltaPhi (phi1, phi2, false);
 return sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
}


/**
 * Returns true iff this eta, phi coordinate lies in the disabled HEC region.
 */
bool InDisabledHEC (const double eta, double phi, const double dr) {
  phi = InTwoPi (phi);
  return 1.5-dr < eta && eta < 3.2+dr &&
         pi-dr < phi && phi < 3*pi/2+dr;
}


/**
 * Returns true iff this eta lies within the EMCal.
 */
bool InEMCal (const float eta) {
  return fabs(eta) < 1.37 || (1.52 < fabs(eta) && fabs(eta) < 2.47);
}


/**
 * Returns true iff this object is within a given radius in the HCal.
 */
bool InHadCal (const float eta, const float R) {
  return fabs(eta) <= 4.9 - R;
}


/**
 * Returns the energy scale correction factor for this electron as derived in MC between pp and Pb+Pb.
 */
float GetZmassSF_MC (const float fcal_et, const float eta) {
  if (fcal_et < 66.402)
    return 1;
  else if (66.402 <= fcal_et && fcal_et < 1378.92)
    return (fabs (eta) > 1 ? 1.0162 : 1.00327);
  else if (1378.92 <= fcal_et && fcal_et < 2995.94)
    return (fabs (eta) > 1 ? 1.02528 : 1.01058);
  else if (2995.94 <= fcal_et && fcal_et < 5000) 
    return (fabs (eta) > 1 ? 1.03133 : 1.01624);
  else
    return 1;
}


/**
 * Returns the energy scale correction factor for the electron as derived in Pb+Pb between data and MC.
 */
float GetZmassSF_PbPb (const float fcal_et, const float eta) {
  if (fcal_et < 66.402)
    return 1;
  else if (66.402 <= fcal_et && fcal_et < 1378.92)
    return (fabs (eta) > 1 ? 0.991796 : 0.996473);
  else if (1378.92 <= fcal_et && fcal_et < 2995.94)
    return (fabs (eta) > 1 ? 0.991246 : 1.00095);
  else if (2995.94 <= fcal_et && fcal_et < 5000) 
    return (fabs (eta) > 1 ? 0.993157 : 0.995561);
  else
    return 1;
}


} // end namespace

#endif
