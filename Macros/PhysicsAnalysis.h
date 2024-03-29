#ifndef __PhysicsAnalysis_h__
#define __PhysicsAnalysis_h__

#include "Params.h"
#include "EventPlaneCalibrator.h"
//#include "CovarianceFit.C"

#include <ArrayTemplates.h>

#include <TEfficiency.h>
#include <TClass.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TBox.h>
#include <TVirtualFitter.h>
#include <TLatex.h>

#include <iostream>
#include <string>

using namespace std;

typedef TGraphAsymmErrors TGAE;

class PhysicsAnalysis {

  protected:
  string name = "";
  bool backgroundSubtracted = false;
  bool sameSignsSubtracted = false;
  bool hasBkg = true;

  vector<TH1*> drawnHists;
  vector<TGAE*> drawnGraphs;

  bool iaaCalculated = false;
  //bool icpCalculated = false;
  bool sigToBkgCalculated = false;
  bool trackMeansCalculated = false;

  TFile* eventWeightsFile = nullptr;
  bool eventWeightsLoaded = false;
  TFile* trkEffFile = nullptr;
  bool effsLoaded   = false;
  TFile* trkPurFile = nullptr;
  bool pursLoaded   = false;
  TFile* histFile   = nullptr;
  bool histsLoaded  = false;
  bool histsScaled  = false;
  bool histsReflected = false;

  public:
  bool histsUnfolded = false;
  bool plotFill       = false; // whether to plot as filled (bar) graph or points w/ errors
  bool plotSignal     = true; // whether to plot background subtracted plots
  bool useAltMarker   = false; // whether to plot as open markers (instead of closed)
  bool plotAsSyst     = false; // whether to plot as a systematic uncertainty (box uncertainties)

  bool isMC           = false;
  bool isBkg          = false; // whether analysis represents a background (UE) track yield.
  bool useImpactParameter = false; // whether to use impact parameter mixing (instead of FCal Sum Et) -- only applicable for Hijing
  bool doPPMBMixing   = false; // whether analysis uses Minimum Bias mixing in pp collisions
  bool subtractPP     = true;  // whether analysis should subtract pp background
  bool unfoldPbPb     = true;  // whether Z events are reconstructed or truth-selected in PbPb
  bool useCentWgts    = false; // whether to reweight this analysis to the data centrality distribution (important for getting correct tracking efficiencies)
  bool useQ2Wgts      = false; // whether to reweight this analysis to the data |q2| distribution
  bool usePsi2Wgts    = false; // whether to reweight this analysis to the data psi2 distribution
  bool is2015Conds    = false; // whether this analysis uses 2015 data (different conditions)
  bool useHITight     = false; // whether to use HITight tracking efficiencies
  bool useHijingEffs  = false; // whether to use tracking efficiencies derived from Hijing
  bool doLeptonRejVar = false; // whether to impose an additional dR cut on tracks away from the leptons
  bool doTrackPurVar  = false; // whether to impose an additional correction based on tracking purity
  bool doTrackEffVar  = false; // whether to use pions-only tracking efficiency variation
  float trkEffNSigma  = 0; // how many sigma to vary the track efficiency by (-1,0,+1 suggested)
  float trkPurNSigma  = 0; // how many sigma to vary the track purity by (-1,0,+1 suggested)

  TString eventWeightsFileName = "DataAnalysis/Nominal/eventWeightsFile.root";
  TString ewExt = "";

  // Reference to background histograms
  PhysicsAnalysis* bkg = nullptr;

  // Analysis checks
  TH1D*   h_fcal_et               = nullptr;
  TH1D*   h_fcal_et_reweighted    = nullptr;
  TH1D**  h_centrality            = Get1DArray <TH1D*> (3);
  TH1D**  h_centrality_reweighted = Get1DArray <TH1D*> (3);

  TH1D**  h_q2                  = Get1DArray <TH1D*> (numFineCentBins);
  TH1D**  h_q2_reweighted       = Get1DArray <TH1D*> (numFineCentBins);
  TH1D**  h_psi2                = Get1DArray <TH1D*> (numFineCentBins);
  TH1D**  h_psi2_reweighted     = Get1DArray <TH1D*> (numFineCentBins);
  TH1D*   h_PbPb_vz             = nullptr;
  TH1D*   h_PbPb_vz_reweighted  = nullptr;
  TH1D*   h_pp_vz               = nullptr;
  TH1D*   h_pp_vz_reweighted    = nullptr;
  TH1D*   h_pp_nch              = nullptr;
  TH1D*   h_pp_nch_reweighted   = nullptr;

  // Event info distributions (for reweighting)
  TH1D*** h_PbPbFCal_weights   = Get2DArray <TH1D*> (3, nPtZBins+1);                  // iSpc, iPtZ
  TH1D**** h_PbPbQ2_weights    = Get3DArray <TH1D*> (3, numFineCentBins, nPtZBins+1); // iSpc, iFineCent, iPtZ
  TH1D**** h_PbPbPsi2_weights  = Get3DArray <TH1D*> (3, numFineCentBins, nPtZBins+1); // iSpc, iFineCent, iPtZ

  // Event plane calibration information
  EventPlaneCalibrator eventPlaneCalibrator;

  // Efficiencies
  TH1D*** h_trk_effs          = Get2DArray <TH1D*> (numTrkCorrCentBins, numEtaTrkBins); // iCent, iEta
  //TF1**   f_trk_effs          = Get1DArray <TH1D*> (numTrkCorrCentBins); // iCent (eta dependence is extrapolated)
  TH2D**  h2_trk_effs         = Get1DArray <TH2D*> (numTrkCorrCentBins); // iCent
  //TEfficiency*** h_trk_effs   = Get2DArray <TEfficiency*> (numTrkCorrCentBins, numEtaTrkBins); // iCent, iEta
  TH2D** h2_num_trk_effs      = Get1DArray <TH2D*> (numTrkCorrCentBins); // iCent
  TH2D** h2_den_trk_effs      = Get1DArray <TH2D*> (numTrkCorrCentBins);

  // Tracking purities
  TH1D*** h_trk_purs          = Get2DArray <TH1D*> (numTrkCorrCentBins, numEtaTrkBins); // iCent, iEta
  TH2D**  h2_trk_purs         = Get1DArray <TH2D*> (numTrkCorrCentBins); // iCent
  //TEfficiency*** h_trk_purs   = Get2DArray <TEfficiency*> (numTrkCorrCentBins, numEtaTrkBins); // iCent, iEta
  TH2D** h2_num_trk_purs      = Get1DArray <TH2D*> (numTrkCorrCentBins); // iCent
  TH2D** h2_den_trk_purs      = Get1DArray <TH2D*> (numTrkCorrCentBins);

  // Bin migration factors
  TF1**** f_trk_pt_ptz_binMigration  = Get3DArray <TF1*> (2, nPtZBins, numBBBCorrCentBins); // iSpc, iPtZ, iCent
  TF1**** f_trk_xhz_ptz_binMigration = Get3DArray <TF1*> (2, nPtZBins, numBBBCorrCentBins); // iSpc, iPtZ, iCent

  // Physics plots
  TH1D*****   h_trk_dphi              = Get4DArray <TH1D*> (3, nPtZBins, maxNPtchBins+2, numCentBins+1); // iSpc, iPtZ, iPtch (+1 for integrated), iCent (+1 for 0-30% integrated)
  TH2D*****   h2_trk_dphi_cov         = Get4DArray <TH2D*> (3, nPtZBins, maxNPtchBins, numCentBins); // iSpc, iPtZ, iPtch, iCent
  TH1D*****   h_trk_dphi_sub          = Get4DArray <TH1D*> (3, nPtZBins, maxNPtchBins+2, numCentBins+1); // iSpc, iPtZ, iPtch (+1 for integrated), iCent (+1 for 0-30% integrated)
  TH1D*****   h_trk_dphi_sig_to_bkg   = Get4DArray <TH1D*> (3, nPtZBins, maxNPtchBins+2, numCentBins+1); // iSpc, iPtZ, iPtch (+1 for integrated), iCent (+1 for 0-30% integrated)

  TH1D*****   h_trk_pt_dphi_raw = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins);   // iSpc, iPtZ, iPhi, iCent
  TH1D****    h_z_counts        = Get3DArray <TH1D*> (3, nPtZBins, numCentBins);               // iSpc, iPtZ, iCent

  TH1D*****  h_trk_pt_dphi            = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH2D*****  h2_trk_pt_dphi_cov       = Get4DArray <TH2D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_trk_pt_dphi_sub        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_trk_pt_dphi_sig_to_bkg = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D*****  h_trk_pt_dphi_iaa        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  //TH1D*****  h_trk_pt_dphi_icp        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent

  TH1D****   h_trk_pt_ptz               = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH2D****   h2_trk_pt_ptz_cov          = Get3DArray <TH2D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_pt_ptz_sub           = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_pt_ptz_sig_to_bkg    = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_pt_ptz_iaa           = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D***    h_trk_pt_ptz_iaa_2015proj  = Get2DArray <TH1D*> (nPtZBins, numCentBins); // iPtZ, iCent
  //TH1D****   h_trk_pt_ptz_icp           = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TF1****    f_trk_fit_pt_ptz           = Get3DArray <TF1*>  (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TGAE***    g_trk_avg_pt_ptz           = Get2DArray <TGAE*> (3, numCentBins); // iSpc, iCent

  TH1D***** h_trk_xhz_dphi            = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH2D***** h2_trk_xhz_dphi_cov       = Get4DArray <TH2D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_trk_xhz_dphi_sub        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_trk_xhz_dphi_sig_to_bkg = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  TH1D***** h_trk_xhz_dphi_iaa        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  //TH1D***** h_trk_xhz_dphi_icp        = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent

  TH1D****   h_trk_xhz_ptz              = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH2D****   h2_trk_xhz_ptz_cov         = Get3DArray <TH2D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_xhz_ptz_sub          = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_xhz_ptz_sig_to_bkg   = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D****   h_trk_xhz_ptz_iaa          = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D***    h_trk_xhz_ptz_iaa_2015proj = Get2DArray <TH1D*> (nPtZBins, numCentBins); // iPtZ, iCent
  //TH1D****   h_trk_xhz_ptz_icp          = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TF1****    f_trk_fit_xhz_ptz          = Get3DArray <TF1*>  (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TGAE***    g_trk_avg_xhz_ptz          = Get2DArray <TGAE*> (3, numCentBins); // iSpc, iCent

  TGAE***    g_avg_ntrk_ptz             = Get2DArray <TGAE*> (3, numCentBins); // iSpc, iCent

  PhysicsAnalysis () { }

  PhysicsAnalysis (const char* _name) {
    name = _name;
    plotFill = false;
  }

  virtual ~PhysicsAnalysis () {
    Delete1DArray (&h_centrality,            3);
    Delete1DArray (&h_centrality_reweighted, 3);
    Delete1DArray (&h_q2,                    numFineCentBins);
    Delete1DArray (&h_q2_reweighted,         numFineCentBins);
    Delete1DArray (&h_psi2,                  numFineCentBins);
    Delete1DArray (&h_psi2_reweighted,       numFineCentBins);

    Delete2DArray (&h_PbPbFCal_weights,  3, nPtZBins+1);
    Delete3DArray (&h_PbPbQ2_weights,    3, numFineCentBins, nPtZBins+1);
    Delete3DArray (&h_PbPbPsi2_weights,  3, numFineCentBins, nPtZBins+1);

    ClearHists ();

    Delete2DArray (&h_trk_effs,      numTrkCorrCentBins, numEtaTrkBins);
    Delete1DArray (&h2_trk_effs,     numTrkCorrCentBins);
    Delete1DArray (&h2_num_trk_effs, numTrkCorrCentBins);
    Delete1DArray (&h2_den_trk_effs, numTrkCorrCentBins);

    Delete2DArray (&h_trk_purs,      numTrkCorrCentBins, numEtaTrkBins);
    Delete1DArray (&h2_trk_purs,     numTrkCorrCentBins);
    Delete1DArray (&h2_num_trk_purs, numTrkCorrCentBins);
    Delete1DArray (&h2_den_trk_purs, numTrkCorrCentBins);

    Delete3DArray (&f_trk_pt_ptz_binMigration,  2, nPtZBins, numBBBCorrCentBins);
    Delete3DArray (&f_trk_xhz_ptz_binMigration, 2, nPtZBins, numBBBCorrCentBins);

    Delete4DArray (&h_trk_dphi,                3, nPtZBins, maxNPtchBins+2, numCentBins+1);
    Delete4DArray (&h2_trk_dphi_cov,           3, nPtZBins, maxNPtchBins, numCentBins);
    Delete4DArray (&h_trk_dphi_sub,            3, nPtZBins, maxNPtchBins+2, numCentBins+1);

    Delete4DArray (&h_trk_pt_dphi_raw,         3, nPtZBins, numPhiBins, numCentBins);
    Delete3DArray (&h_z_counts,                3, nPtZBins, numCentBins);

    Delete4DArray (&h_trk_pt_dphi,             3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (&h2_trk_pt_dphi_cov,        3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (&h_trk_pt_dphi_sub,         3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (&h_trk_pt_dphi_sig_to_bkg,  3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (&h_trk_pt_dphi_iaa,         3, nPtZBins, numPhiBins, numCentBins);
    //Delete4DArray (&h_trk_pt_dphi_icp,         3, nPtZBins, numPhiBins, numCentBins);

    Delete3DArray (&h_trk_pt_ptz,              3, nPtZBins, numCentBins);
    Delete3DArray (&h2_trk_pt_ptz_cov,         3, nPtZBins, numCentBins);
    Delete3DArray (&h_trk_pt_ptz_sub,          3, nPtZBins, numCentBins);
    Delete3DArray (&h_trk_pt_ptz_sig_to_bkg,   3, nPtZBins, numCentBins);
    Delete3DArray (&h_trk_pt_ptz_iaa,          3, nPtZBins, numCentBins);
    //Delete3DArray (&h_trk_pt_ptz_icp,          3, nPtZBins, numCentBins);
    Delete2DArray (&g_trk_avg_pt_ptz,          3, numCentBins);

    Delete4DArray (&h_trk_xhz_dphi,            3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (&h2_trk_xhz_dphi_cov,       3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (&h_trk_xhz_dphi_sub,        3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (&h_trk_xhz_dphi_sig_to_bkg, 3, nPtZBins, numPhiBins, numCentBins);
    Delete4DArray (&h_trk_xhz_dphi_iaa,        3, nPtZBins, numPhiBins, numCentBins);
    //Delete4DArray (&h_trk_xhz_dphi_icp,        3, nPtZBins, numPhiBins, numCentBins);

    Delete3DArray (&h_trk_xhz_ptz,             3, nPtZBins, numCentBins);
    Delete3DArray (&h2_trk_xhz_ptz_cov,        3, nPtZBins, numCentBins);
    Delete3DArray (&h_trk_xhz_ptz_sub,         3, nPtZBins, numCentBins);
    Delete3DArray (&h_trk_xhz_ptz_sig_to_bkg,  3, nPtZBins, numCentBins);
    Delete3DArray (&h_trk_xhz_ptz_iaa,         3, nPtZBins, numCentBins);
    //Delete3DArray (&h_trk_xhz_ptz_icp,         3, nPtZBins, numCentBins);
    Delete2DArray (&g_trk_avg_xhz_ptz,         3, numCentBins);

    if (eventWeightsFile && eventWeightsFile->IsOpen ()) {
      eventWeightsFile->Close ();
      SaferDelete (&eventWeightsFile);
    }
    if (trkEffFile && trkEffFile->IsOpen ()) {
      trkEffFile->Close ();
      SaferDelete (&trkEffFile);
    }
    if (trkPurFile && trkPurFile->IsOpen ()) {
      trkPurFile->Close ();
      SaferDelete (&trkPurFile);
    }
    if (histFile && histFile->IsOpen ()) {
      histFile->Close ();
      SaferDelete (&histFile);
    }
  }


  protected:
  void LabelCorrelations (const short iPtZ, const short iPtch, const short iCent, const bool subBkg);
  void LabelIAA_dCent (const short iCent, const short iPhi, const short iPtZ);

  void GetDrawnObjects ();
  void GetMinAndMax (double &min, double &max, const bool log = false);
  void SetMinAndMax (double min, double max);

  public:
  string Name () { return name; }
  void SetName (string _name) { name = _name; }

  virtual TGAE* GetTGAE (TH1D* h);

  static void SetErrors (TH1D* h, TH2D* h2);
  static void TruncateTGAE (TGAE* g, const short iPtZ, const bool useTrkPt);
  static void TruncateTH1D (TH1D** _h, const short iPtZ, const bool useTrkPt);

  virtual void ClearHists ();
  virtual void CreateHists ();
  virtual void CopyAnalysis (PhysicsAnalysis* a, const bool copySigs = false, const bool copyIAAs = true);
  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true);
  virtual void SaveHists (const char* histFileName = "savedHists.root");
  virtual void SaveResults (const char* saveFileName = "resultsHists.root");
  virtual void ScaleHists ();
  virtual void SetVariances ();
  virtual void SubtractBackground (PhysicsAnalysis* _bkg = nullptr, const bool addUnc = false);
  virtual void UnfoldSubtractedYield ();
  virtual void CombineHists ();
  virtual void MakeSummaryCorrelations ();
  virtual void CalculateSigToBkg ();
  virtual void CalculateTrackMeans (PhysicsAnalysis* nom, TH1D*** h_zpt_ptr, PhysicsAnalysis* cov = nullptr, const double meanVar = 0);
  virtual void TruncatePhysicsPlots ();
  virtual void ProjectLumiIncrease ();

  virtual void InflateStatUnc (const float amount);
  virtual void ApplyRelativeVariation (float**** relVar, const bool upVar = true); // multiplies yield results by relErr in each bin (or divides if not upVar)
  virtual void ConvertToStatVariation (const bool upVar = true, const float nSigma = 1); // adds or subtracts nSigma of statistical errors to analysis

  virtual void Execute (const char* inFileName, const char* outFileName);

  virtual void GenerateEventWeights (const char* weightedSampleInFileName, const char* matchedSampleInFileName, const char* outFileName);
  virtual void LoadEventWeights ();

  virtual void LoadTrackingEfficiencies (const bool doRebin = false); // defaults to HILoose
  virtual double GetTrackingEfficiency (const float fcal_et, double trk_pt, const float trk_eta, const bool isPbPb = true);

  virtual void LoadTrackingPurities (const bool doRebin = false); // defaults to HILoose
  virtual double GetTrackingPurity (const float fcal_et, double trk_pt, const float trk_eta, const bool isPbPb = true);

  void CorrectQ2Vector (float& q2x_a, float& q2y_a, float& q2x_c, float& q2y_c);

  void PrintZYields (const bool latexForm = false);
  void PrintTrackYields (const bool latexForm = false);

  void PlotFCalDists (const bool _treatAsData = true);
  void PlotQ2Dists (const bool _treatAsData = true);
  void PlotQ2Weights (PhysicsAnalysis* a);
  void PlotPsi2Dists (const bool _treatAsData = true);
  void PlotPsi2Weights (PhysicsAnalysis* a);
  void PlotVZDists (const bool _treatAsData = true);
  void PlotNchDists (const bool _treatAsData = true);

  void PlotCorrelations (const short pSpc = 2, const short pPtZ = nPtZBins-1, const bool _subBkg = false);
  void PlotTrackingEfficiencies ();
  void PlotTrackingEfficienciesComparison (PhysicsAnalysis* a = nullptr);
  void PlotTrackingEfficiencies2D ();
  void PlotTrackingPurities ();
  void PlotTrackingPuritiesComparison (PhysicsAnalysis* a = nullptr);
  void PlotTrackingPurities2D ();

  void CalculateIAA (const bool overwrite = false);

  virtual void PlotUnweightedTrkYields (const bool useTrkPt = true, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotAllYields_dPhi (const bool useTrkPt = true, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotAllYields_dPtZ (const bool useTrkPt = true, const short pSpc = 2);
  virtual void PlotAllYields_dPtZ_SpcComp (const bool useTrkPt = true);
  virtual void PlotPPYields_dPtZ_SpcComp (const bool useTrkPt, const short iPtZ = nPtZBins-1);
  virtual void PlotClosureCheck (PhysicsAnalysis* truthComp, const bool useTrkPt = true, const short iSpc = 2);
  virtual void PlotClosureCheck_SigToBkg (PhysicsAnalysis* truthComp, const bool useTrkPt = true, const short iSpc = 2, const char* saveFileName = "");
  virtual void PlotClosureCheck_SigToBkg_BkgNormalization (PhysicsAnalysis* truthComp, const bool useTrkPt = true, const short iSpc = 2, const char* saveFileName = "");
  virtual void PlotAllCovMatrices ();
  virtual void PlotCovMatrix (const bool useTrkPt = true, const short pSpc = 0, const short pPtZ = nPtZBins-1, const short pCent = numCentBins-1);
  virtual void PlotIAA_dPhi (const bool useTrkPt = true, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotIAA_dCent (const bool useTrkPt = true, const short pSpc = 2, const short pPtZ = nPtZBins-1);
  virtual void PlotIAA_dPtZ (const bool useTrkPt = true, const short pSpc = 2);
  virtual void PlotSingleIAA_dPtZ (const bool useTrkPt = true, const short pPtZ = -1, const short iCent = numCentBins-1, const short pSpc = 2);
  virtual void PlotIAA_dPtZ_SpcComp (const bool useTrkPt = true);
  virtual void PlotSignalToBkg (const bool useTrkPt = true, const short iSpc = 2);
  virtual void PlotTrackMeans (const bool useTrkPt = true, const short iSpc = 2);
  virtual void PlotSubYields_dPtZ_Fits (short pSpc = 2);

  virtual void WriteIAAs ();
  virtual void PrintIAA (const bool printErrs, const bool useTrkPt = true, const short iCent = numCentBins-1, const short iPtZ = nPtZBins-1, const short iSpc = 2);
};


void SafeWrite (TObject* tobj, TString nname = "") {
  if (tobj)
    tobj->Write (nname);
}


void PhysicsAnalysis :: GetDrawnObjects () {
  TList* primitives = gPad->GetListOfPrimitives ();
  drawnHists.clear ();
  drawnGraphs.clear ();
  for (int i = 0; i < primitives->GetSize (); i++) {
    TObject *obj = primitives->At (i);
    if (obj->IsA()->InheritsFrom (TH1::Class ())) {
      drawnHists.push_back ((TH1*)obj);
    }
    else if (obj->IsA()->InheritsFrom (TGAE::Class ())) {
      drawnGraphs.push_back ((TGAE*)obj);
    }
  }
}


void PhysicsAnalysis :: GetMinAndMax (double &min, double &max, const bool log) {
  for (TH1* h : drawnHists) {
    const double _max = log ? h->GetMaximum (0) : h->GetMaximum ();
    const double _min = log ? h->GetMinimum (0) : h->GetMinimum ();

    if (_max > max) max = _max;
    if (_min < min) min = _min;
  }
  for (TGAE* g : drawnGraphs) {
    const double* ys = g->GetY ();
    for (int n = 0; n < g->GetN (); n++) {
      if (log && ys[n] <= 0)
        continue;
      if (ys[n] > max) max = ys[n];
      if (ys[n] < min) min = ys[n];
    }
  }
  return;
}


void PhysicsAnalysis :: SetMinAndMax (double min, double max) {
  for (TH1* h : drawnHists) {
    h->GetYaxis ()->SetRangeUser (min, max);
  }
  for (TGAE* g : drawnGraphs) {
    g->GetYaxis ()->SetRangeUser (min, max);
  }
  gPad->Update ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Create TGraphAsymmErrors from a histogram
////////////////////////////////////////////////////////////////////////////////////////////////
TGAE* PhysicsAnalysis :: GetTGAE (TH1D* h) {
  return make_graph (h);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Sets errors in h based on variance matrix h2
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SetErrors (TH1D* h, TH2D* h2) {
  assert (h->GetNbinsX () == h2->GetNbinsX ());
  assert (h2->GetNbinsX () == h2->GetNbinsY ());

  for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
    assert (h2->GetBinContent (iX, iX) >= 0);
    h->SetBinError (iX, sqrt (h2->GetBinContent (iX, iX)));
  }
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Removes extra points stored in a TGAE to avoid showing them
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: TruncateTGAE (TGAE* g, const short iPtZ, const bool useTrkPt) {
  //TGAE* g = (*_g);
  if (!g) return;
  if ((useTrkPt && g->GetN () != maxNPtchBins) || (!useTrkPt && g->GetN () != maxNXhZBins)) {
    cout << "Warning: In PhysicsAnalysis :: TruncateTGAE: " << g->GetName () << " already had its bins trimmed!" << endl;
    return;
  }
  double x = 0, y = 0;
  for (int i = 0; i < g->GetN (); i++) {
    g->GetPoint (i, x, y);
    if (useTrkPt) {
      if (x < pTchBins[iPtZ][0] || pTchBins[iPtZ][nPtchBins[iPtZ]] < x) {
        g->RemovePoint (i);
        i--;
        continue;
      }
    }
    else {
      if (x < xhZBins[iPtZ][0] || xhZBins[iPtZ][nXhZBins[iPtZ]] < x) {
        g->RemovePoint (i);
        i--;
        continue;
      }
    }
  } // end loop over i

  double* bins = (useTrkPt ? pTchBins : xhZBins)[iPtZ];
  int nBins = (useTrkPt ? nPtchBins : nXhZBins)[iPtZ];
  assert (nBins <= g->GetN ());
  if (nBins == g->GetN ())
    return;

  for (int i = 0; i < nBins; i++) {
    g->GetPoint (i, x, y);
    int j = i;
    while (j < g->GetN () && x < bins[i+1]) {
      j++;
      g->GetPoint (j, x, y);
    }

    //if (j == i+1)
    //  continue;

    double sum = 0;
    double yherr = 0;
    double ylerr = 0;
    for (int k = i; k < j; k++) {
      g->GetPoint (k, x, y);
      const double width = g->GetErrorXhigh (k) + g->GetErrorXlow (k);
      sum += y * width;
      yherr += pow (g->GetErrorYhigh (k) * width, 2);
      ylerr += pow (g->GetErrorYlow (k) * width, 2);
    }
    assert (ylerr >= 0 && yherr >= 0);
    x = 0.5*(bins[i]+bins[i+1]);
    sum = sum / (fabs (bins[i+1]-bins[i]));
    yherr = sqrt (yherr) / fabs (bins[i+1]-bins[i]);
    ylerr = sqrt (ylerr) / fabs (bins[i+1]-bins[i]);

    g->SetPoint (i, x, sum); 
    g->SetPointEXhigh (i, 0.5*(bins[i+1]-bins[i]));
    g->SetPointEXlow (i, 0.5*(bins[i+1]-bins[i]));
    g->SetPointEYhigh (i, yherr);
    g->SetPointEYlow (i, ylerr);

    for (int k = i+1; k < j; k++)
      g->RemovePoint (i+1);
  }
  return;
}




/**
 * Combines bins in the provided histogram _h based on the x-axis bins we want for this pT^Z bin.
 */
void PhysicsAnalysis :: TruncateTH1D (TH1D** _h, const short iPtZ, const bool useTrkPt) {
  TH1D* h = (*_h);
  if (!h) return;
  if ((h->GetNbinsX () != maxNPtchBins && useTrkPt) || (h->GetNbinsX () != maxNXhZBins && !useTrkPt)) {
    //cout << "Warning: In PhysicsAnalysis :: TruncateTH1D: " << h->GetName () << " already had its bins trimmed!" << endl;
    return;
  }
  TH1D* hnew = new TH1D ("new", "", useTrkPt ? nPtchBins[iPtZ] : nXhZBins[iPtZ], useTrkPt ? pTchBins[iPtZ] : xhZBins[iPtZ]);
  for (int iX = 1; iX <= hnew->GetNbinsX (); iX++) {
    for (int iY = 1; iY <= h->GetNbinsX (); iY++) {
      if (hnew->GetBinLowEdge (iX) < h->GetBinCenter (iY) && h->GetBinCenter (iY) < hnew->GetBinLowEdge (iX) + hnew->GetBinWidth (iX)) {
        hnew->SetBinContent (iX, (hnew->GetBinContent (iX) * hnew->GetBinWidth (iX) + h->GetBinContent (iY) * h->GetBinWidth (iY)) / hnew->GetBinWidth (iX));
        //hnew->SetBinError (iX, sqrt (pow (hnew->GetBinError (iX) * hnew->GetBinWidth (iX), 2) + pow (h->GetBinError (iY) * h->GetBinWidth (iY), 2)) / hnew->GetBinWidth (iX));
        hnew->SetBinError (iX, h->GetBinError (iY));
      }
    } // end loop over iY
  } // end loop over iX
  TString newName = TString (h->GetName ());
  SaferDelete (_h);
  hnew->SetName (newName);
  (*_h) = hnew;
  return;
}



////////////////////////////////////////////////////////////////////////////////////////////////
// Clears histograms from memory (allows them to be overwritten).
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: ClearHists () {
  // Should clone these histograms
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        SaferDelete (&(h_trk_pt_ptz[iSpc][iPtZ][iCent]));
        SaferDelete (&(h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]));
        SaferDelete (&(h_trk_xhz_ptz[iSpc][iPtZ][iCent]));
        SaferDelete (&(h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]));
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          SaferDelete (&(h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]));
          SaferDelete (&(h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]));
          SaferDelete (&(h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]));
          SaferDelete (&(h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]));
          SaferDelete (&(h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]));
        } // end loop over iPhi
        for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
          SaferDelete (&(h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]));
        }
      } // end loop over iPtZ
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        SaferDelete (&(h_z_counts[iSpc][iPtZ][iCent]));
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

        SaferDelete (&(h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]));
        SaferDelete (&(h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent]));
        SaferDelete (&(h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]));
        SaferDelete (&(h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]));

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          SaferDelete (&(h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]));
          SaferDelete (&(h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]));
          SaferDelete (&(h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]));
          SaferDelete (&(h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]));
        } // end loop over iPhi
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins+1; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
        for (int iPtch = 0; iPtch < maxNPtchBins+2; iPtch++) {
          SaferDelete (&(h_trk_dphi[iSpc][iPtZ][iPtch][iCent]));
          SaferDelete (&(h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]));
        } // end loop over iPtch
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  histsReflected = false;
  backgroundSubtracted = false;
  sigToBkgCalculated = false;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        SaferDelete (&(h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]));
        SaferDelete (&(h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]));

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          SaferDelete (&(h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent]));
          SaferDelete (&(h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]));
        } // end loop over iPhi
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
  iaaCalculated = false;
  

  /*
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        SaferDelete (&(h_trk_pt_ptz_icp[iSpc][iPtZ][iCent]));
        SaferDelete (&(h_trk_xhz_ptz_icp[iSpc][iPtZ][iCent]));

        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          SaferDelete (&(h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]));
          SaferDelete (&(h_trk_pt_dphi_icp[iSpc][iPtZ][iPhi][iCent]));
        } // end loop over iPhi
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
  icpCalculated = false;
  */

  histsLoaded = false;
  histsScaled = false;

  SaferDelete (&histFile);
  return;    

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CreateHists () {
  h_fcal_et = new TH1D (Form ("h_fcal_et_%s", name.c_str ()), "", numSuperFineCentBins-1, superFineCentBins); 
  h_fcal_et->Sumw2 ();
  h_fcal_et_reweighted = new TH1D (Form ("h_fcal_et_reweighted_%s", name.c_str ()), "", numSuperFineCentBins-1, superFineCentBins);
  h_fcal_et_reweighted->Sumw2 ();
  for (short iMBTrig = 0; iMBTrig < 3; iMBTrig++) {
    h_centrality[iMBTrig] = new TH1D (Form ("h_centrality_trig%i_%s", iMBTrig, name.c_str ()), "", 80, 0, 80);
    h_centrality[iMBTrig]->Sumw2 ();
    h_centrality_reweighted[iMBTrig] = new TH1D (Form ("h_centrality_reweighted_trig%i_%s", iMBTrig, name.c_str ()), "", 80, 0, 80);
    h_centrality_reweighted[iMBTrig]->Sumw2 ();
  } // end loop over iMBTrig

  for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
    h_q2[iFineCent]               = new TH1D (Form ("h_q2_iCent%i_%s", iFineCent, name.c_str ()), "", 20, 0, 0.3);
    h_q2[iFineCent]->Sumw2 ();
    h_q2_reweighted[iFineCent]    = new TH1D (Form ("h_q2_reweighted_iCent%i_%s", iFineCent, name.c_str ()), "", 20, 0, 0.3);
    h_q2_reweighted[iFineCent]->Sumw2 ();
    h_psi2[iFineCent]             = new TH1D (Form ("h_psi2_iCent%i_%s", iFineCent, name.c_str ()), "", 8, -pi/2, pi/2);
    h_psi2[iFineCent]->Sumw2 ();
    h_psi2_reweighted[iFineCent]  = new TH1D (Form ("h_psi2_reweighted_iCent%i_%s", iFineCent, name.c_str ()), "", 8, -pi/2, pi/2);
    h_psi2_reweighted[iFineCent]->Sumw2 ();
  } // end loop over iFineCent
  h_PbPb_vz = new TH1D (Form ("h_PbPb_vz_%s", name.c_str ()), "", 50, -200, 200);
  h_PbPb_vz_reweighted = new TH1D (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_vz = new TH1D (Form ("h_pp_vz_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_vz_reweighted = new TH1D (Form ("h_pp_vz_reweighted_%s", name.c_str ()), "", 50, -200, 200);
  h_pp_nch = new TH1D (Form ("h_pp_nch_%s", name.c_str ()), "", 80, -0.5, 160.5);
  h_pp_nch_reweighted = new TH1D (Form ("h_pp_nch_reweighted_%s", name.c_str ()), "", 80, -0.5, 160.5);

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_pt_dphi_raw_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", maxNPtchBins, allPtchBins); // old name: h_z_trk_raw_pt
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_pt_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", maxNPtchBins, allPtchBins); // old name: h_z_trk_pt
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent] = new TH2D (Form ("h2_trk_pt_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", maxNPtchBins, allPtchBins, maxNPtchBins, allPtchBins);
          h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_xhz_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", maxNXhZBins, allXhZBins); // old name: h_z_trk_xzh
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
          h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent] = new TH2D (Form ("h2_trk_xhz_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "", maxNXhZBins, allXhZBins, maxNXhZBins, allXhZBins);
          h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]->Sumw2 ();
        } // end loop over iPhi
        for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent] = new TH1D (Form ("h_trk_dphi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()), "", GetNdPhiBins (0, 0), -pi/2, 3*pi/2); // old name: h_z_trk_phi
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->Sumw2 ();
          h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent] = new TH2D (Form ("h2_trk_dphi_cov_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()), "", GetNdPhiBins (0, 0), -pi/2, 3*pi/2, GetNdPhiBins (0, 0), -pi/2, 3*pi/2);
          h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->Sumw2 ();
        } // end loop over iPtch
        h_trk_pt_ptz[iSpc][iPtZ][iCent] = new TH1D (Form ("h_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", maxNPtchBins, allPtchBins); // old name: h_z_trk_zpt
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->Sumw2 ();
        h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent] = new TH2D (Form ("h2_trk_pt_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", maxNPtchBins, allPtchBins, maxNPtchBins, allPtchBins);
        h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->Sumw2 ();
        h_trk_xhz_ptz[iSpc][iPtZ][iCent] = new TH1D (Form ("h_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", maxNXhZBins, allXhZBins); // old name: h_z_trk_zxzh
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Sumw2 ();
        h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent] = new TH2D (Form ("h2_trk_xhz_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", maxNXhZBins, allXhZBins, maxNXhZBins, allXhZBins);
        h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->Sumw2 ();

        h_z_counts[iSpc][iPtZ][iCent] = new TH1D (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", 3, 0, 3);
        h_z_counts[iSpc][iPtZ][iCent]->Sumw2 ();
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

  histsLoaded = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms as clones from another analysis, where appropriate
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CopyAnalysis (PhysicsAnalysis* a, const bool copySigs, const bool copyIAAs) {
  if (name == "")
    cout << "Warning in PhysicsAnalysis :: CopyAnalysis: name of analysis not set!" << endl;

  ClearHists ();

  //// Should clone these histograms
  //h_fcal_et               = (TH1D*) a->h_fcal_et->Clone (Form ("h_fcal_et_%s", name.c_str ()));
  //h_fcal_et_reweighted    = (TH1D*) a->h_fcal_et_reweighted->Clone (Form ("h_fcal_et_reweighted_%s", name.c_str ()));
  //for (short iMBTrig = 0; iMBTrig < 3; iMBTrig++) {
  //  h_centrality[iMBTrig]             = (TH1D*) a->h_centrality[iMBTrig]->Clone (Form ("h_centrality_trig%i_%s", iMBTrig, name.c_str ()));
  //  h_centrality_reweighted[iMBTrig]  = (TH1D*) a->h_centrality[iMBTrig]->Clone (Form ("h_centrality_reweighted_trig%i_%s", iMBTrig, name.c_str ()));
  //} // end loop over iMBTrig

  //for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
  //  h_q2[iFineCent]               = (TH1D*) a->h_q2[iFineCent]->Clone (Form ("h_q2_iCent%i_%s", iFineCent, name.c_str ()));
  //  h_q2_reweighted[iFineCent]    = (TH1D*) a->h_q2_reweighted[iFineCent]->Clone (Form ("h_q2_reweighted_iCent%i_%s", iFineCent, name.c_str ()));
  //  h_psi2[iFineCent]             = (TH1D*) a->h_psi2[iFineCent]->Clone (Form ("h_psi2_iCent%i_%s", iFineCent, name.c_str ()));
  //  h_psi2_reweighted[iFineCent]  = (TH1D*) a->h_psi2_reweighted[iFineCent]->Clone (Form ("h_psi2_reweighted_iCent%i_%s", iFineCent, name.c_str ()));
  //} // end loop over iFineCent
  //h_PbPb_vz             = (TH1D*) a->h_PbPb_vz->Clone (Form ("h_PbPb_vz_%s", name.c_str ()));
  //h_PbPb_vz_reweighted  = (TH1D*) a->h_PbPb_vz_reweighted->Clone (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()));
  //h_pp_vz               = (TH1D*) a->h_pp_vz->Clone (Form ("h_pp_vz_%s", name.c_str ()));
  //h_pp_vz_reweighted    = (TH1D*) a->h_pp_vz_reweighted->Clone (Form ("h_pp_vz_reweighted_%s", name.c_str ()));
  //h_pp_nch              = (TH1D*) a->h_pp_nch->Clone (Form ("h_pp_nch_%s", name.c_str ()));
  //h_pp_nch_reweighted   = (TH1D*) a->h_pp_nch_reweighted->Clone (Form ("h_pp_nch_reweighted_%s", name.c_str ()));

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]       = (TH1D*) a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]  = (TH2D*) a->h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->Clone (Form ("h2_trk_pt_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]      = (TH1D*) a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh
        h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent] = (TH2D*) a->h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->Clone (Form ("h2_trk_xhz_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh

        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]    = (TH1D*) a->h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_raw_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_raw_pt
        //  h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]        = (TH1D*) a->h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt
        //  h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]   = (TH2D*) a->h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h2_trk_pt_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt
        //  h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]       = (TH1D*) a->h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh
        //  h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]  = (TH2D*) a->h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h2_trk_xhz_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh
        //} // end loop over iPhi
      } // end loop over iPtZ
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_counts[iSpc][iPtZ][iCent] = (TH1D*) a->h_z_counts[iSpc][iPtZ][iCent]->Clone (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent


  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
          if (a->h_trk_dphi[iSpc][iPtZ][iPtch][iCent])
            h_trk_dphi[iSpc][iPtZ][iPtch][iCent] = (TH1D*) a->h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h_trk_dphi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
          if (a->h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent])
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent] = (TH2D*) a->h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h2_trk_dphi_cov_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
        } // end loop over iPtch
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
  MakeSummaryCorrelations ();


  if (copySigs && a->backgroundSubtracted) {

    {
      const short iSpc = 2;
      const char* spc = "comb";
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        for (short iCent = 0; iCent <= numCentBins; iCent++) {
          for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
            if (a->h_trk_dphi[iSpc][iPtZ][iPtch][iCent])
              h_trk_dphi[iSpc][iPtZ][iPtch][iCent] = (TH1D*) a->h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h_trk_dphi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
          } // end loop over iPtch
          if (a->h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent])
            h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent] = (TH1D*) a->h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent]->Clone (Form ("h_trk_dphi_gt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_phi
          if (a->h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent])
            h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent] = (TH1D*) a->h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent]->Clone (Form ("h_trk_dphi_lt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_phi
        } // end loop over iCent
      } // end loop over iPtZ
    } // end loop over iSpc

    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        {
          const short iCent = numCentBins;
          for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
            if (a->h_trk_dphi[iSpc][iPtZ][iPtch][iCent])
              h_trk_dphi[iSpc][iPtZ][iPtch][iCent] = (TH1D*) a->h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h_trk_dphi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
            if (a->h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent])
              h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent] = (TH1D*) a->h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h_trk_dphi_sub_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
          } // end loop over iPtch
          if (a->h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent])
            h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent] = (TH1D*) a->h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent]->Clone (Form ("h_trk_dphi_gt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_phi
          if (a->h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent])
            h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent] = (TH1D*) a->h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent]->Clone (Form ("h_trk_dphi_lt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_phi
        } // end iCent = numCentBins scope

        for (short iCent = 0; iCent <= numCentBins; iCent++) {
          for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
            if (a->h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent])
              h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent] = (TH1D*) a->h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h_trk_dphi_sub_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
          } // end loop over iPtch
          if (a->h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent])
            h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent] = (TH1D*) a->h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent]->Clone (Form ("h_trk_dphi_sub_gt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_phi
          if (a->h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent])
            h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent] = (TH1D*) a->h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent]->Clone (Form ("h_trk_dphi_sub_lt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_phi
        } // end loop over iCent
      } // end loop over iSpc
    } // end loop over iPtZ

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

          if (a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent])
            h_trk_pt_ptz_sub[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt_sub
          if (a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent])
            h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh_sub

          //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          //  if (a->h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent])
          //    h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt_sub
          //  if (a->h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent])
          //    h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh_sub
          //} // end loop over iPhi
        } // end loop over iPtZ
      } // end loop over iCent
    } // end loop over iSpc

    backgroundSubtracted = true;
    bkg = a->bkg;
  }

  if (copySigs && copyIAAs && a->iaaCalculated ) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {

          if (a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent])
            h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt_iaa
          if (a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent])
            h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent] = (TH1D*) a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh_iaa

          //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          //  if (a->h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent])
          //    h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt_iaa
          //  if (a->h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent])
          //    h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent] = (TH1D*) a->h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh_iaa
          //} // end loop over iPhi
        } // end loop over iCent
      } // end loop over iPtZ
    } // end loop over iSpc
    iaaCalculated = true;
  }

  histsLoaded = true;
  histsScaled = true;
  return;    
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  //if (histsLoaded)
  ClearHists ();

  TDirectory* _gDirectory = gDirectory;
  histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "read");
      
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        for (short iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent] = (TH1D*) histFile->Get (Form ("h_trk_dphi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ())); // old name: h_z_trk_phi
          h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent] = (TH2D*) histFile->Get (Form ("h2_trk_dphi_cov_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()));
        } // end loop over iPtch

        h_trk_pt_ptz[iSpc][iPtZ][iCent]       = (TH1D*) histFile->Get (Form ("h_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]  = (TH2D*) histFile->Get (Form ("h2_trk_pt_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]      = (TH1D*) histFile->Get (Form ("h_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh
        h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent] = (TH2D*) histFile->Get (Form ("h2_trk_xhz_ptz_cov_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]    = (TH1D*) histFile->Get (Form ("h_trk_pt_dphi_raw_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_raw_pt
          h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]        = (TH1D*) histFile->Get (Form ("h_trk_pt_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_pt
          h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]   = (TH2D*) histFile->Get (Form ("h2_trk_pt_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]       = (TH1D*) histFile->Get (Form ("h_trk_xhz_dphi_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())); // old name: h_z_trk_xzh
          h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]  = (TH2D*) histFile->Get (Form ("h2_trk_xhz_dphi_cov_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        } // end loop over iPhi
      } // end loop over iPtZ
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_counts[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

  histsLoaded = true;

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        h_z_counts[2][iPtZ][iCent]->Add (h_z_counts[iSpc][iPtZ][iCent]);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]) h_trk_pt_dphi_raw[2][iPtZ][iPhi][iCent]->Add (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over iPhi
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  double newbins[numSuperFineCentBins];
  for (int i = 0; i < numSuperFineCentBins; i++)
    newbins[i] = superFineCentBins[i]/1000.;

  h_fcal_et             = new TH1D (Form ("h_fcal_et_fixedBins_%s", name.c_str ()), "", sizeof (newbins)/sizeof (newbins[0]) - 1, newbins);
  h_fcal_et_reweighted  = new TH1D (Form ("h_fcal_et_reweighted_fixedBins_%s", name.c_str ()), "", sizeof (newbins)/sizeof (newbins[0]) - 1, newbins);

  TH1D* temp = (TH1D*) histFile->Get (Form ("h_fcal_et_%s", name.c_str ()));
  assert (h_fcal_et->GetNbinsX () == temp->GetNbinsX ());
  for (int iX = 1; iX <= temp->GetNbinsX (); iX++)
    h_fcal_et->SetBinContent (iX, temp->GetBinContent (iX));

  temp = (TH1D*) histFile->Get (Form ("h_fcal_et_reweighted_%s", name.c_str ()));
  assert (h_fcal_et_reweighted->GetNbinsX () == temp->GetNbinsX ());
  for (int iX = 1; iX <= temp->GetNbinsX (); iX++)
    h_fcal_et_reweighted->SetBinContent (iX, temp->GetBinContent (iX));

  for (short iMBTrig = 0; iMBTrig < 3; iMBTrig++) {
    h_centrality[iMBTrig]             = (TH1D*) histFile->Get (Form ("h_centrality_trig%i_%s", iMBTrig, name.c_str ()));
    h_centrality_reweighted[iMBTrig]  = (TH1D*) histFile->Get (Form ("h_centrality_reweighted_trig%i_%s", iMBTrig, name.c_str ()));
  } // end loop over iMBTrig

  for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
    h_q2[iFineCent]               = (TH1D*) histFile->Get (Form ("h_q2_iCent%i_%s", iFineCent, name.c_str ()));
    h_q2_reweighted[iFineCent]    = (TH1D*) histFile->Get (Form ("h_q2_reweighted_iCent%i_%s", iFineCent, name.c_str ()));
    h_psi2[iFineCent]             = (TH1D*) histFile->Get (Form ("h_psi2_iCent%i_%s", iFineCent, name.c_str ()));
    h_psi2_reweighted[iFineCent]  = (TH1D*) histFile->Get (Form ("h_psi2_reweighted_iCent%i_%s", iFineCent, name.c_str ()));
  } // end loop over iFineCent
  h_PbPb_vz             = (TH1D*) histFile->Get (Form ("h_PbPb_vz_%s", name.c_str ()));
  h_PbPb_vz_reweighted  = (TH1D*) histFile->Get (Form ("h_PbPb_vz_reweighted_%s", name.c_str ()));
  h_pp_vz               = (TH1D*) histFile->Get (Form ("h_pp_vz_%s", name.c_str ()));
  h_pp_vz_reweighted    = (TH1D*) histFile->Get (Form ("h_pp_vz_reweighted_%s", name.c_str ()));
  h_pp_nch               = (TH1D*) histFile->Get (Form ("h_pp_nch_%s", name.c_str ()));
  h_pp_nch_reweighted    = (TH1D*) histFile->Get (Form ("h_pp_nch_reweighted_%s", name.c_str ()));
  

  if (_finishHists) {
    ScaleHists ();
    SetVariances ();
    MakeSummaryCorrelations ();
  }

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SaveHists (const char* histFileName) {
  if (!histsLoaded)
    return;

  TDirectory* _gDirectory = gDirectory;
  if (!histFile) {
    histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "recreate");
    histFile->cd ();
  }

  cout << "Saving histograms to " << Form ("%s/%s", rootPath.Data (), histFileName) << endl;

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {

        for (short iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
          SafeWrite (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]);
          SafeWrite (h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]);
        } // end loop over iPtch

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          SafeWrite (h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]);
          SafeWrite (h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over iPhi

        SafeWrite (h_trk_pt_ptz[iSpc][iPtZ][iCent]);
        SafeWrite (h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]);
        SafeWrite (h_trk_xhz_ptz[iSpc][iPtZ][iCent]);
        SafeWrite (h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]);

        SafeWrite (h_z_counts[iSpc][iPtZ][iCent]);
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

  SafeWrite (h_fcal_et);
  SafeWrite (h_fcal_et_reweighted);
  for (short iMBTrig = 0; iMBTrig < 3; iMBTrig++) {
    SafeWrite (h_centrality[iMBTrig]);
    SafeWrite (h_centrality_reweighted[iMBTrig]);
  } // end loop over iMBTrig

  for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
    SafeWrite (h_q2[iFineCent]);
    SafeWrite (h_q2_reweighted[iFineCent]);
    SafeWrite (h_psi2[iFineCent]);
    SafeWrite (h_psi2_reweighted[iFineCent]);
  } // end loop over iFineCent

  SafeWrite (h_PbPb_vz);
  SafeWrite (h_PbPb_vz_reweighted);
  SafeWrite (h_pp_vz);
  SafeWrite (h_pp_vz_reweighted);
  SafeWrite (h_pp_nch);
  SafeWrite (h_pp_nch_reweighted);
  
  histFile->Close ();
  histFile = nullptr;
  histsLoaded = false;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SaveResults (const char* saveFileName) {
  if (!histsLoaded)
    return;

  TDirectory* _gDirectory = gDirectory;
  TFile* saveFile = new TFile (Form ("%s/%s", rootPath.Data (), saveFileName), "recreate");
  saveFile->cd ();

  cout << "Saving results to " << Form ("%s/%s", rootPath.Data (), saveFileName) << endl;

  for (short iCent = 0; iCent < numCentBins+1; iCent++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      SafeWrite (h_trk_dphi[2][iPtZ][maxNPtchBins][iCent],      Form ("h_trk_dphi_gt4_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      SafeWrite (h_trk_dphi[2][iPtZ][maxNPtchBins+1][iCent],      Form ("h_trk_dphi_lt4_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      SafeWrite (h_trk_dphi_sub[2][iPtZ][maxNPtchBins][iCent],  Form ("h_trk_dphi_sub_gt4_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      SafeWrite (h_trk_dphi_sub[2][iPtZ][maxNPtchBins+1][iCent],  Form ("h_trk_dphi_sub_lt4_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
    } // end loop over iPtZ
  }

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        for (short iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          SafeWrite (h_trk_dphi[iSpc][iPtZ][iPtch][iCent],      Form ("h_trk_dphi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()));
          SafeWrite (h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent],  Form ("h_trk_dphi_sub_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()));
        } // end loop over iPtZ

        SafeWrite (h_trk_pt_ptz[iSpc][iPtZ][iCent],                Form ("h_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent],            Form ("h_trk_pt_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent],     Form ("h_trk_pt_ptz_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz[iSpc][iPtZ][iCent],               Form ("h_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent],           Form ("h_trk_xhz_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent],    Form ("h_trk_xhz_ptz_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));

        SafeWrite (h_z_counts[iSpc][iPtZ][iCent],                 Form ("h_z_counts_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      } // end loop over iPtZ
    } // end loop over iSpc
   } // end loop over iCent

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        SafeWrite (h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent],            Form ("h_trk_pt_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        SafeWrite (h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent],           Form ("h_trk_xhz_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));

        if (iSpc == 2 && h_trk_pt_ptz_iaa_2015proj[iPtZ][iCent]) SafeWrite (h_trk_pt_ptz_iaa_2015proj[iPtZ][iCent]);
        if (iSpc == 2 && h_trk_xhz_ptz_iaa_2015proj[iPtZ][iCent]) SafeWrite (h_trk_xhz_ptz_iaa_2015proj[iPtZ][iCent]);
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      SafeWrite (g_trk_avg_pt_ptz[iSpc][iCent],   Form ("g_trk_avg_pt_ptz_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      SafeWrite (g_trk_avg_xhz_ptz[iSpc][iCent],  Form ("g_trk_avg_xhz_ptz_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      SafeWrite (g_avg_ntrk_ptz[iSpc][iCent],   Form ("g_avg_ntrk_ptz_%s_iCent%i_%s", spc, iCent, name.c_str ()));
    } // end loop over iSpc
  } // end loop over iCent
  
  saveFile->Close ();
  saveFile = nullptr;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Scale histograms for plotting, calculating signals, etc.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: ScaleHists () {
  if (histsScaled || !histsLoaded)
    return;

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        TH1D* h = h_z_counts[iSpc][iPtZ][iCent];
        const double counts = h->GetBinContent (1);
        const double sumWeights = h->GetBinContent (2);
        const double sumWeightsSq = h->GetBinContent (3);
        if (counts <= 0)  continue;

        TH2D* h2 = nullptr;

        // finalize covariance calculation by normalizing to bin widths and subtracting off product of means
        // then use diagonals of the covariance matrix to determine statistical uncertainties
        // see internal documentation or https://en.wikipedia.org/wiki/Sample_mean_and_covariance for more details
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h = h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent];
          h2 = h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent];
          assert (h2->GetNbinsX () == h->GetNbinsX ());
          h->Scale (1/ (sumWeights * (isBkg && !doPPMBMixing && iCent == 0 ? pi/8. : (phiHighBins[iPhi]-phiLowBins[iPhi]))), "width");
          h2->Scale (1/ pow (isBkg && !doPPMBMixing && iCent == 0 ? pi/8. : (phiHighBins[iPhi]-phiLowBins[iPhi]), 2), "width");
          for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
              h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (sumWeights)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
          h2->Scale (sumWeights / (counts * (sumWeights*sumWeights - sumWeightsSq)));

          h = h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];
          h2 = h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent];
          assert (h2->GetNbinsX () == h->GetNbinsX ());
          h->Scale (1/ (sumWeights * (isBkg && !doPPMBMixing && iCent == 0 ? pi/8. : (phiHighBins[iPhi]-phiLowBins[iPhi]))), "width");
          h2->Scale (1/ pow (isBkg && !doPPMBMixing && iCent == 0 ? pi/8. : (phiHighBins[iPhi]-phiLowBins[iPhi]), 2), "width");
          for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
              h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (sumWeights)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
          h2->Scale (sumWeights / (counts * (sumWeights*sumWeights - sumWeightsSq)));
        } // end loop over iPhi


        h = h_trk_pt_ptz[iSpc][iPtZ][iCent];
        h2 = h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent];
        assert (h2->GetNbinsX () == h->GetNbinsX ());
        h->Scale (1/ (sumWeights * (isBkg && !doPPMBMixing && iCent == 0 ? pi/8. : pi/4.)), "width");
        h2->Scale (1/ pow (isBkg && !doPPMBMixing && iCent == 0 ? pi/8. : pi/4., 2), "width");
        for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
          for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
            h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (sumWeights)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
        h2->Scale (sumWeights / (counts * (sumWeights*sumWeights - sumWeightsSq)));

        h = h_trk_xhz_ptz[iSpc][iPtZ][iCent];
        h2 = h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent];
        assert (h2->GetNbinsX () == h->GetNbinsX ());
        h->Scale (1/ (sumWeights * (isBkg && !doPPMBMixing && iCent == 0 ? pi/8. : pi/4.)), "width");
        h2->Scale (1/ pow (isBkg && !doPPMBMixing && iCent == 0 ? pi/8. : pi/4., 2), "width");
        for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
          for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
            h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (sumWeights)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
        h2->Scale (sumWeights / (counts * (sumWeights*sumWeights - sumWeightsSq)));

        for (short iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
          h = h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
          h2 = h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent];
          assert (h2->GetNbinsX () == h->GetNbinsX ());
          h->Scale (1/sumWeights, "width");
          h2->Scale (1, "width");
          for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
              h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (sumWeights)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
          h2->Scale (sumWeights / (counts * (sumWeights*sumWeights - sumWeightsSq)));
        } // end loop over iPtch
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  histsScaled = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Scale histograms for plotting, calculating signals, etc.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SetVariances () {
  assert (histsLoaded);

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        const float counts = h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1);
        if (counts <= 0)  continue;

        TH1D* h = nullptr;
        TH2D* h2 = nullptr;

        // finalize covariance calculation by normalizing to bin widths and subtracting off product of means
        // then use diagonals of the covariance matrix to determine statistical uncertainties
        // see internal documentation or https://en.wikipedia.org/wiki/Sample_mean_and_covariance for more details
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          h = h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent];
          h2 = h2_trk_pt_dphi_cov[iSpc][iPtZ][iPhi][iCent];
          PhysicsAnalysis :: SetErrors (h, h2);

          h = h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];
          h2 = h2_trk_xhz_dphi_cov[iSpc][iPtZ][iPhi][iCent];
          PhysicsAnalysis :: SetErrors (h, h2);
        } // end loop over iPhi


        h = h_trk_pt_ptz[iSpc][iPtZ][iCent];
        h2 = h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent];
        PhysicsAnalysis :: SetErrors (h, h2);

        h = h_trk_xhz_ptz[iSpc][iPtZ][iCent];
        h2 = h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent];
        PhysicsAnalysis :: SetErrors (h, h2);

        for (short iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
          h = h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
          h2 = h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent];
          PhysicsAnalysis :: SetErrors (h, h2);
        } // end loop over iPtch
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates pT^ch- and centrality-integrated dphi correlations
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: MakeSummaryCorrelations () {

  if (histsReflected) {
    cout << "Correlations already reflected for " << name.c_str () << endl;
    return;
  }

  /**** Create empty histograms for azimuthal plots integrated over pTch ****/
  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

    double* lt4_abbrev_dphi_bins = nullptr;
    short num_lt4_abbrev_dphi_bins = 0;
    if (iPtZ >= 4) {
      lt4_abbrev_dphi_bins = new double[8];
      lt4_abbrev_dphi_bins[0] = 0;
      lt4_abbrev_dphi_bins[1] = 5.*pi/20.;
      lt4_abbrev_dphi_bins[2] = 9.*pi/20.;
      lt4_abbrev_dphi_bins[3] = 13.*pi/20.;
      lt4_abbrev_dphi_bins[4] = 15.*pi/20.;
      lt4_abbrev_dphi_bins[5] = 17.*pi/20.;
      lt4_abbrev_dphi_bins[6] = 19.*pi/20.;
      lt4_abbrev_dphi_bins[7] = pi;
      num_lt4_abbrev_dphi_bins = 7;
    }
    else if (iPtZ == 3) {
      lt4_abbrev_dphi_bins = new double[7];
      lt4_abbrev_dphi_bins[0] = 0;
      lt4_abbrev_dphi_bins[1] = 4.*pi/20.;
      lt4_abbrev_dphi_bins[2] = 8.*pi/20.;
      lt4_abbrev_dphi_bins[3] = 12.*pi/20.;
      lt4_abbrev_dphi_bins[4] = 15.*pi/20.;
      lt4_abbrev_dphi_bins[5] = 18.*pi/20.;
      lt4_abbrev_dphi_bins[6] = pi;
      num_lt4_abbrev_dphi_bins = 6;
    }
    else {
      lt4_abbrev_dphi_bins = new double[4];
      lt4_abbrev_dphi_bins[0] = 0;
      lt4_abbrev_dphi_bins[1] = 8.*pi/20.;
      lt4_abbrev_dphi_bins[2] = 15.*pi/20.;
      lt4_abbrev_dphi_bins[3] = pi;
      num_lt4_abbrev_dphi_bins = 3;
    }

    double* gt4_abbrev_dphi_bins = nullptr;
    short num_gt4_abbrev_dphi_bins = 0;
    if (iPtZ >= 3) {
      gt4_abbrev_dphi_bins = new double[21];
      gt4_abbrev_dphi_bins[0] = 0;
      gt4_abbrev_dphi_bins[1] = 1.*pi/20.;
      gt4_abbrev_dphi_bins[2] = 2.*pi/20.;
      gt4_abbrev_dphi_bins[3] = 3.*pi/20.;
      gt4_abbrev_dphi_bins[4] = 4.*pi/20.;
      gt4_abbrev_dphi_bins[5] = 5.*pi/20.;
      gt4_abbrev_dphi_bins[6] = 6.*pi/20.;
      gt4_abbrev_dphi_bins[7] = 7.*pi/20.;
      gt4_abbrev_dphi_bins[8] = 8.*pi/20.;
      gt4_abbrev_dphi_bins[9] = 9.*pi/20.;
      gt4_abbrev_dphi_bins[10] = 10.*pi/20.;
      gt4_abbrev_dphi_bins[11] = 11.*pi/20.;
      gt4_abbrev_dphi_bins[12] = 12.*pi/20.;
      gt4_abbrev_dphi_bins[13] = 13.*pi/20.;
      gt4_abbrev_dphi_bins[14] = 14.*pi/20.;
      gt4_abbrev_dphi_bins[15] = 15.*pi/20.;
      gt4_abbrev_dphi_bins[16] = 16.*pi/20.;
      gt4_abbrev_dphi_bins[17] = 17.*pi/20.;
      gt4_abbrev_dphi_bins[18] = 18.*pi/20.;
      gt4_abbrev_dphi_bins[19] = 19.*pi/20.;
      gt4_abbrev_dphi_bins[20] = pi;
      num_gt4_abbrev_dphi_bins = 20;
    }
    else {
      gt4_abbrev_dphi_bins = new double[14];
      gt4_abbrev_dphi_bins[0] = 0;
      gt4_abbrev_dphi_bins[1] = 2.*pi/20.;
      gt4_abbrev_dphi_bins[2] = 4.*pi/20.;
      gt4_abbrev_dphi_bins[3] = 6.*pi/20.;
      gt4_abbrev_dphi_bins[4] = 8.*pi/20.;
      gt4_abbrev_dphi_bins[5] = 10.*pi/20.;
      gt4_abbrev_dphi_bins[6] = 12.*pi/20.;
      gt4_abbrev_dphi_bins[7] = 14.*pi/20.;
      gt4_abbrev_dphi_bins[8] = 15.*pi/20.;
      gt4_abbrev_dphi_bins[9] = 16.*pi/20.;
      gt4_abbrev_dphi_bins[10] = 17.*pi/20.;
      gt4_abbrev_dphi_bins[11] = 18.*pi/20.;
      gt4_abbrev_dphi_bins[12] = 19.*pi/20.;
      gt4_abbrev_dphi_bins[13] = pi;
      num_gt4_abbrev_dphi_bins = 13;
    }

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iCent = 0; iCent <= numCentBins; iCent++) {
        SaferDelete (&h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent]);
        h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent] = new TH1D (Form ("h_trk_dphi_gt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", num_gt4_abbrev_dphi_bins, gt4_abbrev_dphi_bins);
        h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent]->Sumw2 ();
        SaferDelete (&h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent]);
        h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent] = new TH1D (Form ("h_trk_dphi_lt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", num_lt4_abbrev_dphi_bins, lt4_abbrev_dphi_bins);
        h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent]->Sumw2 ();

        if (hasBkg) {
          SaferDelete (&h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent]);
          h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent] = new TH1D (Form ("h_trk_dphi_sub_gt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", num_gt4_abbrev_dphi_bins, gt4_abbrev_dphi_bins);
          h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent]->Sumw2 ();
          SaferDelete (&h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent]);
          h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent] = new TH1D (Form ("h_trk_dphi_sub_lt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "", num_lt4_abbrev_dphi_bins, lt4_abbrev_dphi_bins);
          h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent]->Sumw2 ();
        }
      } // end loop over iCent
    } // end loop over iSpc
    delete[] lt4_abbrev_dphi_bins;
    delete[] gt4_abbrev_dphi_bins;
  } // end loop over iPtZ


  /**** Now create 2-4 GeV plots and 4+ GeV plots for electrons and for muons (not for combined plots yet) ****/
  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (int iPtch = 2; iPtch < maxNPtchBins; iPtch++) {
          TH1D* h = (TH1D*) h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
          TH2D* h2 = (TH2D*) h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent];
          TH1D* hsum = h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent];

          for (short iX = 1; iX <= hsum->GetNbinsX (); iX++) {
            std::vector<short> xbins = {};
            for (short iX1 = 1; iX1 <= h->GetNbinsX (); iX1++) {
              const double x = acos (cos (h->GetBinCenter (iX1)));
              if (hsum->GetBinLowEdge (iX) < x && x < hsum->GetBinLowEdge (iX) + hsum->GetBinWidth (iX))
                xbins.push_back (iX1);
            }
            double content = hsum->GetBinContent (iX) * hsum->GetBinWidth (iX);
            double error = pow (hsum->GetBinError (iX) * hsum->GetBinWidth (iX), 2);
            for (short iX1 : xbins) {
              content += h->GetBinContent (iX1)*h->GetBinWidth (iX1);
              for (short iX2 : xbins) {
                error += (h2->GetBinContent (iX1, iX2)) * (h2->GetXaxis ()->GetBinWidth (iX1)) * (h2->GetYaxis ()->GetBinWidth (iX2));
              }
            }
            hsum->SetBinContent (iX, content / hsum->GetBinWidth (iX));
            hsum->SetBinError (iX, sqrt (fabs (error)) / hsum->GetBinWidth (iX));
          }
        } // end loop over iPtch
        for (int iPtch = 1; iPtch < 2; iPtch++) {
          TH1D* h = (TH1D*) h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
          TH2D* h2 = (TH2D*) h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent];
          TH1D* hsum = h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent];

          for (short iX = 1; iX <= hsum->GetNbinsX (); iX++) {
            std::vector<short> xbins = {};
            for (short iX1 = 1; iX1 <= h->GetNbinsX (); iX1++) {
              const double x = acos (cos (h->GetBinCenter (iX1)));
              if (hsum->GetBinLowEdge (iX) < x && x < hsum->GetBinLowEdge (iX) + hsum->GetBinWidth (iX))
                xbins.push_back (iX1);
            }
            double content = hsum->GetBinContent (iX) * hsum->GetBinWidth (iX);
            double error = pow (hsum->GetBinError (iX) * hsum->GetBinWidth (iX), 2);
            for (short iX1 : xbins) {
              content += h->GetBinContent (iX1)*h->GetBinWidth (iX1);
              for (short iX2 : xbins) {
                error += (h2->GetBinContent (iX1, iX2)) * (h2->GetXaxis ()->GetBinWidth (iX1)) * (h2->GetYaxis ()->GetBinWidth (iX2));
              }
            }
            hsum->SetBinContent (iX, content / hsum->GetBinWidth (iX));
            hsum->SetBinError (iX, sqrt (fabs (error)) / hsum->GetBinWidth (iX));
          }
        } // end loop over iPtch
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  histsReflected = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Subtracts mixed event background from track yields
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: SubtractBackground (PhysicsAnalysis* _bkg, const bool addUnc) {
  if (isBkg) {
    cout << "Cannot subtract background on background analysis" << endl;
    return;
  }

  if (backgroundSubtracted) {
    cout << "Background already subtracted for " << name << endl;
    return;
  }

  if (!_bkg) {
    cout << "No background provided! Will not do subtraction." << endl;
    return;
  }

  bkg = _bkg;

  bkg->MakeSummaryCorrelations ();

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

        TH1D* sub = nullptr, *h = nullptr;

        //******** Do subtraction of integrated dPhi plots ********//
        sub = bkg->h_trk_pt_ptz[iSpc][iPtZ][iCent];
        SaferDelete (&(h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]));
        h = (TH1D*) h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_pt_ptz_sub[iSpc][iPtZ][iCent] = h;
        if (iCent != 0 || subtractPP) {
          if (!isBkg && !addUnc)      AddNoErrors (h, sub, -1);
          else if (!isBkg && addUnc)  h->Add (sub, -1);
          else                        h->Reset ();
        }

        sub = bkg->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
        SaferDelete (&(h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]));
        h = (TH1D*) h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent] = h;
        if (iCent != 0 || subtractPP) {
          if (!isBkg && !addUnc)      AddNoErrors (h, sub, -1);
          else if (!isBkg && addUnc)  h->Add (sub, -1);
          else                        h->Reset ();
        }


        ////******** Do subtraction of binned dPhi plots ********//
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  sub = bkg->h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent];
        //  SaferDelete (&(h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]));
        //  h = (TH1D*) h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        //  h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent] = h;
        //  if (iCent != 0 || subtractPP) {
        //    if (!isBkg && !addUnc)      AddNoErrors (h, sub, -1);
        //    else if (!isBkg && addUnc)  h->Add (sub, -1);
        //    else                        h->Reset ();
        //  }

        //  sub = bkg->h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent];
        //  SaferDelete (&(h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]));
        //  h = (TH1D*) h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        //  h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent] = h;
        //  if (iCent != 0 || subtractPP) {
        //    if (!isBkg && !addUnc)      AddNoErrors (h, sub, -1);
        //    else if (!isBkg && addUnc)  h->Add (sub, -1);
        //    else                        h->Reset ();
        //  }
        //} // end loop over iPhi

      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  //******** Do background subtraction of phi distributions ********//
  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

        TH1D* sub = nullptr, *h = nullptr;
        for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
          sub = bkg->h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
          SaferDelete (&(h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]));

          h = (TH1D*) h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->Clone (Form ("h_trk_dphi_sub_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()));
          h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent] = h;
          if (iCent != 0 || subtractPP) {
            if (!isBkg && !addUnc)      AddNoErrors (h, sub, -1);
            else if (!isBkg && addUnc)  h->Add (sub, -1);
            else                        h->Reset ();
          }
        } // end loop over iPtch

        sub = bkg->h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent];
        SaferDelete (&(h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent]));

        h = (TH1D*) h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent]->Clone (Form ("h_trk_dphi_sub_gt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent] = h;
        if (iCent != 0 || subtractPP) {
          if (!isBkg && !addUnc)      AddNoErrors (h, sub, -1);
          else if (!isBkg && addUnc)  h->Add (sub, -1);
          else                        h->Reset ();
        }

        sub = bkg->h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent];
        SaferDelete (&(h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent]));

        h = (TH1D*) h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent]->Clone (Form ("h_trk_dphi_sub_gt4_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent] = h;
        if (iCent != 0 || subtractPP) {
          if (!isBkg && !addUnc)      AddNoErrors (h, sub, -1);
          else if (!isBkg && addUnc)  h->Add (sub, -1);
          else                        h->Reset ();
        }

      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc
  backgroundSubtracted = true;

  UnfoldSubtractedYield ();

  CombineHists ();

  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Applies bin migration factors to subtracted yields
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: UnfoldSubtractedYield () {
  if (histsUnfolded) return;
  assert (histsLoaded);

  TFile* f_binMigrationFile = new TFile (Form ("%s/BinMigrationFactors/binmigration_corrfactors_master.root", rootPath.Data ()), "read");

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : "mumu");
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        if (iCent != 0 && !unfoldPbPb) continue;

        const short iBBBCent = (iCent == 0 ? 0 : GetBBBCorrCentBin (0.5 * (centBins[iCent-1] + centBins[iCent])));

        const char* cent = (iBBBCent == 0 ? "pp" : Form ("CENT%i", numBBBCorrCentBins-iBBBCent-1));

        TF1* f = (TF1*) f_binMigrationFile->Get (Form ("tf1_%s_pt_ZPT%i_%s", spc, iPtZ-2, cent));
        TH1D* h = h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];
        for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
          const float x = h->GetBinCenter (ix);
          float y = h->GetBinContent (ix);
          float yerr = h->GetBinError (ix);
          y = y / f->Eval (x);
          yerr = yerr / f->Eval (x);
          h->SetBinContent (ix, y);
          h->SetBinError (ix, yerr);
        }
        TH2D* h2 = h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent];
        for (int ix = 1; ix <= h2->GetNbinsX (); ix++) {
          for (int iy = 1; iy <= h2->GetNbinsY (); iy++) {
            const double x = h2->GetXaxis ()->GetBinCenter (ix);
            const double y = h2->GetYaxis ()->GetBinCenter (iy);
            double z = h2->GetBinContent (ix, iy);
            z = z / (f->Eval (x) * f->Eval (y));
            h2->SetBinContent (ix, iy, z);
          }
        }

        f = f_trk_xhz_ptz_binMigration[iSpc][iPtZ][iCent] = (TF1*) f_binMigrationFile->Get (Form ("tf1_%s_xh_ZPT%i_%s", spc, iPtZ-2, cent));
        h = h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
        for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
          const float x = h->GetBinCenter (ix);
          float y = h->GetBinContent (ix);
          float yerr = h->GetBinError (ix);
          y = y / f->Eval (x);
          yerr = yerr / f->Eval (x);
          h->SetBinContent (ix, y);
          h->SetBinError (ix, yerr);
        }
        h2 = h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent];
        for (int ix = 1; ix <= h2->GetNbinsX (); ix++) {
          for (int iy = 1; iy <= h2->GetNbinsY (); iy++) {
            const double x = h2->GetXaxis ()->GetBinCenter (ix);
            const double y = h2->GetYaxis ()->GetBinCenter (iy);
            double z = h2->GetBinContent (ix, iy);
            z = z / (f->Eval (x) * f->Eval (y));
            h2->SetBinContent (ix, iy, z);
          }
        }
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  f_binMigrationFile->Close ();
  histsUnfolded = true;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CombineHists () {

  //**** Create empty subtracted histograms for combined channel yield measurement ****//
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
      h_trk_pt_ptz[2][iPtZ][iCent]->Reset ();
      h_trk_xhz_ptz[2][iPtZ][iCent]->Reset ();
      if (hasBkg && backgroundSubtracted) {
        h_trk_pt_ptz_sub[2][iPtZ][iCent]  = new TH1D (Form ("h_trk_pt_ptz_sub_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "", maxNPtchBins, allPtchBins);
        h_trk_xhz_ptz_sub[2][iPtZ][iCent] = new TH1D (Form ("h_trk_xhz_ptz_sub_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "", maxNXhZBins, allXhZBins);
        h_trk_pt_ptz_sub[2][iPtZ][iCent]->Sumw2 ();
        h_trk_xhz_ptz_sub[2][iPtZ][iCent]->Sumw2 ();
      }
      //for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
      //  h_trk_pt_dphi[2][iPtZ][iPhi][iCent]->Reset ();
      //  h_trk_xhz_dphi[2][iPtZ][iPhi][iCent]->Reset ();
      //  if (hasBkg && backgroundSubtracted) {
      //    h_trk_pt_dphi_sub[2][iPtZ][iPhi][iCent]  = new TH1D (Form ("h_trk_pt_dphi_sub_comb_iPtZ%i_iPhi%i_iCent%i_%s", iPtZ, iPhi, iCent, name.c_str ()), "", maxNPtchBins, allPtchBins);
      //    h_trk_xhz_dphi_sub[2][iPtZ][iPhi][iCent] = new TH1D (Form ("h_trk_xhz_dphi_sub_comb_iPtZ%i_iPhi%i_iCent%i_%s", iPtZ, iPhi, iCent, name.c_str ()), "", maxNXhZBins, allXhZBins);
      //    h_trk_pt_dphi_sub[2][iPtZ][iPhi][iCent]->Sumw2 ();
      //    h_trk_xhz_dphi_sub[2][iPtZ][iPhi][iCent]->Sumw2 ();
      //  }
      //} // end loop over iPhi
    } // end loop over iPtZ
  } // end loop over iCent


  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
        h_trk_dphi[2][iPtZ][iPtch][iCent] = new TH1D (Form ("h_trk_dphi_comb_iPtZ%i_iPtch%i_iCent%i_%s", iPtZ, iPtch, iCent, name.c_str ()), "", GetNdPhiBins (0, 0), -pi/2, 3*pi/2);
        h_trk_dphi[2][iPtZ][iPtch][iCent]->Sumw2 ();
        if (hasBkg && backgroundSubtracted) {
          h_trk_dphi_sub[2][iPtZ][iPtch][iCent] = new TH1D (Form ("h_trk_dphi_sub_comb_iPtZ%i_iPtch%i_iCent%i_%s", iPtZ, iPtch, iCent, name.c_str ()), "", GetNdPhiBins (0, 0), -pi/2, 3*pi/2);
          h_trk_dphi_sub[2][iPtZ][iPtch][iCent]->Sumw2 ();
        }
      } // end loop over iPtch
    } // end loop over iCent

    {
      short iCent = numCentBins;
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent] = new TH1D (Form ("h_trk_dphi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()), "", GetNdPhiBins (0, 0), -pi/2, 3*pi/2);
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->Sumw2 ();
          if (hasBkg && backgroundSubtracted) {
            h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent] = new TH1D (Form ("h_trk_dphi_sub_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()), "", GetNdPhiBins (0, 0), -pi/2, 3*pi/2);
            h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]->Sumw2 ();
          }
        } // end loop over iPtch
      } // end loop over iSpc
    } // end iCent = numCentBins scope
  } // end loop over iPtZ

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      const double totalSumWgts = h_z_counts[2][iPtZ][iCent]->GetBinContent (2);

      for (short iSpc = 0; iSpc < 2; iSpc++) {
        // Gets the weighting factor needed for this species.
        // E.g. if there are 2 muon events and 1 electron event,
        // the per-Z yield should be weighted by 2/3 in the muon
        // channel and 1/3 in the electron channel.
        const double channelCounts = h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1);
        const double channelSumWgts = h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (2);
        
        if (channelCounts == 0) {
          cout << "Warning: In PhysicsAnalysis :: CombineHists: Found 0 total Z bosons in this bin, iCent = " << iCent << ", iPtZ = " << iPtZ << ", iSpc = " << iSpc << "; weight is set to 0!" << endl;
          continue;
        }

        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  h_trk_pt_dphi[2][iPtZ][iPhi][iCent]->Add  (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent],  channelSumWgts/totalSumWgts);
        //  h_trk_xhz_dphi[2][iPtZ][iPhi][iCent]->Add (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent], channelSumWgts/totalSumWgts);
        //  if (hasBkg && backgroundSubtracted) {
        //    h_trk_pt_dphi_sub[2][iPtZ][iPhi][iCent]->Add  (h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent],  channelSumWgts/totalSumWgts);
        //    h_trk_xhz_dphi_sub[2][iPtZ][iPhi][iCent]->Add (h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent], channelSumWgts/totalSumWgts);
        //  }
        //} // end loop over iPhi

        h_trk_pt_ptz[2][iPtZ][iCent]->Add  (h_trk_pt_ptz[iSpc][iPtZ][iCent],  channelSumWgts/totalSumWgts);
        h_trk_xhz_ptz[2][iPtZ][iCent]->Add (h_trk_xhz_ptz[iSpc][iPtZ][iCent], channelSumWgts/totalSumWgts);
        if (hasBkg && backgroundSubtracted) {
          h_trk_pt_ptz_sub[2][iPtZ][iCent]->Add  (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent],  channelSumWgts/totalSumWgts);
          h_trk_xhz_ptz_sub[2][iPtZ][iCent]->Add (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent], channelSumWgts/totalSumWgts);
        }
      } // end loop over iSpc
    } // end loop over iPtZ
  } // end loop over iCent

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      const double totalSumWgts = h_z_counts[2][iPtZ][iCent]->GetBinContent (2);

      for (short iSpc = 0; iSpc < 2; iSpc++) {
        const double channelSumWgts = h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (2);
        
        for (int iPtch = 0; iPtch < maxNPtchBins+2; iPtch++) {
          h_trk_dphi[2][iPtZ][iPtch][iCent]->Add (h_trk_dphi[iSpc][iPtZ][iPtch][iCent], channelSumWgts/totalSumWgts);
          if (hasBkg && backgroundSubtracted) {
            h_trk_dphi_sub[2][iPtZ][iPtch][iCent]->Add (h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent], channelSumWgts/totalSumWgts);
          }
        } // end loop over iPtch
      } // end loop over iSpc
    } // end loop over iPtZ
  } // end loop over iCent

  
  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
    const double totalSumWgts = h_z_counts[iSpc][iPtZ][2]->GetBinContent (2) + h_z_counts[iSpc][iPtZ][3]->GetBinContent (2);

      for (short iCent = 2; iCent < numCentBins; iCent++) {
        const double channelSumWgts = h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (2);
      
        for (int iPtch = 0; iPtch < maxNPtchBins+2; iPtch++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][numCentBins]->Add (h_trk_dphi[iSpc][iPtZ][iPtch][iCent], channelSumWgts/totalSumWgts);
          if (hasBkg && backgroundSubtracted) {
            h_trk_dphi_sub[iSpc][iPtZ][iPtch][numCentBins]->Add (h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent], channelSumWgts/totalSumWgts);
          }
        } // end loop over iPtch
      } // end loop over iCent
    } // end loop over iSpc
  } // end loop over iPtZ

  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates signal-to-bkg ratios
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CalculateSigToBkg () {
  if (sigToBkgCalculated) return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

        TH1D* bkgYield = nullptr, *sigYield = nullptr;

        bkgYield = bkg->h_trk_pt_ptz[iSpc][iPtZ][iCent];
        SaferDelete (&(h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent]));
        sigYield = (TH1D*) h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent] = sigYield;
        if (iCent != 0 || subtractPP) sigYield->Divide (bkgYield);

        bkgYield = bkg->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
        SaferDelete (&(h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]));
        sigYield = (TH1D*) h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_sigToBkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent] = sigYield;
        if (iCent != 0 || subtractPP) sigYield->Divide (bkgYield);


        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  bkgYield = bkg->h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent];
        //  SaferDelete (&(h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]));
        //  sigYield = (TH1D*) h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_sigToBkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        //  h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = sigYield;
        //  if (iCent != 0 || subtractPP) sigYield->Divide (bkgYield);

        //  bkgYield = bkg->h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent];
        //  SaferDelete (&(h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]));
        //  sigYield = (TH1D*) h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_sigToBkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        //  h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = sigYield;
        //  if (iCent != 0 || subtractPP) sigYield->Divide (bkgYield);
        //} // end loop over iPhi
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent <= numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
        TH1D* bkgYield = nullptr, *sigYield = nullptr;

        bkgYield = bkg->h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent];
        SaferDelete (&(h_trk_dphi_sig_to_bkg[iSpc][iPtZ][maxNPtchBins][iCent]));
        sigYield = (TH1D*) h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent]->Clone (Form ("h_trk_dphi_sigToBkg_%s_iPtZ%i_gt4_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_dphi_sig_to_bkg[iSpc][iPtZ][maxNPtchBins][iCent] = sigYield;
        if (iCent != 0 || subtractPP) sigYield->Divide (bkgYield);

        bkgYield = bkg->h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent];
        SaferDelete (&(h_trk_dphi_sig_to_bkg[iSpc][iPtZ][maxNPtchBins+1][iCent]));
        sigYield = (TH1D*) h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent]->Clone (Form ("h_trk_dphi_sigToBkg_%s_iPtZ%i_lt4_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        h_trk_dphi_sig_to_bkg[iSpc][iPtZ][maxNPtchBins+1][iCent] = sigYield;
        if (iCent != 0 || subtractPP) sigYield->Divide (bkgYield);
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  sigToBkgCalculated = true;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates mean <xhZ> and <pTch> for each centrality as a function of pT^Z
// (Requires a set of pT^Z spectra to calculate <pT^Z> along the x-axis.)
// Also takes in a background estimate (optional) to incorporate the background statistical
// uncertainty matrices for systematics.
// If meanVar is specified, the bin centers used in the integration will be shifted by
// meanVar * meanErr
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CalculateTrackMeans (PhysicsAnalysis* nom, TH1D*** h_zpt_ptr, PhysicsAnalysis* cov, const double meanVar) {
  if (trackMeansCalculated) {
    cout << "Track means already calculated for " << name << ", returning" << endl; 
    return;
  }

  double zpt_norm, mean_zpt, mean_zpt_err;
  double*** mean_zpts = Get3DArray <double> (3, numCentBins, nPtZBins);
  double*** mean_zpt_errs = Get3DArray <double> (3, numCentBins, nPtZBins);
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : "mumu");

    for (short iCent = 0; iCent < numCentBins; iCent++) {

      //g_trk_avg_pt_ptz[iSpc][iCent] = new TGAE ();
      //g_trk_avg_pt_ptz[iSpc][iCent]->SetName (Form ("g_trk_avg_pt_ptz_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      //g_trk_avg_xhz_ptz[iSpc][iCent] = new TGAE ();
      //g_trk_avg_xhz_ptz[iSpc][iCent]->SetName (Form ("g_trk_avg_xhz_ptz_%s_iCent%i_%s", spc, iCent, name.c_str ()));

      TH1D* h_zpt = h_zpt_ptr[iCent][iSpc];
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        zpt_norm = 0;
        mean_zpt = 0;
        mean_zpt_err = 0;
        const int firstBin = h_zpt->FindBin (zPtBins[iPtZ]);
        const int lastBin = h_zpt->FindBin(zPtBins[iPtZ+1])-1;
        // calculate mean pT^Z
        for (int iX = firstBin; iX <= lastBin; iX++) {
          zpt_norm += h_zpt->GetBinContent (iX) * h_zpt->GetBinWidth (iX);
          mean_zpt += h_zpt->GetBinContent (iX) * h_zpt->GetBinCenter (iX) * h_zpt->GetBinWidth (iX);
        }
        mean_zpt = mean_zpt / zpt_norm;
        // now calculate the derivative vector and assume the covariance matrix is diagonal (events are independent)
        mean_zpt_err = 0;
        for (int iX = firstBin; iX <= lastBin; iX++) {
          mean_zpt_err += pow ((h_zpt->GetBinCenter (iX) - mean_zpt) * h_zpt->GetBinWidth (iX), 2) * pow (h_zpt->GetBinError (iX), 2) / pow (zpt_norm, 2);
        }
        mean_zpt_err = sqrt (mean_zpt_err);

        mean_zpts[iSpc][iCent][iPtZ] = mean_zpt;
        mean_zpt_errs[iSpc][iCent][iPtZ] = mean_zpt_err;
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  double track_norm, track_norm_err, mean_track, mean_track_err;
  double* deriv_vec = Get1DArray <double> (max (maxNPtchBins, maxNXhZBins));
  for (short iSpc : {0, 1}) {
    const char* spc = (iSpc == 0 ? "ee" : "mumu");
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      g_trk_avg_pt_ptz[iSpc][iCent] = new TGAE ();
      g_trk_avg_pt_ptz[iSpc][iCent]->SetName (Form ("g_trk_avg_pt_ptz_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      g_trk_avg_xhz_ptz[iSpc][iCent] = new TGAE ();
      g_trk_avg_xhz_ptz[iSpc][iCent]->SetName (Form ("g_trk_avg_xhz_ptz_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      g_avg_ntrk_ptz[iSpc][iCent] = new TGAE ();
      g_avg_ntrk_ptz[iSpc][iCent]->SetName (Form ("g_avg_ntrk_ptz_%s_iCent%i_%s", spc, iCent, name.c_str ()));

      TGAE* g = nullptr;
      TH1D* h = nullptr;
      TH2D* h2 = nullptr; // (co)variance matrix

      g = g_trk_avg_pt_ptz[iSpc][iCent];
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        h = nom->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]; // normalized per dpT^ch
        h2 = (cov == nullptr ? nom->h2_trk_pt_ptz_cov : cov->h2_trk_pt_ptz_cov)[iSpc][iPtZ][iCent];
        assert (h->GetNbinsX () == h2->GetNbinsX ());

        TF1* fit = new TF1 ("fit", "[0] * pow (x, [1])", pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
        fit->SetParameter (0, 1);
        fit->SetParameter (1, -5);
        h->Fit (fit, "RN0Q");

        double* fitParams = new double[2];

        fitParams[0] = fit->GetParameter (0);
        fitParams[1] = fit->GetParameter (1);

        SaferDelete (&fit);

        fit = new TF1 (Form ("f_trk_fit_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "[0] * pow (x, [1])", pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
        f_trk_fit_pt_ptz[iSpc][iPtZ][iCent] = fit;
        fit->SetParameter (0, fitParams[0]);
        fit->SetParameter (1, fitParams[1]);
        h->Fit (fit, "RN0Q");

        delete[] fitParams;

        const double gamma = fit->GetParameter (1);
        const double gamma1 = gamma+1.;
        const double gamma2 = gamma+2.;
        const double gammaErr = fit->GetParError (1);

        double meanPtch[maxNPtchBins] = {};
        double meanPtchErr[maxNPtchBins] = {};
        for (int iX = 0; iX < h->GetNbinsX (); iX++) {
          //double num = 0, numErr = 0, den = 0, denErr = 0;

          const double x1 = allPtchBins[iX];
          const double x2 = allPtchBins[iX+1];

          const double x2x1g1 = pow (x2, gamma1) - pow (x1, gamma1);
          const double x2x1g2 = pow (x2, gamma2) - pow (x1, gamma2);

          meanPtch[iX] = (gamma1/gamma2) * (x2x1g2 / x2x1g1);
          meanPtchErr[iX] = fabs ((x2x1g2 * x2x1g1 + log (x1/x2) * gamma1*gamma2* (x2-x1) * pow (x2, gamma1)*pow (x1, gamma1)) / (pow (gamma2 * x2x1g1, 2))) * fabs (gammaErr);
        } // end loop over iX

        track_norm = 0;
        mean_track = 0;
        // calculate mean pT^ch
        for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
          if (h->GetBinCenter (iX) < 2) continue;
          track_norm += h->GetBinContent (iX) * h->GetBinWidth (iX);
          //mean_track += h->GetBinContent (iX) * (h->GetBinCenter (iX) + meanVar * h->GetBinWidth (iX)) * h->GetBinWidth (iX);
          mean_track += h->GetBinContent (iX) * (meanPtch[iX-1] + meanVar * meanPtchErr[iX-1]) * h->GetBinWidth (iX);
        } // end loop over iX
        mean_track = mean_track / track_norm;

        // Now calculate derivative vector...
        for (int iX = 1; iX <= h->GetNbinsX (); iX++)
          deriv_vec[iX-1] = ((meanPtch[iX-1] - mean_track) * h->GetBinWidth (iX)) / track_norm;
        // ... and take its quadratic form with the error matrix to get the final variance
        mean_track_err = 0;
        for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
          if (h->GetBinCenter (iX) < 2) continue;
          for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
            if (h->GetBinCenter (iY) < 2) continue;
            mean_track_err += deriv_vec[iX-1] * h2->GetBinContent (iX, iY) * deriv_vec[iY-1];
          } // end loop over iY
        } // end loop over iX
        mean_track_err = sqrt (mean_track_err);

        g->SetPoint (g->GetN (), mean_zpts[iSpc][iCent][iPtZ], mean_track);
        g->SetPointError (g->GetN () - 1, mean_zpt_errs[iSpc][iCent][iPtZ], mean_zpt_errs[iSpc][iCent][iPtZ], mean_track_err, mean_track_err);
      } // end loop over iPtZ


      g = g_avg_ntrk_ptz[iSpc][iCent];
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        h = nom->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]; // normalized per dpT^ch
        h2 = (cov == nullptr ? nom->h2_trk_pt_ptz_cov : cov->h2_trk_pt_ptz_cov)[iSpc][iPtZ][iCent];
        assert (h->GetNbinsX () == h2->GetNbinsX ());

        track_norm = 0;
        // calculate mean pT^ch
        for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
          if (h->GetBinCenter (iX) < 2) continue;
          track_norm += h->GetBinContent (iX) * h->GetBinWidth (iX);
        } // end loop over iX

        // Now calculate uncertainty on total yield
        track_norm_err = 0;
        for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
          if (h->GetBinCenter (iX) < 2) continue;
          for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
            if (h->GetBinCenter (iY) < 2) continue;
            track_norm_err += (h2->GetBinContent (iX, iY)) * (h2->GetXaxis ()->GetBinWidth (iX)) * (h2->GetYaxis ()->GetBinWidth (iY));
          } // end loop over iY
        } // end loop over iX
        assert (track_norm_err > 0);
        track_norm_err = sqrt (track_norm_err);

        g = g_avg_ntrk_ptz[iSpc][iCent];
        g->SetPoint (g->GetN (), mean_zpts[iSpc][iCent][iPtZ], track_norm*(pi/4.));
        g->SetPointError (g->GetN () - 1, mean_zpt_errs[iSpc][iCent][iPtZ], mean_zpt_errs[iSpc][iCent][iPtZ], track_norm_err*(pi/4.), track_norm_err*(pi/4.));
      } // end loop over iPtZ


      g = g_trk_avg_xhz_ptz[iSpc][iCent];
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        h = nom->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
        h2 = (cov == nullptr ? nom->h2_trk_xhz_ptz_cov : cov->h2_trk_xhz_ptz_cov)[iSpc][iPtZ][iCent];
        assert (h->GetNbinsX () == h2->GetNbinsX ());

        TF1* fit = new TF1 ("fit", "[0] * pow (x, [1])", xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
        fit->SetParameter (0, 1);
        fit->SetParameter (1, -5);
        h->Fit (fit, "RN0Q");

        double* fitParams = new double[2];

        fitParams[0] = fit->GetParameter (0);
        fitParams[1] = fit->GetParameter (1);

        SaferDelete (&fit);

        fit = new TF1 (Form ("f_trk_fit_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "[0] * pow (x, [1])", xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
        f_trk_fit_xhz_ptz[iSpc][iPtZ][iCent] = fit;
        fit->SetParameter (0, fitParams[0]);
        fit->SetParameter (1, fitParams[1]);
        h->Fit (fit, "RN0Q");

        delete[] fitParams;

        const double gamma = fit->GetParameter (1);
        const double gamma1 = gamma+1.;
        const double gamma2 = gamma+2.;
        const double gammaErr = fit->GetParError (1);

        double meanXhZ[maxNXhZBins] = {};
        double meanXhZErr[maxNXhZBins] = {};
        for (int iX = 0; iX < h->GetNbinsX (); iX++) {
          //double num = 0, numErr = 0, den = 0, denErr = 0;

          const double x1 = allXhZBins[iX];
          const double x2 = allXhZBins[iX+1];

          const double x2x1g1 = pow (x2, gamma1) - pow (x1, gamma1);
          const double x2x1g2 = pow (x2, gamma2) - pow (x1, gamma2);

          meanXhZ[iX] = (gamma1/gamma2) * (x2x1g2 / x2x1g1);
          meanXhZErr[iX] = fabs ((x2x1g2 * x2x1g1 + log (x1/x2) * gamma1*gamma2* (x2-x1) * pow (x2, gamma1)*pow (x1, gamma1)) / (pow (gamma2 * x2x1g1, 2))) * fabs (gammaErr);
        } // end loop over iX
       
        track_norm = 0;
        mean_track = 0;
        // calculate mean x_hZ
        for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
          if (h->GetBinCenter (iX) < 1./15.) continue;
          track_norm += h->GetBinContent (iX) * h->GetBinWidth (iX);
          mean_track += h->GetBinContent (iX) * (meanXhZ[iX-1] + meanVar * meanXhZErr[iX-1]) * h->GetBinWidth (iX);
        } // end loop over iX
        mean_track = mean_track / track_norm;

        // now calculate derivative vector...
        for (int iX = 1; iX <= h->GetNbinsX (); iX++)
          deriv_vec[iX-1] = ((meanXhZ[iX-1] - mean_track) * h->GetBinWidth (iX)) / track_norm;
        // ... and take its quadratic form with the error matrix to get the final variance
        mean_track_err = 0;
        for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
          if (h->GetBinCenter (iX) < 1./15.) continue;
          for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
            if (h->GetBinCenter (iY) < 1./15.) continue;
            mean_track_err += deriv_vec[iX-1] * h2->GetBinContent (iX, iY) * deriv_vec[iY-1];
          } // end loop over iY
        } // end loop over iX
        mean_track_err = sqrt (mean_track_err);

        g->SetPoint (g->GetN (), mean_zpts[iSpc][iCent][iPtZ], mean_track);
        g->SetPointError (g->GetN () - 1, mean_zpt_errs[iSpc][iCent][iPtZ], mean_zpt_errs[iSpc][iCent][iPtZ], mean_track_err, mean_track_err);

      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  Delete1DArray (&deriv_vec, max (maxNPtchBins, maxNXhZBins));

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    g_trk_avg_pt_ptz[2][iCent] = new TGAE ();
    g_trk_avg_pt_ptz[2][iCent]->SetName (Form ("g_trk_avg_pt_ptz_comb_iCent%i_%s", iCent, name.c_str ()));
    g_trk_avg_xhz_ptz[2][iCent] = new TGAE ();
    g_trk_avg_xhz_ptz[2][iCent]->SetName (Form ("g_trk_avg_xhz_ptz_comb_iCent%i_%s", iCent, name.c_str ()));
    g_avg_ntrk_ptz[2][iCent] = new TGAE ();
    g_avg_ntrk_ptz[2][iCent]->SetName (Form ("g_avg_ntrk_ptz_comb_iCent%i_%s", iCent, name.c_str ()));

    TGAE* g = nullptr;

    g = g_trk_avg_pt_ptz[2][iCent];

    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
      mean_track = 0;
      mean_track_err = 0;
      double xval = 0, yval = 0;
      double sumWeights = 0;
      for (int iSpc : {0, 1}) {
        TGAE* gspc = g_trk_avg_pt_ptz[iSpc][iCent];
        const double channelWeight = pow (gspc->GetErrorY (iPtZ-2), -2);
        gspc->GetPoint (iPtZ-2, xval, yval);
        mean_track += yval * channelWeight;
        sumWeights += channelWeight;
      } // end loop over iSpc
      mean_track = mean_track / sumWeights;
      mean_track_err = sqrt (1/sumWeights);

      g->SetPoint (g->GetN (), mean_zpts[2][iCent][iPtZ], mean_track);
      g->SetPointError (g->GetN () - 1, mean_zpt_errs[2][iCent][iPtZ], mean_zpt_errs[2][iCent][iPtZ], mean_track_err, mean_track_err);
    } // end loop over iPtZ

    g = g_avg_ntrk_ptz[2][iCent];

    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
      mean_track = 0;
      mean_track_err = 0;
      double xval = 0, yval = 0;
      double sumWeights = 0;
      for (int iSpc : {0, 1}) {
        TGAE* gspc = g_avg_ntrk_ptz[iSpc][iCent];
        const double channelWeight = pow (gspc->GetErrorY (iPtZ-2), -2);
        gspc->GetPoint (iPtZ-2, xval, yval);
        mean_track += yval * channelWeight;
        sumWeights += channelWeight;
      } // end loop over iSpc
      mean_track = mean_track / sumWeights;
      mean_track_err = sqrt (1/sumWeights);

      g->SetPoint (g->GetN (), mean_zpts[2][iCent][iPtZ], mean_track);
      g->SetPointError (g->GetN () - 1, mean_zpt_errs[2][iCent][iPtZ], mean_zpt_errs[2][iCent][iPtZ], mean_track_err, mean_track_err);
    } // end loop over iPtZ

    g = g_trk_avg_xhz_ptz[2][iCent];
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
      mean_track = 0;
      mean_track_err = 0;
      double xval = 0, yval = 0;
      double sumWeights = 0;
      for (int iSpc : {0, 1}) {
        TGAE* gspc = g_trk_avg_xhz_ptz[iSpc][iCent];
        const double channelWeight = pow (gspc->GetErrorY (iPtZ-2), -2);
        gspc->GetPoint (iPtZ-2, xval, yval);
        mean_track += yval * channelWeight;
        sumWeights += channelWeight;
      }
      mean_track = mean_track / sumWeights;
      mean_track_err = sqrt (1/sumWeights);

      g->SetPoint (g->GetN (), mean_zpts[2][iCent][iPtZ], mean_track);
      g->SetPointError (g->GetN () - 1, mean_zpt_errs[2][iCent][iPtZ], mean_zpt_errs[2][iCent][iPtZ], mean_track_err, mean_track_err);
    } // end loop over iPtZ
  } // end loop over iCent

  Delete3DArray (&mean_zpts, 3, numCentBins, nPtZBins);
  Delete3DArray (&mean_zpt_errs, 3, numCentBins, nPtZBins);

  trackMeansCalculated = true;

  TruncatePhysicsPlots ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Removes excess bins from plots which no longer need to be considered.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: TruncatePhysicsPlots () {
  cout << "Truncating extra bins for " << name << endl;
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        //for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  TruncateTH1D (&(h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]), iPtZ, true);
        //  TruncateTH1D (&(h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]), iPtZ, true);
        //  TruncateTH1D (&(h_trk_pt_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]), iPtZ, true);
        //  TruncateTH1D (&(h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent]), iPtZ, true);
        //  TruncateTH1D (&(h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]), iPtZ, false);
        //  TruncateTH1D (&(h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]), iPtZ, false);
        //  TruncateTH1D (&(h_trk_xhz_dphi_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]), iPtZ, false);
        //  TruncateTH1D (&(h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]), iPtZ, false);
        //}

        TruncateTH1D (&(h_trk_pt_ptz[iSpc][iPtZ][iCent]), iPtZ, true);
        TruncateTH1D (&(h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]), iPtZ, true);
        TruncateTH1D (&(h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent]), iPtZ, true);
        TruncateTH1D (&(h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]), iPtZ, true);
        TruncateTH1D (&(h_trk_xhz_ptz[iSpc][iPtZ][iCent]), iPtZ, false);
        TruncateTH1D (&(h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]), iPtZ, false);
        TruncateTH1D (&(h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]), iPtZ, false);
        TruncateTH1D (&(h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]), iPtZ, false);
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Inflates statistical uncertainty in final results histograms by some amount
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: ProjectLumiIncrease () {

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h_pp = h_trk_pt_ptz_sub[2][iPtZ][0];
    for (short iCent = 1; iCent < numCentBins; iCent++) {

      const float n_ee = h_z_counts[0][iPtZ][iCent]->GetBinContent (2);
      const float n_mumu = h_z_counts[1][iPtZ][iCent]->GetBinContent (2);
      TH1D* h_ee, *h_mumu, *h_comb;
      h_ee = (TH1D*) h_trk_pt_ptz_sub[0][iPtZ][iCent]->Clone ("h_ee");
      h_mumu = (TH1D*) h_trk_pt_ptz_sub[1][iPtZ][iCent]->Clone ("h_mumu");
      h_comb = (TH1D*) h_trk_pt_ptz_sub[2][iPtZ][iCent]->Clone ("h_comb");
      h_comb->Reset ();

      h_ee->Scale (n_ee);
      h_mumu->Scale (n_mumu);

      const float n_ee_proj = n_ee * (2.19 / 1.7);
      const float n_mumu_proj = n_mumu * (1.89 / 1.4);

      for (int iX = 1; iX <= h_ee->GetNbinsX (); iX++) {
        h_ee->SetBinContent (iX, h_ee->GetBinContent (iX) * n_ee_proj / n_ee);
        h_ee->SetBinError (iX, h_ee->GetBinError (iX) * sqrt (n_ee_proj / n_ee));
      }
      for (int iX = 1; iX <= h_mumu->GetNbinsX (); iX++) {
        h_mumu->SetBinContent (iX, h_mumu->GetBinContent (iX) * n_mumu_proj / n_mumu);
        h_mumu->SetBinError (iX, h_mumu->GetBinError (iX) * sqrt (n_mumu_proj / n_mumu));
      }

      h_comb->Add (h_ee);
      h_comb->Add (h_mumu);
      h_comb->Scale (1. / (n_ee_proj+n_mumu_proj));

      h_trk_pt_ptz_iaa_2015proj[iPtZ][iCent] = (TH1D*) h_trk_pt_ptz_iaa[2][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_iaa_2015proj_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      TH1D* h = h_trk_pt_ptz_iaa_2015proj[iPtZ][iCent];

      for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
        const float iaa = h_comb->GetBinContent (iX) / h_pp->GetBinContent (iX);
        const float iaa_err = iaa * sqrt (pow (h_comb->GetBinError (iX) / h_comb->GetBinContent (iX), 2) + pow (h_pp->GetBinError (iX) / h_pp->GetBinContent (iX), 2));
        h->SetBinError (iX, iaa_err);
      }
    } // end loop over iCent

    h_pp = h_trk_xhz_ptz_sub[2][iPtZ][0];
    for (short iCent = 1; iCent < numCentBins; iCent++) {

      const float n_ee = h_z_counts[0][iPtZ][iCent]->GetBinContent (2);
      const float n_mumu = h_z_counts[1][iPtZ][iCent]->GetBinContent (2);
      TH1D* h_ee, *h_mumu, *h_comb;
      h_ee = (TH1D*) h_trk_xhz_ptz_sub[0][iPtZ][iCent]->Clone ("h_ee");
      h_mumu = (TH1D*) h_trk_xhz_ptz_sub[1][iPtZ][iCent]->Clone ("h_mumu");
      h_comb = (TH1D*) h_trk_xhz_ptz_sub[2][iPtZ][iCent]->Clone ("h_comb");
      h_comb->Reset ();

      h_ee->Scale (n_ee);
      h_mumu->Scale (n_mumu);

      const float n_ee_proj = n_ee * (2.19 / 1.7);
      const float n_mumu_proj = n_mumu * (1.89 / 1.4);

      for (int iX = 1; iX <= h_ee->GetNbinsX (); iX++) {
        h_ee->SetBinContent (iX, h_ee->GetBinContent (iX) * n_ee_proj / n_ee);
        h_ee->SetBinError (iX, h_ee->GetBinError (iX) * sqrt (n_ee_proj / n_ee));
      }
      for (int iX = 1; iX <= h_mumu->GetNbinsX (); iX++) {
        h_mumu->SetBinContent (iX, h_mumu->GetBinContent (iX) * n_mumu_proj / n_mumu);
        h_mumu->SetBinError (iX, h_mumu->GetBinError (iX) * sqrt (n_mumu_proj / n_mumu));
      }

      h_comb->Add (h_ee);
      h_comb->Add (h_mumu);
      h_comb->Scale (1. / (n_ee_proj+n_mumu_proj));

      h_trk_xhz_ptz_iaa_2015proj[iPtZ][iCent] = (TH1D*) h_trk_xhz_ptz_iaa[2][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_iaa_2015proj_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      TH1D* h = h_trk_xhz_ptz_iaa_2015proj[iPtZ][iCent];

      for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
        const float iaa = h_comb->GetBinContent (iX) / h_pp->GetBinContent (iX);
        const float iaa_err = iaa * sqrt (pow (h_comb->GetBinError (iX) / h_comb->GetBinContent (iX), 2) + pow (h_pp->GetBinError (iX) / h_pp->GetBinContent (iX), 2));
        h->SetBinError (iX, iaa_err);
      }
    } // end loop over iCent
  } // end loop over iPtZ
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Inflates statistical uncertainty in final results histograms by some amount
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: InflateStatUnc (const float amount) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        for (TH1D* h : {h_trk_pt_ptz[iSpc][iPtZ][iCent], h_trk_xhz_ptz[iSpc][iPtZ][iCent], h_trk_pt_ptz_sub[iSpc][iPtZ][iCent], h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent], h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent], h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent]}) {
          if (!h)
            continue;
          for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
            h->SetBinError (ix, h->GetBinError (ix) * (1+amount));
          }
        }
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Applies a constant (x-axis independent) systematic to the analysis, e.g. a luminosity unc.
// Only applies to pT dependent histograms, i.e. will need to subtract background still.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: ApplyRelativeVariation (float**** relVar, const bool upVar) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

      // Hadron yield systematics, signal & signal+bkg levels
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][iPhi][iCent] : 1/relVar[iSpc][iPtZ][iPhi][iCent]);
        //  h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][iPhi][iCent] : 1/relVar[iSpc][iPtZ][iPhi][iCent]);
        //} // end loop over iPhi

        h_trk_pt_ptz[iSpc][iPtZ][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][numPhiBins][iCent] : 1/relVar[iSpc][iPtZ][numPhiBins][iCent]);
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][numPhiBins][iCent] : 1/relVar[iSpc][iPtZ][numPhiBins][iCent]);
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Converts analysis to a systematic variation by adding or subtracting statistical errors
// Only converts relevant histograms, e.g. for minbias will only do this to track yields
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: ConvertToStatVariation (const bool upVar, const float nSigma) {
  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {

      // Hadron yield systematics, signal & signal+bkg levels
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  AddStatVar (h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
        //  AddStatVar (h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);

        //  AddStatVar (h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
        //  AddStatVar (h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent], upVar, nSigma);
        //} // end loop over iPhi

        AddStatVar (h_trk_pt_ptz[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_xhz_ptz[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_pt_ptz_sub[iSpc][iPtZ][iCent], upVar, nSigma);
        AddStatVar (h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent], upVar, nSigma);
      } // end loop over iCent
     
    } // end loop over iSpc

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (short iPtch = 0; iPtch < maxNPtchBins+2; iPtch++) {
          AddStatVar (h_trk_dphi[iSpc][iPtZ][iPtch][iCent], upVar, nSigma);
          //AddStatVar (h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent], upVar, nSigma);
        } // end loop over iPtch
      } // end loop over iCent
    } // end loop over iSpc
  } // end loop over iPtZ
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
// Designed to be overloaded. The default here is for analyzing data.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  cout << "Arguments provided: " << endl;
  cout << "inFileName = " << inFileName << endl;
  cout << "outFileName = " << outFileName << endl;

  LoadEventWeights ();
  //eventPlaneCalibrator = EventPlaneCalibrator (Form ("%s/FCalCalibration/Nominal/data18hi.root", rootPath.Data ()));

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1, fcal_weight = 1, q2_weight = 1, psi2_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0;
  float z_pt = 0, z_y = 0, z_phi = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  float l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];

  int***    trks_counts   = Get3DArray <int> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights1 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights2 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  int**     trks_counts_inPhi   = Get2DArray <int> (maxNPtchBins, 40);
  float**   trks_weights1_inPhi = Get2DArray <float> (maxNPtchBins, 40);
  float**   trks_weights2_inPhi = Get2DArray <float> (maxNPtchBins, 40);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("isEE",         &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("vz",           &vz);
    PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
    PbPbTree->SetBranchAddress ("z_y",          &z_y);
    PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
    PbPbTree->SetBranchAddress ("z_m",          &z_m);
    PbPbTree->SetBranchAddress ("l1_pt",        &l1_pt);
    PbPbTree->SetBranchAddress ("l1_eta",       &l1_eta);
    PbPbTree->SetBranchAddress ("l1_phi",       &l1_phi);
    PbPbTree->SetBranchAddress ("l1_charge",    &l1_charge);
    PbPbTree->SetBranchAddress ("l1_trk_pt",    &l1_trk_pt);
    PbPbTree->SetBranchAddress ("l1_trk_eta",   &l1_trk_eta);
    PbPbTree->SetBranchAddress ("l1_trk_phi",   &l1_trk_phi);
    PbPbTree->SetBranchAddress ("l2_pt",        &l2_pt);
    PbPbTree->SetBranchAddress ("l2_eta",       &l2_eta);
    PbPbTree->SetBranchAddress ("l2_phi",       &l2_phi);
    PbPbTree->SetBranchAddress ("l2_charge",    &l2_charge);
    PbPbTree->SetBranchAddress ("l2_trk_pt",    &l2_trk_pt);
    PbPbTree->SetBranchAddress ("l2_trk_eta",   &l2_trk_eta);
    PbPbTree->SetBranchAddress ("l2_trk_phi",   &l2_trk_phi);
    PbPbTree->SetBranchAddress ("ntrk",         &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",       &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",      &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",      &trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      if (fabs (vz) > 150) continue;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

      const short iCent = GetCentBin (fcal_et);
      if (iCent < 1 || iCent > numCentBins-1) continue;

      const short iFineCent = GetFineCentBin (fcal_et);
      if (iFineCent < 1 || iFineCent > numFineCentBins-1) continue;

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

      // do a reweighting procedure
      {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)  dphi = pi - dphi;
        if (useCentWgts)  fcal_weight = h_PbPbFCal_weights[iSpc][iPtZ]->GetBinContent (h_PbPbFCal_weights[iSpc][iPtZ]->FindBin (fcal_et));
        if (useQ2Wgts)    q2_weight   = h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->FindBin (q2));
        if (usePsi2Wgts)  psi2_weight = h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->FindBin (dphi));

        event_weight *= fcal_weight * q2_weight * psi2_weight;
      }

      if (event_weight == 0) continue;

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);

      h_q2[iFineCent]->Fill (q2);
      h_q2_reweighted[iFineCent]->Fill (q2, event_weight);
      h_psi2[iFineCent]->Fill (psi2);
      h_psi2_reweighted[iFineCent]->Fill (psi2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.02 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.02)) continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true);
        if (trkEff == 0 || trkPur == 0) continue;
        const float trkWeight = trkPur / trkEff;

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (allPtchBins[0] <= trkpt) {
          iPtch = 0;
          while (iPtch < maxNPtchBins && allPtchBins[iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < maxNPtchBins) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (0, 0) && (-pi/2.)+(2.*pi/GetNdPhiBins (0, 0))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (allXhZBins[0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < maxNXhZBins && allXhZBins[iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < maxNPtchBins) {
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < maxNXhZBins) {
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if (3*pi/4 <= dPhi) {
          if (iPtch != -1 && iPtch < maxNPtchBins) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < maxNXhZBins) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks

      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
        for (int idPhi1 = 0; idPhi1 < GetNdPhiBins (0, 0); idPhi1++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi1+1, sqrt (pow (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinError (idPhi1+1), 2) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi1])));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (0, 0); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi1

      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < max (maxNPtchBins, maxNXhZBins); j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < maxNPtchBins; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

    } // end loop over Pb+Pb tree
    cout << "Done primary Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("event_weight", &event_weight);
    ppTree->SetBranchAddress ("isEE",         &isEE);
    ppTree->SetBranchAddress ("vz",           &vz);
    ppTree->SetBranchAddress ("z_pt",         &z_pt);
    ppTree->SetBranchAddress ("z_y",          &z_y);
    ppTree->SetBranchAddress ("z_phi",        &z_phi);
    ppTree->SetBranchAddress ("z_m",          &z_m);
    ppTree->SetBranchAddress ("l1_pt",        &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",       &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",       &l1_phi);
    ppTree->SetBranchAddress ("l1_charge",    &l1_charge);
    ppTree->SetBranchAddress ("l1_trk_pt",    &l1_trk_pt);
    ppTree->SetBranchAddress ("l1_trk_eta",   &l1_trk_eta);
    ppTree->SetBranchAddress ("l1_trk_phi",   &l1_trk_phi);
    ppTree->SetBranchAddress ("l2_pt",        &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",       &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",       &l2_phi);
    ppTree->SetBranchAddress ("l2_charge",    &l2_charge);
    ppTree->SetBranchAddress ("l2_trk_pt",    &l2_trk_pt);
    ppTree->SetBranchAddress ("l2_trk_eta",   &l2_trk_eta);
    ppTree->SetBranchAddress ("l2_trk_phi",   &l2_trk_phi);
    ppTree->SetBranchAddress ("ntrk",         &ntrk);
    ppTree->SetBranchAddress ("trk_pt",       &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",      &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",      &trk_phi);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      if (fabs (vz) > 150) continue;

      if (event_weight == 0) continue;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.02 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.02)) continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false);
        if (trkEff == 0 || trkPur == 0) continue;
        const float trkWeight = trkPur / trkEff;

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (allPtchBins[0] <= trkpt) {
          iPtch = 0;
          while (iPtch < maxNPtchBins && allPtchBins[iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < maxNPtchBins) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (0, 0) && (-pi/2.)+(2.*pi/GetNdPhiBins (0, 0))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (allXhZBins[0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < maxNXhZBins && allXhZBins[iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < maxNPtchBins) {
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < maxNXhZBins) {
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if (3*pi/4 <= dPhi) {
          if (iPtch != -1 && iPtch < maxNPtchBins) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < maxNXhZBins) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks

      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
        for (int idPhi1 = 0; idPhi1 < GetNdPhiBins (0, 0); idPhi1++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi1+1, sqrt (pow (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinError (idPhi1+1), 2) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi1])));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (0, 0); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi1

      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < max (maxNPtchBins, maxNXhZBins); j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < maxNPtchBins; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

    } // end loop over pp tree
    cout << "Done primary pp loop." << endl;
  }

  Delete3DArray (&trks_counts, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights1, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights2, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete2DArray (&trks_counts_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights1_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights2_inPhi, maxNPtchBins, 40);

  SaveHists (outFileName);

  if (inFile) inFile->Close ();
  SaferDelete (&inFile);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Generates weights between a weighted sample and a matched sample.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: GenerateEventWeights (const char* weightedSampleInFilePattern, const char* matchedSampleInFilePattern, const char* outFileName) {

  const int nQ2Bins = 20;
  const double* q2Bins = linspace (0, 0.3, nQ2Bins);
  const int nPsi2Bins = 8;
  const double* psi2Bins = linspace (0, pi/2, nPsi2Bins);

  TH1D* h_fcal_et_dist[2][3][nPtZBins];
  TH1D* h_q2_dist[2][3][numFineCentBins][nPtZBins];
  TH1D* h_psi2_dist[2][3][numFineCentBins][nPtZBins];

  for (int iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      h_fcal_et_dist[0][iSpc][iPtZ] = new TH1D (Form ("h_weighted_fcal_et_dist_%s_iPtZ%i_%s", spc, iPtZ, name.c_str ()), "", numSuperFineCentBins-1, superFineCentBins);
      h_fcal_et_dist[0][iSpc][iPtZ]->Sumw2 ();
      h_fcal_et_dist[1][iSpc][iPtZ] = new TH1D (Form ("h_fcal_et_dist_%s_iPtZ%i_%s", spc, iPtZ, name.c_str ()), "", numSuperFineCentBins-1, superFineCentBins);
      h_fcal_et_dist[1][iSpc][iPtZ]->Sumw2 ();

      for (int iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
        h_q2_dist[0][iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_weighted_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nQ2Bins, q2Bins);
        h_q2_dist[0][iSpc][iFineCent][iPtZ]->Sumw2 ();
        h_q2_dist[1][iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nQ2Bins, q2Bins);
        h_q2_dist[1][iSpc][iFineCent][iPtZ]->Sumw2 ();

        h_psi2_dist[0][iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_weighted_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nPsi2Bins, psi2Bins);
        h_psi2_dist[0][iSpc][iFineCent][iPtZ]->Sumw2 ();
        h_psi2_dist[1][iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nPsi2Bins, psi2Bins);
        h_psi2_dist[1][iSpc][iFineCent][iPtZ]->Sumw2 ();
      }
    }
  }


  bool isEE = false;
  float event_weight = 0, fcal_et = 0, z_pt = 0, z_phi = 0, q2 = 0, psi2 = 0;

  TChain* PbPbTree = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  PbPbTree->Add (weightedSampleInFilePattern);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  PbPbTree->SetBranchAddress ("isEE",         &isEE);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
  PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);

  int nEvts = PbPbTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    PbPbTree->GetEntry (iEvt);

    const short iSpc = (isEE ? 0 : 1);

    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;

    h_fcal_et_dist[0][iSpc][iPtZ]->Fill (fcal_et);
  }
  cout << "Done 1st Pb+Pb loop over weighted sample." << endl;

  PbPbTree->Reset ();
  PbPbTree->Add (matchedSampleInFilePattern);

  PbPbTree->SetBranchAddress ("isEE",         &isEE);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
  PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);

  nEvts = PbPbTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    PbPbTree->GetEntry (iEvt);

    const short iSpc = (isEE ? 0 : 1);

    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;

    h_fcal_et_dist[1][iSpc][iPtZ]->Fill (fcal_et);
  }
  cout << "Done 1st Pb+Pb loop over matched events." << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Finalize FCal weighting histograms
  ////////////////////////////////////////////////////////////////////////////////////////////////
  for (short iSample : {0, 1})
    for (short iSpc = 0; iSpc < 2; iSpc++)
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++)
        h_fcal_et_dist[iSample][2][iPtZ]->Add (h_fcal_et_dist[iSample][iSpc][iPtZ]);

  for (short iSample : {0, 1})
    for (short iSpc = 0; iSpc < 3; iSpc++)
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++)
        h_fcal_et_dist[iSample][iSpc][iPtZ]->Scale (1/h_fcal_et_dist[iSample][iSpc][iPtZ]->Integral ());

  for (short iSpc = 0; iSpc < 3; iSpc++)
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++)
      h_fcal_et_dist[1][iSpc][iPtZ]->Divide (h_fcal_et_dist[0][iSpc][iPtZ]);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree for additional q2, psi2 weights
  ////////////////////////////////////////////////////////////////////////////////////////////////
  PbPbTree->Reset ();
  PbPbTree->Add (weightedSampleInFilePattern);

  PbPbTree->SetBranchAddress ("isEE",         &isEE);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
  PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);

  nEvts = PbPbTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    PbPbTree->GetEntry (iEvt);

    const short iSpc = (isEE ? 0 : 1);

    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;
    const short iFineCent = GetFineCentBin (fcal_et);
    if (iFineCent < 1 || iFineCent > numFineCentBins-1)
      continue;

    event_weight = h_fcal_et_dist[1][iSpc][iPtZ]->GetBinContent (h_fcal_et_dist[1][iSpc][iPtZ]->FindFixBin (fcal_et));

    h_q2_dist[0][iSpc][iFineCent][iPtZ]->Fill (q2, event_weight);

    float dphi = DeltaPhi (z_phi, psi2, false);
    if (dphi > pi/2)
      dphi = pi - dphi;
    h_psi2_dist[0][iSpc][iFineCent][iPtZ]->Fill (dphi);
  }
  cout << "Done 2nd Pb+Pb loop over weighted sample." << endl;

  PbPbTree->Reset ();
  PbPbTree->Add (matchedSampleInFilePattern);

  PbPbTree->SetBranchAddress ("isEE",         &isEE);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
  PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);

  nEvts = PbPbTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    PbPbTree->GetEntry (iEvt);

    const short iSpc = (isEE ? 0 : 1);

    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;
    const short iFineCent = GetFineCentBin (fcal_et);
    if (iFineCent < 1 || iFineCent > numFineCentBins-1)
      continue;

    h_q2_dist[1][iSpc][iFineCent][iPtZ]->Fill (q2);

    float dphi = DeltaPhi (z_phi, psi2, false);
    if (dphi > pi/2)
      dphi = pi - dphi;
    h_psi2_dist[1][iSpc][iFineCent][iPtZ]->Fill (dphi);
  }
  cout << "Done 2nd Pb+Pb loop over matched events." << endl;


  delete PbPbTree;


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Finalize weighting histograms
  ////////////////////////////////////////////////////////////////////////////////////////////////
  for (short iSample : {0, 1}) {
    for (short iSpc = 0; iSpc < 2; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
          h_q2_dist[iSample][2][iFineCent][iPtZ]->Add (h_q2_dist[iSample][iSpc][iFineCent][iPtZ]);
          h_psi2_dist[iSample][2][iFineCent][iPtZ]->Add (h_psi2_dist[iSample][iSpc][iFineCent][iPtZ]);
        }
      }
    }
  }

  for (short iSample : {0, 1}) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
          h_q2_dist[iSample][iSpc][iFineCent][iPtZ]->Scale (1/h_q2_dist[iSample][iSpc][iFineCent][iPtZ]->Integral ());
          h_psi2_dist[iSample][iSpc][iFineCent][iPtZ]->Scale (1/h_psi2_dist[iSample][iSpc][iFineCent][iPtZ]->Integral ());
        }
      }
    }
  }

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
        h_q2_dist[1][iSpc][iFineCent][iPtZ]->Divide (h_q2_dist[0][iSpc][iFineCent][iPtZ]);
        h_psi2_dist[1][iSpc][iFineCent][iPtZ]->Divide (h_psi2_dist[0][iSpc][iFineCent][iPtZ]);
      }
    }
  }


  if (eventWeightsFile && eventWeightsFile->IsOpen ()) {
    eventWeightsFile->Close ();
    eventWeightsLoaded = false;
  }

  eventWeightsFile = new TFile (outFileName, "recreate");
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      SafeWrite (h_fcal_et_dist[1][iSpc][iPtZ]);
      for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
        SafeWrite (h_q2_dist[1][iSpc][iFineCent][iPtZ]);
        SafeWrite (h_psi2_dist[1][iSpc][iFineCent][iPtZ]);
      }
    }
  }
  eventWeightsFile->Close ();
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Load the event weights into memory
//////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LoadEventWeights () {
  if (eventWeightsLoaded)
    return;

  if (!useCentWgts && !useQ2Wgts && !usePsi2Wgts)
    return;

  const TString ext = (ewExt == "" ? TString (name) : ewExt);

  TDirectory* _gDirectory = gDirectory;

  cout << "Loading event weights from " << rootPath.Data () << "/" << eventWeightsFileName.Data () << endl;
  eventWeightsFile = new TFile (Form ("%s/%s", rootPath.Data (), eventWeightsFileName.Data ()), "read");

  //if (useCentWgts)  h_PbPbFCal_weights = (TH1D*) eventWeightsFile->Get (Form ("h_PbPbFCal_weights_%s", name.c_str ()));
  //for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
  //  if (useQ2Wgts)    h_PbPbQ2_weights[iFineCent] = (TH1D*) eventWeightsFile->Get (Form ("h_PbPbQ2_weights_iCent%i_%s", iFineCent, name.c_str ()));
  //  if (usePsi2Wgts)  h_PbPbPsi2_weights[iFineCent] = (TH1D*) eventWeightsFile->Get (Form ("h_PbPbPsi2_weights_iCent%i_%s", iFineCent, name.c_str ()));
  //}
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (useCentWgts)  h_PbPbFCal_weights[iSpc][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_fcal_et_dist_%s_iPtZ%i_%s", spc, iPtZ, ext.Data ()));
      for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
        if (useQ2Wgts)    h_PbPbQ2_weights[iSpc][iFineCent][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, ext.Data ()));
        if (usePsi2Wgts)  h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, ext.Data ()));
      }
    }
  }

  cout << "Event weights loaded." << endl;
  eventWeightsLoaded = true;

  _gDirectory-> cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load the tracking efficiencies into memory
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LoadTrackingEfficiencies (const bool doRebin) {
  if (effsLoaded)
    return;

  TDirectory* _gDirectory = gDirectory;

  TString _effDir = "Nominal";
  if (useHijingEffs)
    _effDir = "Hijing";
  else if (useHITight)
    _effDir = "Variations/TrackHITightWPVariation";
  else if (doTrackEffVar)
    _effDir = "Variations/TrackEffPionsVariation";
  cout << Form ("Reading tracking efficiencies from %s/TrackingEfficiencies/%s/trackingEfficiencies_%s.root", rootPath.Data (), _effDir.Data (), is2015Conds ? "15" : "18") << endl;
  trkEffFile = new TFile (Form ("%s/TrackingEfficiencies/%s/trackingEfficiencies_%s.root", rootPath.Data (), _effDir.Data (), is2015Conds ? "15" : "18"), "read");

  for (int iCent = 0; iCent < numTrkCorrCentBins; iCent++) {
  //for (int iCent = 0; iCent < numFineCentBins; iCent++) {
    h2_num_trk_effs[iCent] = (TH2D*) trkEffFile->Get (Form ("h_truth_matched_reco_tracks_iCent%i", iCent));
    h2_den_trk_effs[iCent] = (TH2D*) trkEffFile->Get (Form ("h_truth_tracks_iCent%i", iCent));

    //if (iCent > 0) {
    //  h2_num_trk_effs[iCent]->RebinX (2);
    //  h2_num_trk_effs[iCent]->RebinY (2);
    //  h2_den_trk_effs[iCent]->RebinX (2);
    //  h2_den_trk_effs[iCent]->RebinY (2);
    //}

    //h2_trk_effs[iCent] = new TEfficiency (*(h2_num_trk_effs[iCent]), *(h2_den_trk_effs[iCent]));
    h2_trk_effs[iCent] = (TH2D*) h2_num_trk_effs[iCent]->Clone (Form ("h2_trk_eff_iCent%i", iCent));
    h2_trk_effs[iCent]->Divide (h2_den_trk_effs[iCent]);

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      //TH1D* num = (TH1D*) ((TEfficiency*) trkEffFile->Get (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta)))->GetCopyPassedHisto ();
      //TH1D* den = (TH1D*) ((TEfficiency*) trkEffFile->Get (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta)))->GetCopyTotalHisto ();

      TH1D* num = (TH1D*) trkEffFile->Get (Form ("h_trk_eff_num_iCent%i_iEta%i", iCent, iEta));
      TH1D* den = (TH1D*) trkEffFile->Get (Form ("h_trk_eff_den_iCent%i_iEta%i", iCent, iEta));

      if (doRebin) {
        RebinSomeBins (&num, maxNPtchBins, allPtchBins);
        RebinSomeBins (&den, maxNPtchBins, allPtchBins);
      }

      //if (iCent > 0) {
      //  num->Rebin (2);
      //  den->Rebin (2);
      //}

      //h_trk_effs[iCent][iEta] =  new TEfficiency (*num, *den);
      h_trk_effs[iCent][iEta] = (TH1D*) num->Clone (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta));
      h_trk_effs[iCent][iEta]->Divide (den);
      //h_trk_effs[iCent][iEta]->SetDirectory (_gDirectory);

      //delete num;
      //delete den;
      //h_trk_effs[iCent][iEta]->SetName ("h_trk_eff_iCent%i_iEta%i");
      
      //h_trk_effs[iCent][iEta] = (TEfficiency*) trkEffFile->Get (Form ("h_trk_eff_iCent%i_iEta%i", iCent, iEta));
    }
  }

  effsLoaded = true;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the appropriate tracking efficiency for this track and centrality.
////////////////////////////////////////////////////////////////////////////////////////////////
double PhysicsAnalysis :: GetTrackingEfficiency (const float fcal_et, double trk_pt, const float trk_eta, const bool isPbPb) {
  if (!effsLoaded)
    LoadTrackingEfficiencies ();

  short iCent = 0;
  if (isPbPb && !useImpactParameter) {
    while (iCent < numTrkCorrCentBins) {
      if (fcal_et < trkCorrCentBins[iCent])
        break;
      else
        iCent++;
    }
    if (iCent == numTrkCorrCentBins) // force ultra-central events to behave like 0-10% Pb+Pb
      iCent--;
    if (iCent < 1 || iCent > numTrkCorrCentBins-1)
      return 0;
  }
  else if (isPbPb && useImpactParameter) {
    iCent = GetIPCentBin (fcal_et); // fcal_et variable should actually be impact parameter if useImpactParameter is true
    if (iCent < 1 || iCent > numTrkCorrCentBins-1)
      return 0;
  }

  TH2D* t = h2_trk_effs[iCent];

  //const int bin = t->FindFixBin (trk_eta, trk_pt);
  //if (bin < 1 || t->GetNbinsX () * t->GetNbinsY () < bin)
  //  return 0;

  const int xbin = t->GetXaxis ()->FindFixBin (trk_eta);
  const int ybin = t->GetYaxis ()->FindFixBin (trk_pt);
  if (xbin < 1 || t->GetXaxis ()->GetNbins () < xbin)
    return 0;
  if (ybin < 1)
    return 0;
  else if (t->GetYaxis ()->GetNbins () < ybin)
    trk_pt = t->GetYaxis ()->GetBinCenter (t->GetYaxis ()->GetNbins ());

  const double eff = t->GetBinContent (xbin, ybin);

  //return t->GetBinContent (t->FindFixBin (trk_eta, trk_pt)) + trkEffNSigma * t->GetBinError (t->FindFixBin (trk_eta, trk_pt));
  return eff + trkEffNSigma * t->GetBinError (xbin, ybin);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load the tracking purities into memory
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LoadTrackingPurities (const bool doRebin) {
  if (pursLoaded)
    return;

  TDirectory* _gDirectory = gDirectory;

  TString _purDir = "Nominal";
  if (useHITight)
    _purDir = "Variations/TrackHITightWPVariation";

  cout << Form ("Reading tracking purities from %s/TrackingPurities/%s/trackingPurities_%s.root", rootPath.Data (), _purDir.Data (), is2015Conds ? "15" : "18") << endl;
  trkPurFile = new TFile (Form ("%s/TrackingPurities/%s/trackingPurities_%s.root", rootPath.Data (), _purDir.Data (), is2015Conds ? "15" : "18"), "read");

  if (!trkPurFile || !trkPurFile->IsOpen ()) {
    cout << "Error in PhysicsAnalysis.h:: LoadTrackingPurities can not find file for " << name << endl;
    return;
  }

  for (int iCent = 0; iCent < numTrkCorrCentBins; iCent++) {
    h2_num_trk_purs[iCent] = (TH2D*) trkPurFile->Get (Form ("h2_primary_reco_tracks_iCent%i", iCent));
    h2_den_trk_purs[iCent] = (TH2D*) trkPurFile->Get (Form ("h2_reco_tracks_iCent%i", iCent));

    h2_num_trk_purs[iCent]->RebinX (4);
    h2_den_trk_purs[iCent]->RebinX (4);

    h2_trk_purs[iCent] = (TH2D*) h2_num_trk_purs[iCent]->Clone (Form ("h2_trk_pur_iCent%i", iCent));
    h2_trk_purs[iCent]->Divide (h2_den_trk_purs[iCent]);

    //for (int iX = 1; iX <= h2_trk_purs[iCent]->GetNbinsX (); iX++) {
    //  for (int iY = 1; iY <= h2_trk_purs[iCent]->GetNbinsY (); iY++) {
    //    const float num = h2_trk_purs[iCent]->GetBinContent (iX, iY);
    //    const float den = h2_den_trk_purs[iCent]->GetBinContent (iX, iY);
    //    if (den != 0 && num != 0)
    //      h2_trk_purs[iCent]->SetBinContent (iX, iY, num / den);
    //  }
    //}

    if (doTrackPurVar) {
      for (int ix = 1; ix <= h2_trk_purs[iCent]->GetNbinsX (); ix++) {
        for (int iy = 1; iy <= h2_trk_purs[iCent]->GetNbinsY (); iy++) {
          float fakeRate = 1 - h2_trk_purs[iCent]->GetBinContent (ix, iy);
          fakeRate = fakeRate + trkPurNSigma * 0.25 * fakeRate;
          h2_trk_purs[iCent]->SetBinContent (ix, iy, 1-fakeRate);
        }
      }
    }

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TH1D* num = (TH1D*) trkPurFile->Get (Form ("h_primary_reco_tracks_iCent%i_iEta%i", iCent, iEta));
      TH1D* den = (TH1D*) trkPurFile->Get (Form ("h_reco_tracks_iCent%i_iEta%i", iCent, iEta));

      if (doRebin) {
        RebinSomeBins (&num, maxNPtchBins, allPtchBins);
        RebinSomeBins (&den, maxNPtchBins, allPtchBins);
      }

      h_trk_purs[iCent][iEta] = (TH1D*) num->Clone (Form ("h_trk_pur_iCent%i_iEta%i", iCent, iEta));
      //h_trk_purs[iCent][iEta]->Divide (den);

      for (int ix = 1; ix <= h_trk_purs[iCent][iEta]->GetNbinsX (); ix++) {
        const float passes = num->GetBinContent (ix);
        const float trials = den->GetBinContent (ix);
        if (passes == 0 || trials == 0)
          continue;
        h_trk_purs[iCent][iEta]->SetBinContent (ix, passes/trials);
        h_trk_purs[iCent][iEta]->SetBinError (ix, sqrt ((passes/trials)*(1-(passes/trials)) * (*den->GetSumw2 ())[ix] / pow (trials, 2)));
        //h_trk_purs[iCent][iEta]->SetBinError (ix, sqrt ((passes+1)*(passes+2) / ((trials+2)*(trials+3)) - pow (passes+1, 2) / pow (trials+2, 2)));
      }

      if (doTrackPurVar) {
        for (int ix = 1; ix <= h_trk_purs[iCent][iEta]->GetNbinsX (); ix++) {
          float fakeRate = 1 - h_trk_purs[iCent][iEta]->GetBinContent (ix);
          fakeRate = fakeRate + trkPurNSigma * 0.25 * fakeRate;
          h_trk_purs[iCent][iEta]->SetBinContent (ix, 1-fakeRate);
        }
      }
    }
  }

  pursLoaded = true;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the appropriate tracking purity factor for this track and centrality.
////////////////////////////////////////////////////////////////////////////////////////////////
double PhysicsAnalysis :: GetTrackingPurity (const float fcal_et, double trk_pt, const float trk_eta, const bool isPbPb) {
  if (!pursLoaded)
    LoadTrackingPurities ();

  short iCent = 0;
  if (isPbPb && !useImpactParameter) {
    while (iCent < numTrkCorrCentBins) {
      if (fcal_et < trkCorrCentBins[iCent])
        break;
      else
        iCent++;
    }
    if (iCent == numTrkCorrCentBins)
      iCent--;
    if (iCent < 1 || iCent > numTrkCorrCentBins-1)
      return 0;
  }
  else if (isPbPb && useImpactParameter) {
    iCent = GetIPCentBin (fcal_et); // fcal_et variable is actually impact parameter if it is being used
    if (iCent < 1 || iCent > numTrkCorrCentBins-1)
      return 0;
  }

  TH2D* t = h2_trk_purs[iCent];

  //const int bin = t->FindFixBin (trk_eta, trk_pt);
  //if (bin < 1 || t->GetNbinsX () * t->GetNbinsY () < bin)
  //  return 0;

  const int xbin = t->GetXaxis ()->FindFixBin (trk_eta);
  const int ybin = t->GetYaxis ()->FindFixBin (trk_pt);
  if (xbin < 1 || t->GetXaxis ()->GetNbins () < xbin)
    return 0;
  if (ybin < 1)
    return 0;
  else if (t->GetYaxis ()->GetNbins () < ybin)
    trk_pt = t->GetYaxis ()->GetBinCenter (t->GetYaxis ()->GetNbins ());

  const double pur = t->GetBinContent (xbin, ybin);
  //const double pur = t->GetBinContent (t->FindFixBin (trk_eta, trk_pt));

  return pur;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Takes references to the q2 vector entries in each detector side and corrects them.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CorrectQ2Vector (float& q2x_a, float& q2y_a, float& q2x_c, float& q2y_c) {
  float temp = eventPlaneCalibrator.CalibrateQ2XA (q2x_a, q2y_a);
  q2y_a = eventPlaneCalibrator.CalibrateQ2YA (q2x_a, q2y_a);
  q2x_a = temp;

  temp = eventPlaneCalibrator.CalibrateQ2XC (q2x_c, q2y_c);
  q2y_c = eventPlaneCalibrator.CalibrateQ2YC (q2x_c, q2y_c);
  q2x_c = temp;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Prints yield of Z's that meet the event selection criteria in each centrality
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PrintZYields (const bool latexForm) {

  if (latexForm) {
    cout << "\t\t\t\\ptz Range & Decay channel & \\pp ";
    for (short iCent = 1; iCent < numCentBins; iCent++)
      cout << Form ("& %i--%i\\%% ", (int)centCuts[iCent], (int)centCuts[iCent-1]);

    cout << "& \\PbPb total \\\\ \\hline \\hline" << endl;
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      cout << "\t\t\t\\multirow{2}{*}";
      if (iPtZ == nPtZBins-1) cout << Form ("{$> \\SI{%g}{\\GeV}$}& \\Zee ", zPtBins[iPtZ]);
      else                    cout << Form ("{%g--$\\SI{%g}{\\GeV}$}& \\Zee ", zPtBins[iPtZ], zPtBins[iPtZ+1]);

      double total = 0;
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        cout << Form ("& %g ", h_z_counts[0][iPtZ][iCent]->GetBinContent (1));

        if (iCent > 0) total += h_z_counts[0][iPtZ][iCent]->GetBinContent (1);
      } // end loop over iCent
      cout << Form ("& %g \\\\", total) << endl;

      cout << "\t\t\t& \\Zmm ";
      total = 0;
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        cout << Form ("& %g ", h_z_counts[1][iPtZ][iCent]->GetBinContent (1));

        if (iCent > 0) total += h_z_counts[1][iPtZ][iCent]->GetBinContent (1);
      } // end loop over iCent
      cout << Form ("& %g \\\\ \\hline", total) << endl;
    } // end loop over iPtZ
  }
  else {
    cout << "ptz Range\tDecay channel\tpp\t";
    for (short iCent = 1; iCent < numCentBins; iCent++)
      cout << Form ("%i--%i%%\t", (int)centCuts[iCent], (int)centCuts[iCent-1]);
    cout << "PbPb total" << endl;

    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (iPtZ == nPtZBins-1) cout << Form ("> %g GeV\tZee\t\t", zPtBins[iPtZ]);
      else                    cout << Form ("%g--%g GeV\tZee\t\t", zPtBins[iPtZ], zPtBins[iPtZ+1]);

      double total = 0;
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        cout << Form ("%g\t", h_z_counts[0][iPtZ][iCent]->GetBinContent (1));

        if (iCent > 0) total += h_z_counts[0][iPtZ][iCent]->GetBinContent (1);
      } // end loop over iCent
      cout << Form ("%g", total) << endl;

      cout << "\t\tZmumu\t\t";
      total = 0;
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        cout << Form ("%g\t", h_z_counts[1][iPtZ][iCent]->GetBinContent (1));

        if (iCent > 0) total += h_z_counts[1][iPtZ][iCent]->GetBinContent (1);
      } // end loop over iCent
      cout << Form ("%g", total) << endl;
    } // end loop over iPtZ
  }
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Prints yield of Z's that meet the event selection criteria in each centrality
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PrintTrackYields (const bool latexForm) {

  if (latexForm) {
    cout << "IMPLEMENT ME!!!" << endl;
  }
  else {
    cout << "ptz Range\t";
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (iPtZ == nPtZBins-1) cout << Form ("> %g GeVt\t", zPtBins[iPtZ]);
      else                    cout << Form ("%g--%g GeV\t", zPtBins[iPtZ], zPtBins[iPtZ+1]);
    }
    cout << endl;

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        double x, avgNtrks, avgNtrksErr;
        g_avg_ntrk_ptz[2][iCent]->GetPoint (iPtZ-2, x, avgNtrks);
        avgNtrksErr = g_avg_ntrk_ptz[2][iCent]->GetErrorY (iPtZ-2);
        if (iCent != 0) 
          cout << Form ("%i--%i%%\t\t", (int)centCuts[iCent], (int)centCuts[iCent-1]);
        else
          cout << "pp\t\t";
        //cout << FormatMeasurement (avgNTrks, avgNTrksErr, 2);
        cout << avgNtrks << "+/-" << avgNtrksErr;
        if (iPtZ == nPtZBins-1) cout << endl;
        else cout << "\t";
      }
    }
  }
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot FCal distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotFCalDists (const bool _treatAsData) {
  const char* canvasName = "c_fcal_et";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 600, 600);
    gDirectory->Add (c);
  }

  c->cd ();

  c->SetLogy ();

  h_fcal_et_reweighted->Scale (1, "width");
  if (!_treatAsData) {
    h_fcal_et_reweighted->SetLineColor (kAzure-1);
  }
  else {
    h_fcal_et_reweighted->SetLineColor (kBlack);
  }

  h_fcal_et_reweighted->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
  h_fcal_et_reweighted->GetYaxis ()->SetTitle ("dN_{evt} / d#Sigma#it{E}_{T} [GeV^{-1}]");

  h_fcal_et_reweighted->GetYaxis ()->SetRangeUser (5e1, 2e7);

  h_fcal_et_reweighted->Draw (canvasExists ? "same hist" : "hist");

  if (!canvasExists) {
    myText (0.22, 0.27, kBlack, "Z-tagged events", 0.04);
    myText (0.22, 0.21, kAzure-1, "Mixed minimum bias", 0.04);

    myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
    myText (0.22, 0.81, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.04);
  }

  c->SaveAs (Form ("%s/FCalDist.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Q2 distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotQ2Dists (const bool _treatAsData) {
  const char* canvasName = "c_q2";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1200, 1200);
    gDirectory->Add (c);
    c->Divide (3, 3);
  }

  c->cd ();

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    c->cd (numCentBins-iCent);
    gPad->SetLogy ();

    double min = 1e30, max = 0;
    GetDrawnObjects ();
    GetMinAndMax (min, max, true);
    SetMinAndMax (min, max);

    TH1D* h = h_q2[iCent];
    if (h->GetMinimum (0) < min)
      min = h->GetMinimum (0);
    if (h->GetMaximum () > max)
      max = h->GetMaximum ();

    h->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    if (!_treatAsData) {
      h->SetLineColor (colors[iCent]);
    }
    else {
      h->SetLineColor (kBlack);
    }

    h->GetXaxis ()->SetTitle ("#it{q}_{2}");
    h->GetYaxis ()->SetTitle ("Counts");

    h->Draw (!canvasExists ? "hist" : "same hist");

    myText (0.51, 0.88, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

    myText (0.51, 0.74, colors[iCent], "Mixed minimum bias", 0.05);
    myText (0.51, 0.80, kBlack, "Z-tagged Data", 0.05);

    SetMinAndMax (0.5*min, 2*max);

  }

  myText (0.52, 0.31, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.52, 0.23, kBlack, "Pb+Pb, 5.02 TeV", 0.05);
  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);

  c->SaveAs (Form ("%s/q2_Mixing/Q2Dists.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Q2 weights
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotQ2Weights (PhysicsAnalysis* a) {
  const char* canvasName = "c_q2_weights";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1200, 1200);
    gDirectory->Add (c);
    c->Divide (2, 2);
  }

  c->cd ();

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    c->cd (numCentBins-iCent);

    TH1D* h = (TH1D*) h_q2_reweighted[iCent]->Clone ();

    const float hint = h->Integral ();
    const float aint = a->h_q2_reweighted[iCent]->Integral ();
    h->Divide (a->h_q2_reweighted[iCent]);

    if (hint != 0)
      h->Scale (aint / hint);

    h->GetYaxis ()->SetRangeUser (0.5, 1.5);

    h->SetMarkerColor (colors[iCent]);
    h->SetLineColor (colors[iCent]);

    h->GetXaxis ()->SetTitle ("#left|#it{q}_{2}#right|");
    h->GetYaxis ()->SetTitle ("Event Weight");

    h->Draw (!canvasExists ? "e1" : "same e1");

    myText (0.61, 0.88, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
  }

  c->cd (1);

  myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.22, 0.84, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);

  c->SaveAs (Form ("%s/q2_Mixing/Q2Weights.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Psi2 distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotPsi2Dists (const bool _treatAsData) {
  const char* canvasName = "c_psi2";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1200, 1200);
    gDirectory->Add (c);
    c->Divide (3, 3);
  }

  c->cd ();

  for (short iCent = 1; iCent < numFineCentBins; iCent++) {
    c->cd (numFineCentBins-iCent);
    gPad->SetLogy ();

    double min = 1e30, max = 0;
    GetDrawnObjects ();
    GetMinAndMax (min, max, true);
    SetMinAndMax (min, max);

    TH1D* h = h_psi2[iCent];
    if (h->GetMinimum (0) < min)
      min = h->GetMinimum (0);
    if (h->GetMaximum () > max)
      max = h->GetMaximum ();

    h->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    if (!_treatAsData) {
      h->SetLineColor (colors[iCent]);
    }
    else {
      h->SetLineColor (kBlack);
    }

    h->GetXaxis ()->SetTitle ("#psi_{2}");
    h->GetYaxis ()->SetTitle ("Counts");

    h->Draw (!canvasExists ? "hist" : "same hist");

    myText (0.51, 0.88, kBlack, Form ("%i-%i%%", (int)fineCentCuts[iCent], (int)fineCentCuts[iCent-1]), 0.06);

    myText (0.51, 0.74, colors[iCent], "Mixed minimum bias", 0.05);
    myText (0.51, 0.80, kBlack, "Z-tagged Data", 0.05);

    SetMinAndMax (0.5*min, 2*max);

  }

  myText (0.52, 0.31, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.52, 0.23, kBlack, "Pb+Pb, 5.02 TeV", 0.05);
  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);

  c->SaveAs (Form ("%s/q2_Mixing/Psi2Dists.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Psi2 weights
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotPsi2Weights (PhysicsAnalysis* a) {
  const char* canvasName = "c_psi2_weights";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1200, 1200);
    gDirectory->Add (c);
    c->Divide (2, 2);
  }

  c->cd ();

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    c->cd (numCentBins-iCent);

    TH1D* h = (TH1D*) h_psi2_reweighted[iCent]->Clone ();

    const float hint = h->Integral ();
    const float aint = a->h_psi2_reweighted[iCent]->Integral ();
    h->Divide (a->h_psi2_reweighted[iCent]);

    if (hint != 0)
      h->Scale (aint / hint);

    h->GetYaxis ()->SetRangeUser (0.5, 1.5);

    h->SetMarkerColor (colors[iCent]);
    h->SetLineColor (colors[iCent]);

    h->GetXaxis ()->SetTitle ("#phi_{Z} - #Psi_{2}");
    h->GetYaxis ()->SetTitle ("Event Weight");

    h->Draw (!canvasExists ? "e1" : "same e1");

    myText (0.61, 0.88, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
  }

  c->cd (1);

  myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.22, 0.84, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
  //myText (0.36, 0.22, kBlack, "Z-tagged data", 0.06);
  //myText (0.36, 0.22, kBlack, "Minimum bias", 0.06);

  c->SaveAs (Form ("%s/q2_Mixing/Psi2Weights.pdf", plotPath.Data ()));

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot VZ distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotVZDists (const bool _treatAsData) {
  const char* canvasName = "c_vz";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1600, 600);
    gDirectory->Add (c);
    c->Divide (2, 1);
    c->cd ();
  }
  c->cd (1);

  gPad->SetLogy ();

  h_PbPb_vz->SetLineColor (kBlack);

  h_PbPb_vz->GetXaxis ()->SetTitle ("#it{v}_{z} [mm]");
  h_PbPb_vz->GetYaxis ()->SetTitle ("Counts");

  h_PbPb_vz->Draw (canvasExists ? "same hist" : "hist");

  c->cd (2);

  gPad->SetLogy ();

  h_pp_vz->SetLineColor (kBlack);

  h_pp_vz->GetXaxis ()->SetTitle ("#it{v}_{z} [mm]");
  h_pp_vz->GetYaxis ()->SetTitle ("Counts");

  h_pp_vz->Draw (canvasExists ? "same hist" : "hist");

  c->cd (1);

  if (!_treatAsData) {
    h_PbPb_vz_reweighted->SetLineColor (kBlue);

    h_PbPb_vz_reweighted->GetXaxis ()->SetTitle ("#it{v}_{z} [mm]");
    h_PbPb_vz_reweighted->GetYaxis ()->SetTitle ("Counts");

    h_PbPb_vz_reweighted->Draw ("same hist");
  }

  myText (0.18, 0.88, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
  //myText (0.18, 0.81, kBlack, "Z-tagged data", 0.04);
  myText (0.18, 0.81, kBlack, "Minimum bias", 0.04);

  c->cd (2);

  if (!_treatAsData) {
    h_pp_vz_reweighted->SetLineColor (kBlue);

    h_pp_vz_reweighted->GetXaxis ()->SetTitle ("#it{v}_{z} [mm]");
    h_pp_vz_reweighted->GetYaxis ()->SetTitle ("Counts");

    h_pp_vz_reweighted->Draw ("same hist");

    myText (0.75, 0.88, kBlack, "Unweighted", 0.04);
    myText (0.75, 0.81, kBlue, "Reweighted", 0.04);
  }

  myText (0.18, 0.88, kBlack, "#it{pp}, 5.02 TeV", 0.04);
  //myText (0.18, 0.81, kBlack, "Z-tagged data", 0.04);
  myText (0.18, 0.81, kBlack, "Minimum bias", 0.04);

  c->SaveAs (Form ("%s/VZDist.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot VZ distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotNchDists (const bool _treatAsData) {
  const char* canvasName = "c_nch";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
  }

  c->cd ();

  gPad->SetLogy ();

  if (_treatAsData)
    h_pp_nch->SetLineColor (kBlack);
  else
    h_pp_nch->SetLineColor (kGray+1);

  h_pp_nch->Scale (1/h_pp_nch->Integral ());
  h_pp_nch->GetXaxis ()->SetTitle ("N_{ch}");
  h_pp_nch->GetYaxis ()->SetTitle ("dN_{evt} / N_{evt}");

  h_pp_nch->Draw (canvasExists ? "same hist" : "hist");

  if (!_treatAsData) {
    h_pp_nch_reweighted->Scale (1/h_pp_nch_reweighted->Integral ());
    h_pp_nch_reweighted->SetLineColor (kRed+1);

    h_pp_nch_reweighted->GetXaxis ()->SetTitle ("N_{ch}");
    h_pp_nch_reweighted->GetYaxis ()->SetTitle ("dN_{evt} / N_{evt");

    h_pp_nch_reweighted->Draw ("same hist");
  }

  myText (0.68, 0.88, kBlack, "#it{pp}, 5.02 TeV", 0.04);
  myText (0.68, 0.81, kBlack, "Z-tagged data", 0.04);
  myText (0.68, 0.74, kGray+1, "Minimum bias", 0.04);
  myText (0.68, 0.67, kRed+1, "MB Reweighted", 0.04);

  c->SaveAs (Form ("%s/NchDist.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dPhi - pTch 2d projections
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotCorrelations (const short pSpc, const short pPtZ, const bool _subBkg) {

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue;

      const char* canvasName = Form ("c_deltaPhi_pTch_iPtZ%i_%s", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 600, 200*numCentBins);
        gDirectory->Add (c);
        c->cd ();
        c->Divide (2, numCentBins);
      }

      c->cd ();

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (bool subBkg : {false, true}) {

          c->cd ((2*iCent)+1 + (int)subBkg);

          if (!canvasExists) {
            gPad->SetTopMargin (0.024);
            gPad->SetBottomMargin (0.12);
            gPad->SetRightMargin (0.01);
            gPad->SetLeftMargin (0.12);

            GetDrawnObjects ();

            double min = -5, max = 10;
            if (!subBkg) {
              GetMinAndMax (min, max, true);
              for (short iPtch = 0; iPtch < std::min (3, nPtchBins[iPtZ]); iPtch++) {
                TH1D* h = (!subBkg ? h_trk_dphi[iSpc][iPtZ][iPtch][iCent] : h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]);
                min = fmin (min, h->GetMinimum ());
                max = fmax (max, h->GetMaximum ());
              } // end loop over iPtch
              if (max != 2.7) max = (max <= 0 ? 1 : 1.2*max);
              SetMinAndMax (min, max);
            }

            TH1D* htemp = new TH1D (Form ("_hplot_iCent%i_subBkg%i", iCent, subBkg ? 1 : 0), "", 1, -pi/2, 3*pi/2);

            TAxis* xax = htemp->GetXaxis ();
            TAxis* yax = htemp->GetYaxis ();

            yax->SetRangeUser (min, max);

            htemp->SetLineWidth (0);
            htemp->SetMarkerSize (0);

            xax->SetTitle ("#Delta#phi");
            yax->SetTitle ("dY / d#Delta#phi");

            xax->SetTitleOffset (0.6);
            yax->SetTitleOffset (0.8);
            xax->SetTitleSize (0.08);
            yax->SetTitleSize (0.07);
            xax->SetLabelSize (0.06);
            yax->SetLabelSize (0.06);

            htemp->SetLineWidth (0);
            htemp->SetMarkerSize (0);
  
            htemp->DrawCopy ("hist ][");
            SaferDelete (&htemp);

            TLine* line1 = new TLine (3*pi/4, min, 3*pi/4, (iCent == 0 ? 0.8:1.)*max);
            TLine* line2 = new TLine (5*pi/4, min, 5*pi/4, (iCent == 0 ? 0.8:1.)*max);

            line1->SetLineStyle (2);
            line2->SetLineStyle (2);
            line1->SetLineWidth (2);
            line2->SetLineWidth (2);
            line1->SetLineColor (kBlack);
            line2->SetLineColor (kBlack);

            line1->Draw ("same");
            line2->Draw ("same");

            if (subBkg) {
              TLine* line3 = new TLine (-pi/2, 0, 3*pi/2, 0);
              line3->SetLineStyle (2);
              line3->SetLineWidth (2);
              line3->SetLineColor (kBlack);
              line3->Draw ("same");
            }
 
            for (short iPtch = 0; iPtch < std::min (3, nPtchBins[iPtZ]); iPtch++)
              LabelCorrelations (iPtZ, iPtch, iCent, subBkg);
          }

          if (subBkg && !_subBkg) continue;

          if (plotFill) {
            for (short iPtch = 0; iPtch < std::min (3, nPtchBins[iPtZ]); iPtch++) {
              TH1D* h = (TH1D*) (!subBkg ? h_trk_dphi[iSpc][iPtZ][iPtch][iCent] : h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent])->Clone ();

              h->SetLineColorAlpha (fillColors[iPtch], 0.8);
              h->SetLineWidth (10);
              h->SetMarkerSize (0);

              h->DrawCopy ("same hist");
            } // end loop over iPtch
            gPad->RedrawAxis ();
          } else {
            for (short iPtch = 0; iPtch < std::min (3, nPtchBins[iPtZ]); iPtch++) {
              TGAE* g = GetTGAE (!subBkg ? h_trk_dphi[iSpc][iPtZ][iPtch][iCent] : h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent]);
              ResetXErrors (g);

              const Style_t markerStyle = markerStyles[iPtch];

              g->SetMarkerStyle (markerStyle);
              g->SetLineColor (colors[iPtch]);
              g->SetMarkerColor (colors[iPtch]);

              g->Draw ("P");
            } // end loop over iPtch
          }
        }
      } // end loop over iCent

      c->SaveAs (Form ("%s/DeltaPhiCorrelations/deltaPhi_pTch_iPtZ%i_%s.pdf", plotPath.Data (), iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for track dPhi distributions
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelCorrelations (const short iPtZ, const short iPtch, const short iCent, const bool subBkg) {
  if (iCent == 0) {
    if (iPtch == 0) {
      if (!subBkg) {
        //myText (0.2, 0.91, kBlack, "#bf{#it{ATLAS}} Simulation", 0.07);
        myText (0.2, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.077);
        myText (0.2, 0.79, kBlack, "#it{pp}, 5.02 TeV", 0.07);
        if (iPtZ == nPtZBins-1) myText (0.2, 0.70, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.07);
        else                    myText (0.2, 0.70, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.07);
        //myMarkerTextNoLine (0.24, 0.64, kBlack, kFullCircle, "Z-tagged Data", 1.25, 0.07);
        //myOnlyBoxText (0.24, 0.55, 1.2, fillColors[0], kBlack, 1, "Minimum Bias", 0.07, 1001, 1);
        myText (0.63, 0.90, kBlack, "Before subtraction", 0.07);
      }
      else {
        myText (0.17, 0.88, kBlack, "UE", 0.07);
        myText (0.26, 0.88, kBlack, "#it{Z}-tagged", 0.07);
        myText (0.65, 0.90, kBlack, "After subtraction", 0.07);
      }
    }
    if (subBkg) {
      const float pt_lo = pTchBins[iPtZ][iPtch];
      const float pt_hi = pTchBins[iPtZ][iPtch+1];
      //if (iPtch == 0)
      //  myText (0.3, 0.93, kBlack, "#it{p}_{T}^{ ch} [GeV]", 0.04);
      myOnlyBoxText      (0.25, 0.81-0.075*(iPtch), 1.2, fillColors[iPtch], kBlack, 1, "", 0.07, 1001, 0.8);
      myMarkerTextNoLine (0.37, 0.81-0.075*(iPtch), colors[iPtch], markerStyles[iPtch], "", 1.4, 0.06);
      myText             (0.40, 0.79-0.075*(iPtch),  kBlack, Form ("%.0f < #it{p}_{T}^{ch} < %.0f GeV", pt_lo, pt_hi), 0.06);
    }
  }
  else if (!subBkg && iPtch == 0 && iCent != 0) {
    myText (0.2, 0.90, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.07);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingEfficiencies () {
  if (!effsLoaded) LoadTrackingEfficiencies ();

  const char* canvasName = "c_trk_effs";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 400);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 1);
  }
  c->cd ();

  //TGraphErrors* g0[numCentBins];
  //TGraphErrors* g1[numCentBins];
  //TGraphErrors* g2[numCentBins];
  //TGraphErrors* g3[numCentBins];

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    gPad->SetLogx ();
    gPad->SetLeftMargin (0.18);
    gPad->SetRightMargin (0.01);

    //g0[iCent] = new TGraphErrors (numEtaTrkBins);
    //g1[iCent] = new TGraphErrors (numEtaTrkBins);
    //g2[iCent] = new TGraphErrors (numEtaTrkBins);
    //g3[iCent] = new TGraphErrors (numEtaTrkBins);

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      //TEfficiency* eff = h_trk_effs[iCent][iEta];
      TGAE* eff = GetTGAE (h_trk_effs[iCent][iEta]);

      eff->SetLineColor (colors[iEta]);
      eff->SetMarkerColor (colors[iEta]);
      eff->SetMarkerStyle (useAltMarker ? kOpenCircle : kFullCircle);
      eff->SetMarkerSize (0.5);

      eff->SetTitle (";#it{p}_{T} [GeV];Weighted Reco. Eff.");
      eff->GetXaxis ()->SetRangeUser (0.5, 60);
      eff->GetYaxis ()->SetRangeUser (0.3, 1.12);

      eff->GetXaxis ()->SetTitleSize (0.07);
      eff->GetYaxis ()->SetTitleSize (0.07);
      eff->GetXaxis ()->SetTitleOffset (0.7 * eff->GetXaxis ()->GetTitleOffset ());
      eff->GetYaxis ()->SetTitleOffset (0.9 * eff->GetYaxis ()->GetTitleOffset ());

      eff->GetXaxis ()->SetMoreLogLabels ();

      eff->Draw (!canvasExists && iEta == 0 ? "AP" : "P");
      //eff->Draw (iEta == 0 ? "APL" : "LP same");

      //gPad->Update ();

      //eff->GetPaintedGraph ()->GetXaxis ()->SetRangeUser (0.5, 65);
      //eff->GetPaintedGraph ()->GetYaxis ()->SetRangeUser (0.3, 1.08);

      ////TH1D* confInts = (TH1D*) h_trk_effs[iCent][iEta]->Clone ("confInts");
      //TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2 + [3]*(log(x))^3", 1, 65);
      ////TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2", 2, 65);
      //fit->SetParameter (0, 1);
      //fit->SetParameter (1, 0);
      //fit->SetParameter (2, 0);
      ////fit->SetParameter (3, 0);
      //h_trk_effs[iCent][iEta]->Fit (fit, "RN0Q");

      //g0[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (0));
      //g1[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (1));
      //g2[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (2));
      //g3[iCent]->SetPoint (iEta, 0.5*(etaTrkBins[iEta]+etaTrkBins[iEta+1]), fit->GetParameter (3));

      //g0[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (0));
      //g1[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (1));
      //g2[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (2));
      //g3[iCent]->SetPointError (iEta, 0.5*(etaTrkBins[iEta+1]-etaTrkBins[iEta]), fit->GetParError (3));

      //fit->SetLineColor (colors[iEta]);
      //fit->SetFillColor (colors[iEta]);
      //fit->SetLineWidth (1);
      //fit->Draw ("same");
      //(TVirtualFitter::GetFitter ())->GetConfidenceIntervals (confInts, 0.68);
      //confInts->SetMarkerSize (0);
      //confInts->SetFillColorAlpha (colors[iEta], 0.3);
      //confInts->SetLineColor (colors[iEta]);
      //confInts->DrawCopy ("e3 same");
      //delete confInts;

      if (iCent == 0)
        myMarkerTextNoLine (0.5, 0.50-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
      if (iEta == 0) {
        if (iCent == 0) myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
        else            myText (0.22, 0.86, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
      }
    }
    TLine* l = new TLine (0.5, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kPink-8);
    l->Draw ("same");
  }

  c->SaveAs (Form ("%s/TrackingEfficiencies/TrackingEfficiencies.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies with comparison to some alternative
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingEfficienciesComparison (PhysicsAnalysis* a) {
  if (!effsLoaded) LoadTrackingEfficiencies ();

  const char* canvasName = "c_trk_effs_comp";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 800);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 2);
  }
  c->cd ();

  //TGraphErrors* g0[numCentBins];
  //TGraphErrors* g1[numCentBins];
  //TGraphErrors* g2[numCentBins];
  //TGraphErrors* g3[numCentBins];

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    gPad->SetLogx ();
    gPad->SetLeftMargin (0.18);
    gPad->SetRightMargin (0.01);

    //g0[iCent] = new TGraphErrors (numEtaTrkBins);
    //g1[iCent] = new TGraphErrors (numEtaTrkBins);
    //g2[iCent] = new TGraphErrors (numEtaTrkBins);
    //g3[iCent] = new TGraphErrors (numEtaTrkBins);

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      //TEfficiency* eff = h_trk_effs[iCent][iEta];
      TGAE* eff = GetTGAE (h_trk_effs[iCent][iEta]);

      eff->SetLineColor (colors[iEta]);
      eff->SetMarkerColor (colors[iEta]);
      eff->SetMarkerStyle (useAltMarker ? kOpenCircle : kFullCircle);
      eff->SetMarkerSize (0.5);

      eff->SetTitle (";#it{p}_{T} [GeV];Weighted Reco. Eff.");
      eff->GetXaxis ()->SetRangeUser (0.5, 60);
      eff->GetYaxis ()->SetRangeUser (0.3, 1.12);

      eff->GetXaxis ()->SetTitleSize (0.07);
      eff->GetYaxis ()->SetTitleSize (0.07);
      eff->GetXaxis ()->SetTitleOffset (0.7 * eff->GetXaxis ()->GetTitleOffset ());
      eff->GetYaxis ()->SetTitleOffset (0.9 * eff->GetYaxis ()->GetTitleOffset ());

      eff->GetXaxis ()->SetMoreLogLabels ();

      eff->Draw (!canvasExists && iEta == 0 ? "AP" : "P");

      if (iEta == 0) {
        if (iCent == 0) myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
        else            myText (0.22, 0.86, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
        if (iCent == 1) {

          if (a && a->doTrackEffVar) {
            myMarkerTextNoLine (0.36, 0.16, kBlack, kFullCircle, "Inclusive hadrons", 1.2, 0.06);
            myMarkerTextNoLine (0.36, 0.10, kBlack, kOpenCircle, "Pions only", 1.2, 0.06);
          }
          else if (a && a->useHITight) {
            myMarkerTextNoLine (0.36, 0.16, kBlack, kFullCircle, "HILoose tracks", 1.2, 0.06);
            myMarkerTextNoLine (0.36, 0.10, kBlack, kOpenCircle, "HITight tracks", 1.2, 0.06);
          }
        }
      }
      if (iCent == 0) myMarkerTextNoLine (0.5, 0.34-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
    }
    TLine* l = new TLine (0.5, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kPink-8);
    l->Draw ("same");

    if (!a) continue;
    else    a->LoadTrackingEfficiencies ();

    gPad->SetBottomMargin (0.02);
    c->cd (iCent+5);
    gPad->SetLogx ();
    gPad->SetTopMargin (0.02);
    gPad->SetLeftMargin (0.18);
    gPad->SetRightMargin (0.01);
    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TH1D* h = (TH1D*) a->h_trk_effs[iCent][iEta]->Clone ("temp");
      h->Divide (h_trk_effs[iCent][iEta]);

      TGAE* eff = GetTGAE (h);
      SaferDelete (&h);

      eff->SetLineColor (colors[iEta]);
      eff->SetMarkerColor (colors[iEta]);
      eff->SetMarkerStyle (useAltMarker ? kOpenCircle : kFullCircle);
      eff->SetMarkerSize (0.5);

      if (a && a->doTrackEffVar) {
        eff->SetTitle (";#it{p}_{T} [GeV];Pions / Inclusive hadrons");
        eff->GetYaxis ()->SetRangeUser (0.89, 1.11);
      }
      else if (a && a->useHITight) {
        eff->SetTitle (";#it{p}_{T} [GeV];HITight / HILoose");
        eff->GetYaxis ()->SetRangeUser (0.7, 1.11);
      }
      eff->GetXaxis ()->SetRangeUser (0.5, 60);

      eff->GetXaxis ()->SetTitleSize (0.07);
      eff->GetYaxis ()->SetTitleSize (0.07);
      eff->GetXaxis ()->SetTitleOffset (0.7 * eff->GetXaxis ()->GetTitleOffset ());
      eff->GetYaxis ()->SetTitleOffset (0.9 * eff->GetYaxis ()->GetTitleOffset ());

      eff->GetYaxis ()->CenterTitle ();

      eff->GetXaxis ()->SetMoreLogLabels ();

      eff->Draw (!canvasExists && iEta == 0 ? "AP" : "P");
    }

    l->Draw ("same");
  }

  if (a) {
    bool temp = a->useAltMarker;
    a->useAltMarker = true;
    a->PlotTrackingEfficienciesComparison ();
    a->useAltMarker = temp;

    if (a->doTrackEffVar) c->SaveAs (Form ("%s/TrackingEfficiencies/TrackingEfficienciesParticleComp.pdf", plotPath.Data ()));
    else if (a->useHITight) c->SaveAs (Form ("%s/TrackingEfficiencies/TrackingEfficienciesIDComp.pdf", plotPath.Data ()));
    else c->SaveAs (Form ("%s/TrackingEfficiencies/TrackingEfficienciesComparison.pdf", plotPath.Data ()));
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking efficiencies as a 2D histogram
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingEfficiencies2D () {
  if (!effsLoaded) LoadTrackingEfficiencies ();

  const char* canvasName = "c_trk_effs_2d";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 400);
    FormatTH2Canvas (c, true);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 1);
  }
  c->cd ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    //if (iCent > 0) continue;
    gPad->SetLogy ();

    //TEfficiency* eff = h2_trk_effs[iCent];
    TH2D* eff = h2_trk_effs[iCent];
    eff->GetZaxis ()->SetRangeUser (0.3, 1.00);
    eff->GetZaxis ()->SetTitle ("Weighted Reco. Eff.");
    eff->GetZaxis ()->SetTitleOffset (1.4 * eff->GetZaxis ()->GetTitleOffset ());

    eff->GetYaxis ()->SetMoreLogLabels ();

    eff->Draw ("colz");
    //eff->GetPaintedHistogram ()->GetYaxis ()->SetRangeUser (2, 10);
    gPad->Update ();

    if (iCent == 0) myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
    else            myText (0.22, 0.86, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
  }

  c->SaveAs (Form ("%s/TrackingEfficiencies/TrackingEfficiencies2D.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking purities
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingPurities () {
  if (!pursLoaded) LoadTrackingPurities ();

  const char* canvasName = "c_trk_purs";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 400);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 1);
  }
  c->cd ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    gPad->SetLogx ();
    gPad->SetLeftMargin (0.18);
    gPad->SetRightMargin (0.01);

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TGAE* pur = GetTGAE (h_trk_purs[iCent][iEta]);

      pur->SetLineColor (colors[iEta]);
      pur->SetMarkerColor (colors[iEta]);
      pur->SetMarkerStyle (useAltMarker ? kOpenCircle : kFullCircle);
      pur->SetMarkerSize (0.5);

      pur->SetTitle (";#it{p}_{T} [GeV];Primary Track Fraction");
      pur->GetXaxis ()->SetRangeUser (0.5, 60);
      pur->GetYaxis ()->SetRangeUser (0.90, 1.030);

      pur->GetXaxis ()->SetTitleSize (0.07);
      pur->GetYaxis ()->SetTitleSize (0.07);
      pur->GetXaxis ()->SetTitleOffset (0.7 * pur->GetXaxis ()->GetTitleOffset ());
      pur->GetYaxis ()->SetTitleOffset (0.9 * pur->GetYaxis ()->GetTitleOffset ());

      pur->GetXaxis ()->SetMoreLogLabels ();

      pur->Draw (!canvasExists && iEta == 0 ? "AP" : "P");

      //TH1D* confInts = (TH1D*) h_trk_purs[iCent][iEta]->Clone ("confInts");
      //TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2 + [3]*(log(x))^3 + [4]*(log(x))^4", 1, 65);
      //TF1* fit = new TF1 ("fit", "[0] + [1]*log(x) + [2]*(log(x))^2", 2, 65);
      //fit->SetParameter (0, 1);
      //fit->SetParameter (1, 0);
      //fit->SetParameter (2, 0);
      //fit->SetParameter (3, 0);
      //fit->FixParameter (4, 0);
      //h_trk_purs[iCent][iEta]->Fit (fit, "RN0Q");
      //fit->SetLineColor (colors[iEta]);
      //fit->SetFillColor (colors[iEta]);
      //fit->SetLineWidth (1);
      //fit->Draw ("same");
      //(TVirtualFitter::GetFitter ())->GetConfidenceIntervals (confInts, 0.68);
      //confInts->SetMarkerSize (0);
      //confInts->SetFillColorAlpha (colors[iEta], 0.3);
      //confInts->SetLineColor (colors[iEta]);
      //confInts->DrawCopy ("e3 same");
      //delete confInts;

      if (iCent == 0)
        myMarkerTextNoLine (0.3, 0.50-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
      if (iEta == 0) {
        if (iCent == 0) myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
        else            myText (0.22, 0.86, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
      }
    }
    TLine* l = new TLine (0.5, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kPink-8);
    l->Draw ("same");
  }

  c->SaveAs (Form ("%s/TrackingPurities/TrackingPurities.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking purities with comparison to some alternative
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingPuritiesComparison (PhysicsAnalysis* a) {
  if (!pursLoaded) LoadTrackingPurities ();

  const char* canvasName = "c_trk_purs_comp";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 800);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 2);
  }
  c->cd ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    gPad->SetLogx ();
    gPad->SetLeftMargin (0.18);
    gPad->SetRightMargin (0.01);

    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TGAE* pur = GetTGAE (h_trk_purs[iCent][iEta]);

      pur->SetLineColor (colors[iEta]);
      pur->SetMarkerColor (colors[iEta]);
      pur->SetMarkerStyle (useAltMarker ? kOpenCircle : kFullCircle);
      pur->SetMarkerSize (0.5);

      pur->SetTitle (";#it{p}_{T} [GeV];Primary Track Fraction");
      pur->GetXaxis ()->SetRangeUser (0.5, 60);
      pur->GetYaxis ()->SetRangeUser (0.90, 1.030);

      pur->GetXaxis ()->SetTitleSize (0.07);
      pur->GetYaxis ()->SetTitleSize (0.07);
      pur->GetXaxis ()->SetTitleOffset (0.7 * pur->GetXaxis ()->GetTitleOffset ());
      pur->GetYaxis ()->SetTitleOffset (0.9 * pur->GetYaxis ()->GetTitleOffset ());

      pur->GetXaxis ()->SetMoreLogLabels ();

      pur->Draw (!canvasExists && iEta == 0 ? "AP" : "P");

      if (iEta == 0) {
        if (iCent == 0) myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
        else            myText (0.22, 0.86, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
        if (iCent == 1) {
          myMarkerTextNoLine (0.36, 0.16, kBlack, kFullCircle, "HILoose tracks", 1.2, 0.06);
          myMarkerTextNoLine (0.36, 0.10, kBlack, kOpenCircle, "HITight tracks", 1.2, 0.06);
        }
      }
      if (iCent == 0) myMarkerTextNoLine (0.5, 0.34-0.06*iEta, colors[iEta], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.06);
    }
    TLine* l = new TLine (0.5, 1, 60, 1);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kPink-8);
    l->Draw ("same");

    if (!a) continue;
    else    a->LoadTrackingPurities ();

    gPad->SetBottomMargin (0.02);
    c->cd (iCent+5);
    gPad->SetLogx ();
    gPad->SetTopMargin (0.02);
    gPad->SetLeftMargin (0.18);
    gPad->SetRightMargin (0.01);
    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      TH1D* h = (TH1D*) a->h_trk_purs[iCent][iEta]->Clone ("temp");
      h->Divide (h_trk_purs[iCent][iEta]);

      TGAE* pur = GetTGAE (h);
      SaferDelete (&h);

      pur->SetLineColor (colors[iEta]);
      pur->SetMarkerColor (colors[iEta]);
      pur->SetMarkerStyle (useAltMarker ? kOpenCircle : kFullCircle);
      pur->SetMarkerSize (0.5);

      pur->SetTitle (";#it{p}_{T} [GeV];HITight / HILoose");
      pur->GetXaxis ()->SetRangeUser (0.5, 60);
      pur->GetYaxis ()->SetRangeUser (0.94, 1.06);

      pur->GetXaxis ()->SetTitleSize (0.07);
      pur->GetYaxis ()->SetTitleSize (0.07);
      pur->GetXaxis ()->SetTitleOffset (0.7 * pur->GetXaxis ()->GetTitleOffset ());
      pur->GetYaxis ()->SetTitleOffset (0.9 * pur->GetYaxis ()->GetTitleOffset ());

      pur->GetYaxis ()->CenterTitle ();

      pur->GetXaxis ()->SetMoreLogLabels ();

      pur->Draw (!canvasExists && iEta == 0 ? "AP" : "P");
    }

    l->Draw ("same");
  }

  if (a) {
    bool temp = a->useAltMarker;
    a->useAltMarker = true;
    a->PlotTrackingPuritiesComparison ();
    a->useAltMarker = temp;

    if (a->useHITight) c->SaveAs (Form ("%s/TrackingPurities/TrackingPuritiesIDComp.pdf", plotPath.Data ()));
    else c->SaveAs (Form ("%s/TrackingPurities/TrackingPuritiesComparison.pdf", plotPath.Data ()));
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots tracking purities as a 2D histogram
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackingPurities2D () {
  if (!pursLoaded)
    LoadTrackingPurities ();

  const char* canvasName = "c_trk_purs_2d";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1800, 400);
    FormatTH2Canvas (c, true);
    gDirectory->Add (c);
    c->cd ();
    c->Divide (4, 1);
  }
  c->cd ();

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd (iCent+1);
    //if (iCent > 0) continue;
    gPad->SetLogy ();

    //TEfficiency* pur = h2_trk_purs[iCent];
    TH2D* pur = h2_trk_purs[iCent];
    pur->GetZaxis ()->SetRangeUser (0.90, 1.00);
    pur->GetZaxis ()->SetTitle ("Primary Track Fraction");
    pur->GetZaxis ()->SetTitleOffset (1.4 * pur->GetZaxis ()->GetTitleOffset ());

    pur->GetYaxis ()->SetMoreLogLabels ();

    pur->Draw ("colz");
    //pur->GetPaintedHistogram ()->GetYaxis ()->SetRangeUser (2, 10);
    gPad->Update ();

    if (iCent == 0) myText (0.22, 0.86, kBlack, "#it{pp}", 0.072);
    else            myText (0.22, 0.86, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.072);
  }

  c->SaveAs (Form ("%s/TrackingPurities/TrackingPurities2D.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot pTch distributions binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotUnweightedTrkYields (const bool useTrkPt, const short pSpc, const short pPtZ) {
  const int axisTextSize = 23;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_UnweightedTrkYields_%s_dPhi_iPtZ%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 400*numCentBins, 400);
        gDirectory->Add (c);
      }
      c->cd ();

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        c->cd ();

        const char* padName = Form ("p_UnweightedTrkYields_%s_iPtZ%i_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, iCent, spc);

        TPad* pad = nullptr;
        if (!canvasExists) {
          pad = new TPad (padName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, 1);

          gDirectory->Add (pad);

          pad->SetTopMargin (0.04);
          pad->SetBottomMargin (0.20);
          pad->SetLeftMargin (0.17);
          pad->SetRightMargin (0.06);
          pad->Draw ();
        }
        else {
          pad = dynamic_cast<TPad*> (gDirectory->Get (padName));
        }

        pad->cd ();
        GetDrawnObjects ();
        bool plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
        gPad->SetLogx ();
        gPad->SetLogy ();

        double min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          TH1D* h = (useTrkPt ? h_trk_pt_dphi_raw : h_trk_xhz_dphi)[iSpc][iPtZ][iPhi][iCent];
          min = fmin (min, h->GetMinimum (0));
          max = fmax (max, h->GetMaximum ());
        } // end loop over iPhi
        min = (min > 0 ? (canvasExists ? 0.5 : 1)*min : 0.1);
        max = (max > 0 ? (canvasExists ? 2 : 1)*max : 1);
        SetMinAndMax (min, max);

        if (plotFill) {
          for (int iPhi = 0; iPhi < 1; iPhi++) {
          //for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
            TH1D* h = (useTrkPt ? h_trk_pt_dphi_raw : h_trk_xhz_dphi)[iSpc][iPtZ][iPhi][iCent];

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetMarkerSize (0);
            h->SetLineColor (kBlack);
            h->SetLineWidth (0);

            useTrkPt ? h->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : h->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
            h->GetYaxis ()->SetRangeUser (min, max);

            h->GetXaxis ()->SetMoreLogLabels ();

            h->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
            h->GetYaxis ()->SetTitle (useTrkPt ? "Y (#it{p}_{T}, #Delta#phi)" : "Y (#it{x}_{hZ}, #Delta#phi)");

            h->GetXaxis ()->SetTitleFont (43);
            h->GetXaxis ()->SetTitleSize (axisTextSize);
            h->GetXaxis ()->SetLabelFont (43);
            h->GetXaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleFont (43);
            h->GetYaxis ()->SetTitleSize (axisTextSize);
            h->GetYaxis ()->SetLabelFont (43);
            h->GetYaxis ()->SetLabelSize (axisTextSize);

            h->GetYaxis ()->SetTitleOffset (1.1);

            h->DrawCopy (plotNewAxes && iPhi == 0 ? "bar" : "bar same");
            //h->DrawCopy (plotNewAxes && iPhi == numPhiBins-1 ? "bar" : "bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");
          }
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            const Style_t markerStyle = (useAltMarker ? (iPhi == 0 ? kOpenSquare : kOpenCircle) : (iPhi == 0 ? kFullSquare : kFullCircle));
            TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi_raw : h_trk_xhz_dphi)[iSpc][iPtZ][iPhi][iCent]);
            RecenterGraph (g);

            if (!plotAsSyst) {
              ResetXErrors (g);
              g->SetMarkerStyle (markerStyle);
              g->SetMarkerColor (colors[iPhi]);
              g->SetLineColor (colors[iPhi]);
              g->SetMarkerSize (1);
              g->SetLineWidth (2);
            } else {
              g->SetMarkerSize (0); 
              g->SetLineWidth (1);
              g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);
            }

            useTrkPt ? g->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
            g->GetYaxis ()->SetRangeUser (min, max);

            g->GetXaxis ()->SetMoreLogLabels ();

            g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
            g->GetYaxis ()->SetTitle (useTrkPt ? "Y (#it{p}_{T}, #Delta#phi)" : "Y (#it{x}_{hZ}, #Delta#phi)");

            g->GetXaxis ()->SetTitleFont (43);
            g->GetXaxis ()->SetTitleSize (axisTextSize);
            g->GetXaxis ()->SetLabelFont (43);
            g->GetXaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleFont (43);
            g->GetYaxis ()->SetTitleSize (axisTextSize);
            g->GetYaxis ()->SetLabelFont (43);
            g->GetYaxis ()->SetLabelSize (axisTextSize);

            g->GetYaxis ()->SetTitleOffset (1.1);

            if (!plotAsSyst) {
              string drawString = string (!canvasExists && iPhi == 0 ? "AP" : "P");
              g->Draw (drawString.c_str ());
            } else {
              string drawString = string (!canvasExists && iPhi == 0 ? "A5P" : "5P");
              ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
              g->Draw ("2P");
            }
          } // end loop over iPhi
        }

        if (!canvasExists) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            if (iCent == 0)
              myText (0.22, 0.26, kBlack, "#it{pp}, 5.02 TeV", 0.06);
            else
              myText (0.22, 0.26, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

            if (iCent == 0)
              myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.068);
            else if (iCent == numCentBins-1) {
              if (iPhi == 1) {
                TVirtualPad* cPad = gPad; // store current pad
                myText (0.653, 0.91, kBlack, "#Delta#phi", 0.054);
                TBox* b = TBoxNDC (0.598-0.018, 0.91-0.06*numPhiBins-0.016, 0.598+0.018, 0.91-0.06*numPhiBins+0.016);
                b->SetFillColorAlpha (fillColors[0], fillAlpha);
                b->Draw ("l");
                cPad->cd ();

                myText (0.65, 0.91-0.06*numPhiBins, kBlack, "MinBias", 0.054);
              }

              const char* lo = GetPiString (phiLowBins[iPhi]);
              const char* hi = GetPiString (phiHighBins[iPhi]);

              myMarkerTextNoLine (0.62, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.5, 0.054);
              myText (0.65, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.054);
            }
          }
        }
      } // end loop over iCent
      
      c->SaveAs (Form ("%s/UnweightedTrkYields/allYields_pTch_iPtZ%i_%s.pdf", plotPath.Data (), iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotAllYields_dPhi (const bool useTrkPt, const short pSpc, const short pPtZ) {
  const double padRatio = 0.9; // ratio of size of upper pad to middle & lower pads. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (padRatio+1.0);
  const int axisTextSize = 23;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_AllYields_%s_iPtZ%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 395*numCentBins, 685);
        gDirectory->Add (c);
      }

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        c->cd ();

        const char* topPadName = Form ("p_AllYields_top_%s_iPtZ%i_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, iCent, spc);
        const char* bottomPadName = Form ("p_AllYields_bottom_%s_iPtZ%i_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, iCent, spc);

        TPad* topPad = nullptr, *bottomPad = nullptr;
        if (!canvasExists) {
          topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
          bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

          gDirectory->Add (topPad);
          gDirectory->Add (bottomPad);

          topPad->SetTopMargin (0.04);
          topPad->SetBottomMargin (0.03);
          topPad->SetLeftMargin (0.23);
          topPad->SetRightMargin (0.03);
          bottomPad->SetTopMargin (0.02);
          bottomPad->SetBottomMargin (0.23);
          bottomPad->SetLeftMargin (0.23);
          bottomPad->SetRightMargin (0.03);
          topPad->Draw ();
          bottomPad->Draw ();
        }
        else {
          topPad = dynamic_cast<TPad*> (gDirectory->Get (topPadName));
          bottomPad = dynamic_cast<TPad*> (gDirectory->Get (bottomPadName));
        }

        topPad->cd ();
        gPad->SetLogx ();
        gPad->SetLogy ();

        double min = 1e30, max = 0;
        if (useTrkPt) {
          if (iCent == 0) { min = 1e-7; max = 2e1; }
          else if (iCent == 1) { min = 1e-6; max = 2e2; }
          else if (iCent == 2) { min = 5e-5; max = 1e3; }
          else if (iCent == 3) { min = 1e-4; max = 2e3; }
        }
        else {
          if (iCent == 0) { min = 2e-5; max = 2e3; }
          else if (iCent == 1) { min = 1e-3; max = 2e4; }
          else if (iCent == 2) { min = 1e-3; max = 5e4; }
          else if (iCent == 3) { min = 1e-3; max = 4e5; }
        }

        if (!canvasExists) {
          TH1D* htemp = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);

          TAxis* xax = htemp->GetXaxis ();
          TAxis* yax = htemp->GetYaxis ();

          useTrkPt ? xax->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : xax->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
          yax->SetRangeUser (min, max);

          xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
          yax->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x}_{hZ} d#Delta#phi");

          xax->SetTitleSize (0);
          xax->SetLabelSize (0);

          yax->SetTitleFont (43);
          yax->SetTitleSize (axisTextSize);
          yax->SetLabelFont (43);
          yax->SetLabelSize (axisTextSize);

          yax->SetTitleOffset (1.8 * yax->GetTitleOffset ());

          htemp->SetLineWidth (0);
          htemp->SetMarkerSize (0);
  
          htemp->DrawCopy ("hist ][");
          SaferDelete (&htemp);

          if (iCent == 0) {
            myText (0.25, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.06);
            myText (0.25, 0.14, kBlack, "L_{int} = 260 pb^{-1}", 0.06);
            myText (0.25, 0.06, kBlack, "Before subtraction", 0.060);
          }
          else {
            if (iCent == 1) {
              myText (0.25, 0.22, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.06);
              myText (0.25, 0.14, kBlack, "L_{int} = 1.4-1.7 nb^{-1}", 0.06);
            }
            myText (0.25, 0.06, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
          }

          if (iCent == 2 && iSpc == 0)      myText (0.72, 0.87, kBlack, "#it{Z} #rightarrow #it{e}#it{e}", 0.06);
          else if (iCent == 2 && iSpc == 1) myText (0.72, 0.87, kBlack, "#it{Z} #rightarrow #it{#mu}#it{#mu}", 0.06);

          if (iCent == 0) myText (0.485, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.065);
          else if (iCent == 1) {
            if (iPtZ == nPtZBins-1) myText (0.50, 0.87, kBlack, Form ("#it{p}_{T}^{Z} = %g+ GeV", zPtBins[iPtZ]), 0.06);
            else                    myText (0.50, 0.87, kBlack, Form ("#it{p}_{T}^{Z} = %g-%g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.06);
          }
          else if (iCent == numCentBins-2) {
            myText (0.495, 0.88, kBlack, "MB", 0.06);
            myText (0.61, 0.88, kBlack, "|#Delta#phi|", 0.06);
            for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
              myLineText (0.58, 0.82-0.06*(iPhi-1), colors[iPhi], iPhi+1, "", 2.0, 0.054);
              const char* lo = GetPiString (phiLowBins[iPhi]);
              const char* hi = GetPiString (phiHighBins[iPhi]);
              myText (0.61, 0.81-0.06*(iPhi-1), kBlack, Form ("(%s, %s)", lo, hi), 0.054);
            } // end loop over iPhi
          }
          else if (iCent == numCentBins-1) {
            myText (0.48, 0.88, kBlack, "Z-tag", 0.06);
            myText (0.61, 0.88, kBlack, "|#Delta#phi|", 0.06);
            for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
              myMarkerTextNoLine (0.58, 0.82-0.06*(iPhi-1), colors[iPhi], kFullCircle, "", 1.5, 0.054); // for plotting data vs bkg.
              myMarkerTextNoLine (0.58, 0.82-0.06*(iPhi-1), kBlack, kOpenCircle, "", 1.5, 0.054); // for plotting data vs bkg.
              const char* lo = GetPiString (phiLowBins[iPhi]);
              const char* hi = GetPiString (phiHighBins[iPhi]);
              myText (0.61, 0.81-0.06*(iPhi-1), kBlack, Form ("(%s, %s)", lo, hi), 0.054);
            } // end loop over iPhi
          }
        }

        if (plotFill) {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          //for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
            TH1D* h = (useTrkPt ? h_trk_pt_dphi : h_trk_xhz_dphi)[iSpc][iPtZ][iPhi][iCent];

            //h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            //h->SetMarkerSize (0);
            h->SetLineColor (colors[iPhi]);
            h->SetLineStyle (iPhi+1);
            h->SetLineWidth (1);

            h->Draw ("hist ][ same");
          }
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
            TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi : h_trk_xhz_dphi)[iSpc][iPtZ][iPhi][iCent]);
            RecenterGraph (g);
            ResetXErrors (g);
            deltaize (g, 0.08*(iPhi-1 - 0.5*(numPhiBins-2)), true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
            ResetXErrors (g);

            if (plotAsSyst) {
              SetConstantXErrors (g, 0.032, true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
              g->SetMarkerSize (0); 
              g->SetLineWidth (1);
              g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);

              ((TGAE*) g->Clone ())->Draw ("5P");
            }
            else {
              Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);
              double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.0 : 1.5);

              g->SetMarkerStyle (markerStyle);
              g->SetMarkerSize (markerSize);
              g->SetMarkerColor (colors[iPhi]);
              g->SetLineColor (colors[iPhi]);
              g->SetLineWidth (2);

              ((TGAE*) g->Clone ())->Draw ("P");

              markerStyle = FullToOpenMarker (markerStyle);

              g->SetMarkerStyle (markerStyle);
              g->SetMarkerSize (markerSize);
              g->SetLineWidth (0);
              g->SetMarkerColor (kBlack);
  
              ((TGAE*) g->Clone ())->Draw ("P");
            }

            SaferDelete (&g);
          } // end loop over iPhi
        }


        bottomPad->cd ();
        gPad->SetLogx ();
        gPad->SetLogy ();

        if (useTrkPt) {
          min = 3e-3;
          max = 20;
        }
        else {
          min = 1e-2;
          max = 8e2;
        }

        if (!canvasExists) {
          TH1D* htemp = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);

          TAxis* xax = htemp->GetXaxis ();
          TAxis* yax = htemp->GetYaxis ();

          useTrkPt ? xax->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : xax->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
          yax->SetRangeUser (min, max);

          xax->SetMoreLogLabels ();

          xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
          yax->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x}_{hZ} d#Delta#phi");

          xax->SetTitleFont (43);
          xax->SetTitleSize (axisTextSize);
          xax->SetLabelFont (43);
          xax->SetLabelSize (axisTextSize);

          yax->SetTitleFont (43);
          yax->SetTitleSize (axisTextSize);
          yax->SetLabelFont (43);
          yax->SetLabelSize (axisTextSize);

          xax->SetTitleOffset (1.8 * xax->GetTitleOffset ());
          yax->SetTitleOffset (1.8 * yax->GetTitleOffset ());

          htemp->SetLineWidth (0);
          htemp->SetMarkerSize (0);
  
          htemp->DrawCopy ("hist ][");
          SaferDelete (&htemp);

          if (iCent == 0) myText (0.25, 0.24, kBlack, "After subtraction", 0.060/padRatio);
        }

        if (!plotSignal)
          continue;

        if (plotFill) {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            TH1D* h = (useTrkPt ? h_trk_pt_dphi_sub : h_trk_xhz_dphi_sub)[iSpc][iPtZ][iPhi][iCent];

            if (!h) continue;

            h->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            h->SetLineColor (kBlack);
            h->SetMarkerSize (0);
            h->SetLineWidth (0);
            h->SetMarkerStyle (kFullCircle);

            h->DrawCopy ("bar same");
            h->SetLineWidth (1);
            h->Draw ("hist same");
          } // end loop over iPhi
          gPad->RedrawAxis ();
        } else {
          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
            TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi_sub : h_trk_xhz_dphi_sub)[iSpc][iPtZ][iPhi][iCent]);
            RecenterGraph (g);
            ResetXErrors (g);
            deltaize (g, 0.08*(iPhi-1 - 0.5*(numPhiBins-2)), true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
            ResetXErrors (g);

            if (plotAsSyst) {
              SetConstantXErrors (g, 0.032, true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
              g->SetMarkerSize (0); 
              g->SetLineWidth (1);
              g->SetLineColor (colors[iPhi]);
              g->SetFillColorAlpha (fillColors[iPhi], 0.3);

              ((TGAE*) g->Clone ())->Draw ("5P");
            }
            else {
              Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);
              double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.0 : 1.5);

              g->SetMarkerStyle (markerStyle);
              g->SetMarkerSize (markerSize);
              g->SetMarkerColor (colors[iPhi]);
              g->SetLineColor (colors[iPhi]);
              g->SetLineWidth (2);

              ((TGAE*) g->Clone ())->Draw ("P");

              markerStyle = FullToOpenMarker (markerStyle);

              g->SetMarkerStyle (markerStyle);
              g->SetMarkerSize (markerSize);
              g->SetLineWidth (0);
              g->SetMarkerColor (kBlack);
  
              ((TGAE*) g->Clone ())->Draw ("P");
            }

            SaferDelete (&g);
          } // end loop over iPhi
        }
      } // end loop over iCent
      
      c->SaveAs (Form ("%s/TrkYields/allYields_%s_iPtZ%i_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in Z Pt
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotAllYields_dPtZ (const bool useTrkPt, const short pSpc) {
  const double padRatio = 0.9; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (padRatio+1.0);
  const int axisTextSize = 23;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_AllYields_%s_dPtZ_%s", useTrkPt ? "pTch" : "xhZ", spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 395*numCentBins, 685);
      gDirectory->Add (c);
    }

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      c->cd ();

      const char* topPadName = Form ("p_AllYields_top_%s_dPtZ_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", iCent, spc);
      const char* bottomPadName = Form ("p_AllYields_bottom_%s_dPtZ_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", iCent, spc);

      TPad* topPad = nullptr, *bottomPad = nullptr;
      if (!canvasExists) {
        topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
        bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

        gDirectory->Add (topPad);
        gDirectory->Add (bottomPad);

        topPad->SetTopMargin (0.04);
        topPad->SetBottomMargin (0.03);
        topPad->SetLeftMargin (0.23);
        topPad->SetRightMargin (0.03);
        bottomPad->SetTopMargin (0.02);
        bottomPad->SetBottomMargin (0.23);
        bottomPad->SetLeftMargin (0.23);
        bottomPad->SetRightMargin (0.03);
        topPad->Draw ();
        bottomPad->Draw ();
     }
     else {
        topPad = dynamic_cast<TPad*> (gDirectory->Get (topPadName));
        bottomPad = dynamic_cast<TPad*> (gDirectory->Get (bottomPadName));
      }

      topPad->cd ();
      GetDrawnObjects ();
      gPad->SetLogx ();
      gPad->SetLogy ();

      double min = 1e30, max = 0;
      if (useTrkPt) {
        if (iCent == 0)       { min = 1e-5; max = 4e1; }
        else if (iCent == 1)  { min = 1e-5; max = 4e2; }
        else if (iCent == 2)  { min = 5e-5; max = 2e3; }
        else if (iCent == 3)  { min = 8e-5; max = 6e3; }
      }
      else {
        if (iCent == 0)       { min = 2e-5; max = 2e3; }
        else if (iCent == 1)  { min = 1e-2; max = 5e4; }
        else if (iCent == 2)  { min = 1e-2; max = 2e5; }
        else if (iCent == 3)  { min = 1e-2; max = 9e5; }
      }

      if (!canvasExists) {
        TH1D* htemp = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);

        TAxis* xax = htemp->GetXaxis ();
        TAxis* yax = htemp->GetYaxis ();

        yax->SetRangeUser (min, max);

        xax->SetMoreLogLabels ();

        yax->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x} d#Delta#phi");

        xax->SetTitleSize (0);
        xax->SetLabelSize (0);

        yax->SetTitleFont (43);
        yax->SetTitleSize (axisTextSize);
        yax->SetLabelFont (43);
        yax->SetLabelSize (axisTextSize);

        yax->SetTitleOffset (1.8 * yax->GetTitleOffset ());

        htemp->SetLineWidth (0);
        htemp->SetMarkerSize (0);
  
        htemp->DrawCopy ("hist ][");
        SaferDelete (&htemp);

        if (iCent == 0) {
          myText (0.25, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.06);
          myText (0.25, 0.14, kBlack, "L_{int} = 260 pb^{-1}", 0.06);
          myText (0.25, 0.06, kBlack, "Before subtraction", 0.060);
        }
        else {
          if (iCent == 1) {
            myText (0.25, 0.22, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.06);
            myText (0.25, 0.14, kBlack, "L_{int} = 1.4-1.7 nb^{-1}", 0.06);
          }
          myText (0.25, 0.06, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
        }

        if (iCent == 2 && iSpc == 0)      myText (0.72, 0.87, kBlack, "#it{Z} #rightarrow #it{e}#it{e}", 0.06);
        else if (iCent == 2 && iSpc == 1) myText (0.72, 0.87, kBlack, "#it{Z} #rightarrow #it{#mu}#it{#mu}", 0.06);

        if (iCent == 0)       myText (0.485, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.065);
        else if (iCent == 1)  myText (0.485, 0.87, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.06);
        else if (iCent == numCentBins-2) {
          myText (0.59, 0.88, kBlack, "MB", 0.06);
          myText (0.69, 0.88, kBlack, "#it{p}_{T}^{Z} [GeV]", 0.06);
          for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
            myLineText (0.67, 0.82-0.06*(iPtZ-2), finalColors[iPtZ-1], iPtZ, "", 2.0, 0.054) ;
            if (iPtZ == nPtZBins-1) myText (0.69, 0.81-0.06*(iPtZ-2), kBlack, Form ("%g+", zPtBins[iPtZ]), 0.054);
            else                    myText (0.69, 0.81-0.06*(iPtZ-2), kBlack, Form ("%g-%g", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.054);
          } // end loop over iPtZ
        }
        else if (iCent == numCentBins-1) {
          //myText (0.38, 0.90, kBlack, "MB", 0.06);
          myText (0.55, 0.88, kBlack, "Z-tag", 0.06);
          myText (0.69, 0.88, kBlack, "#it{p}_{T}^{Z} [GeV]", 0.06);
          for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
            myMarkerTextNoLine (0.65, 0.82-0.06*(iPtZ-2), finalColors[iPtZ-1], markerStyles[iPtZ-2], "", 1.5, 0.054); // for plotting data vs bkg.
            myMarkerTextNoLine (0.65, 0.82-0.06*(iPtZ-2), kBlack, FullToOpenMarker (markerStyles[iPtZ-2]), "", 1.5, 0.054); // for plotting data vs bkg.
            if (iPtZ == nPtZBins-1) myText (0.69, 0.81-0.06*(iPtZ-2), kBlack, Form ("%g+", zPtBins[iPtZ]), 0.054);
            else                    myText (0.69, 0.81-0.06*(iPtZ-2), kBlack, Form ("%g-%g", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.054);
          } // end loop over iPtZ
        }
      }

      if (plotFill) {
        for (int iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
          if (!subtractPP && iCent == 0) continue;

          TH1D* h = (useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[iSpc][iPtZ][iCent];

          //h->SetFillColorAlpha (fillColors[iPtZ-2], fillAlpha);
          //h->SetMarkerSize (0);
          h->SetLineColor (finalColors[iPtZ-1]);
          h->SetLineStyle (iPtZ);
          h->SetLineWidth (1);

          h->Draw ("hist ][ same");
        }
        gPad->RedrawAxis ();
      }
      else {
        for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[iSpc][iPtZ][iCent]);
          RecenterGraph (g);
          ResetXErrors (g);
          deltaize (g, 0.08*(iPtZ-2 - 0.5*(nPtZBins-3)), true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
          ResetXErrors (g);

          if (plotAsSyst) {
            SetConstantXErrors (g, 0.032, true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (finalColors[iPtZ-1]);
            g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

            ((TGAE*) g->Clone ())->Draw ("5P");
          }
          else {
            Style_t markerStyle = markerStyles[iPtZ-2];
            double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.0 : 1.5);

            g->SetMarkerStyle (markerStyle);
            g->SetMarkerSize (markerSize);
            g->SetMarkerColor (finalColors[iPtZ-1]);
            g->SetLineColor (finalColors[iPtZ-1]);
            g->SetLineWidth (2);

            ((TGAE*) g->Clone ())->Draw ("P");

            if (IsFullMarker (markerStyle)) {
              markerStyle = FullToOpenMarker (markerStyle);

              g->SetMarkerStyle (markerStyle);
              g->SetMarkerSize (markerSize);
              g->SetLineWidth (0);
              g->SetMarkerColor (kBlack);
  
              ((TGAE*) g->Clone ())->Draw ("P");
            }
          }

          SaferDelete (&g);
        } // end loop over iPtZ
      }


      bottomPad->cd ();
      GetDrawnObjects ();
      gPad->SetLogx ();
      gPad->SetLogy ();

      if (useTrkPt) { min = 1e-3; max = 2e1; }
      else          { min = 1e-2; max = 8e2; }

      if (!canvasExists) {
        TH1D* htemp = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);

        TAxis* xax = htemp->GetXaxis ();
        TAxis* yax = htemp->GetYaxis ();

        yax->SetRangeUser (min, max);

        xax->SetMoreLogLabels ();

        xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
        yax->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x} d#Delta#phi");

        xax->SetTitleFont (43);
        xax->SetTitleSize (axisTextSize);
        xax->SetLabelFont (43);
        xax->SetLabelSize (axisTextSize);

        yax->SetTitleFont (43);
        yax->SetTitleSize (axisTextSize);
        yax->SetLabelFont (43);
        yax->SetLabelSize (axisTextSize);

        xax->SetTitleOffset (1.8 * xax->GetTitleOffset ());
        yax->SetTitleOffset (1.8 * yax->GetTitleOffset ());

        htemp->SetLineWidth (0);
        htemp->SetMarkerSize (0);
  
        htemp->DrawCopy ("hist ][");
        SaferDelete (&htemp);

        if (iCent == 0) myText (0.25, 0.24, kBlack, "After subtraction", 0.060/padRatio);
      }

      if (!plotSignal)
        continue;

      if (plotFill) {
        for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          TH1D* h = (useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[iSpc][iPtZ][iCent];

          if (!h) continue;

          h->SetFillColorAlpha (finalColors[iPtZ-1], fillAlpha);
          h->SetLineColor (kBlack);
          h->SetMarkerSize (0);
          h->SetLineWidth (0);
          h->SetMarkerStyle (kFullCircle);

          h->DrawCopy ("bar same");
          h->SetLineWidth (1);
          h->Draw ("hist same");
        } // end loop over iPtZ
        gPad->RedrawAxis ();
      } else {
        for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[iSpc][iPtZ][iCent]);
          RecenterGraph (g);
          ResetXErrors (g);
          deltaize (g, 0.08*(iPtZ-2 - 0.5*(nPtZBins-3)), true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
          ResetXErrors (g);

          if (plotAsSyst) {
            SetConstantXErrors (g, 0.032, true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (finalColors[iPtZ-1]);
            g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

            ((TGAE*) g->Clone ())->Draw ("5P");
          }
          else {
            Style_t markerStyle = markerStyles[iPtZ-2];
            double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.0 : 1.5);

            g->SetMarkerStyle (markerStyle);
            g->SetMarkerSize (markerSize);
            g->SetMarkerColor (finalColors[iPtZ-1]);
            g->SetLineColor (finalColors[iPtZ-1]);
            g->SetLineWidth (2);

            ((TGAE*) g->Clone ())->Draw ("P");

            if (IsFullMarker (markerStyle)) {
              markerStyle = FullToOpenMarker (markerStyle);

              g->SetMarkerStyle (markerStyle);
              g->SetMarkerSize (markerSize);
              g->SetLineWidth (0);
              g->SetMarkerColor (kBlack);
  
              ((TGAE*) g->Clone ())->Draw ("P");
            }
          }

          SaferDelete (&g);
        } // end loop over iPtZ
      }
    } // end loop over iCent
    
    c->SaveAs (Form ("%s/TrkYields/allYields_%s_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", spc));
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in Z Pt compared with truth yields
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotClosureCheck (PhysicsAnalysis* truthComp, const bool useTrkPt, const short iSpc) {
  const double padRatio = 0.9; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (padRatio+1.0);
  const int axisTextSize = 23;

  const char* canvasName = Form ("c_ClosureCheck_%s_dPtZ", useTrkPt ? "pTch" : "xhZ");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 395*numCentBins, 685);
    gDirectory->Add (c);
  }

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    const char* topPadName = Form ("p_ClosureCheck_top_%s_dPtZ_iCent%i", useTrkPt ? "pTch" : "xhZ", iCent);
    const char* bottomPadName = Form ("p_ClosureCheck_bottom_%s_dPtZ_iCent%i", useTrkPt ? "pTch" : "xhZ", iCent);

    TPad* topPad = nullptr, *bottomPad = nullptr;
    if (!canvasExists) {
      topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
      bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

      gDirectory->Add (topPad);
      gDirectory->Add (bottomPad);

      topPad->SetTopMargin (0.04);
      topPad->SetBottomMargin (0.03);
      topPad->SetLeftMargin (0.23);
      topPad->SetRightMargin (0.03);
      bottomPad->SetTopMargin (0.02);
      bottomPad->SetBottomMargin (0.20);
      bottomPad->SetLeftMargin (0.23);
      bottomPad->SetRightMargin (0.03);
      topPad->Draw ();
      bottomPad->Draw ();
   }
   else {
      topPad = dynamic_cast<TPad*> (gDirectory->Get (topPadName));
      bottomPad = dynamic_cast<TPad*> (gDirectory->Get (bottomPadName));
    }

    topPad->cd ();
    GetDrawnObjects ();
    gPad->SetLogx ();
    gPad->SetLogy ();

    double min = 1e30, max = 0;
    if (useTrkPt) { min = 1e-3; max = 2e1; }
    else          { min = 1e-2; max = 1.2e3; }

    if (!canvasExists) {
      TH1D* htemp = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);

      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      yax->SetRangeUser (min, max);

      xax->SetMoreLogLabels ();

      yax->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x} d#Delta#phi");

      xax->SetTitleSize (0);
      xax->SetLabelSize (0);

      yax->SetTitleFont (43);
      yax->SetTitleSize (axisTextSize);
      yax->SetLabelFont (43);
      yax->SetLabelSize (axisTextSize);

      yax->SetTitleOffset (1.8 * yax->GetTitleOffset ());

      htemp->SetLineWidth (0);
      htemp->SetMarkerSize (0);
  
      htemp->DrawCopy ("hist ][");
      SaferDelete (&htemp);

      if (iCent == 0) myText (0.25, 0.14, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.06);
      else {
        if (iCent == 1) myText (0.25, 0.14, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.06);
        myText (0.25, 0.06, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
      }

      if (iCent == 0)       myText (0.245, 0.87, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.065);
      else if (iCent == 1)  myText (0.485, 0.87, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.06);
      else if (iCent == numCentBins-2) {
        myText (0.57, 0.88, kBlack, "Truth", 0.06);
        myText (0.69, 0.88, kBlack, "#it{p}_{T}^{Z} [GeV]", 0.06);
        for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          myLineText (0.67, 0.82-0.06*(iPtZ-2), finalColors[iPtZ-1], iPtZ, "", 2.0, 0.054) ;
          if (iPtZ == nPtZBins-1) myText (0.69, 0.81-0.06*(iPtZ-2), kBlack, Form ("%g+", zPtBins[iPtZ]), 0.054);
          else                    myText (0.69, 0.81-0.06*(iPtZ-2), kBlack, Form ("%g-%g", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.054);
        } // end loop over iPtZ
      }
      else if (iCent == numCentBins-1) {
        //myText (0.38, 0.90, kBlack, "MB", 0.06);
        myText (0.55, 0.88, kBlack, "Reco.", 0.06);
        myText (0.69, 0.88, kBlack, "#it{p}_{T}^{Z} [GeV]", 0.06);
        for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          myMarkerTextNoLine (0.65, 0.82-0.06*(iPtZ-2), finalColors[iPtZ-1], markerStyles[iPtZ-2], "", 1.5, 0.054); // for plotting data vs bkg.
          myMarkerTextNoLine (0.65, 0.82-0.06*(iPtZ-2), kBlack, FullToOpenMarker (markerStyles[iPtZ-2]), "", 1.5, 0.054); // for plotting data vs bkg.
          if (iPtZ == nPtZBins-1) myText (0.69, 0.81-0.06*(iPtZ-2), kBlack, Form ("%g+", zPtBins[iPtZ]), 0.054);
          else                    myText (0.69, 0.81-0.06*(iPtZ-2), kBlack, Form ("%g-%g", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.054);
        } // end loop over iPtZ
      }
    }

    if (!plotAsSyst) {
      for (int iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
        TH1D* h = (useTrkPt ? truthComp->h_trk_pt_ptz : truthComp->h_trk_xhz_ptz)[2][iPtZ][iCent];

        //h->SetFillColorAlpha (fillColors[iPtZ-2], fillAlpha);
        //h->SetMarkerSize (0);
        h->SetLineColor (finalColors[iPtZ-1]);
        h->SetLineStyle (iPtZ);
        h->SetLineWidth (1);

        h->Draw ("hist ][ same");
      }
      gPad->RedrawAxis ();
    }
    
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      TGAE* g = GetTGAE ((useTrkPt ? (iCent == 0 ? h_trk_pt_ptz : h_trk_pt_ptz_sub) : (iCent == 0 ? h_trk_xhz_ptz : h_trk_xhz_ptz_sub))[iSpc][iPtZ][iCent]);
      RecenterGraph (g);
      ResetXErrors (g);
      deltaize (g, 0.08*(iPtZ-2 - 0.5*(nPtZBins-3)), true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
      ResetXErrors (g);

      if (plotAsSyst) {
        SetConstantXErrors (g, 0.032, true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
        g->SetMarkerSize (0); 
        g->SetLineWidth (1);
        g->SetLineColor (finalColors[iPtZ-1]);
        g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

        if (iCent != 0) ((TGAE*) g->Clone ())->Draw ("5P");
      }
      else {
        Style_t markerStyle = markerStyles[iPtZ-2];
        double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.0 : 1.5);

        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (markerSize);
        g->SetMarkerColor (finalColors[iPtZ-1]);
        g->SetLineColor (finalColors[iPtZ-1]);
        g->SetLineWidth (2);

        ((TGAE*) g->Clone ())->Draw ("P");

        if (IsFullMarker (markerStyle)) {
          markerStyle = FullToOpenMarker (markerStyle);

          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (markerSize);
          g->SetLineWidth (0);
          g->SetMarkerColor (kBlack);
  
          ((TGAE*) g->Clone ())->Draw ("P");
        }
      }

      SaferDelete (&g);
    } // end loop over iPtZ
    


    bottomPad->cd ();
    GetDrawnObjects ();
    gPad->SetLogx ();

    if (useTrkPt) { min = 0.6; max = 1.4; }
    else          { min = 0.6; max = 1.4; }

    if (!canvasExists) {
      TH1D* htemp = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);

      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      //yax->SetRangeUser (min, max);
      yax->SetRangeUser (iCent == 0 ? 0.9 : 0.6, iCent == 0 ? 1.1 : 1.4);
      //yax->SetRangeUser (iCent == 0 ? -2 : -2, iCent == 0 ? 1 : 1);
      //yax->SetRangeUser (iCent == 0 ? 0.9 : 0.75, iCent == 0 ? 1.1 : 1.25);

      xax->SetMoreLogLabels ();

      xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
      yax->SetTitle ("Reco. level / Pythia Truth");
      yax->CenterTitle ();

      xax->SetTitleFont (43);
      xax->SetTitleSize (axisTextSize);
      xax->SetLabelFont (43);
      xax->SetLabelSize (axisTextSize);

      yax->SetTitleFont (43);
      yax->SetTitleSize (axisTextSize);
      yax->SetLabelFont (43);
      yax->SetLabelSize (axisTextSize);

      xax->SetTitleOffset (1.8 * xax->GetTitleOffset ());
      yax->SetTitleOffset (1.8 * yax->GetTitleOffset ());

      htemp->SetLineWidth (0);
      htemp->SetMarkerSize (0);
  
      htemp->DrawCopy ("hist ][");
      SaferDelete (&htemp);

      const float yerr = (iCent == 0 ? 0.02 : 0.05);

      TGAE* g = new TGAE ();
      g->SetPoint (g->GetN (), (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1][0], 1.);
      g->SetPointEYhigh (g->GetN () - 1, yerr);
      g->SetPointEYlow (g->GetN () - 1, yerr);
      g->SetPoint (g->GetN (), (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1][(useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1]], 1.);
      g->SetPointEYhigh (g->GetN () - 1, yerr);
      g->SetPointEYlow (g->GetN () - 1, yerr);
      g->SetFillColorAlpha (kGray+1, 0.4);
      g->Draw ("3");

      TLine* dashes = new TLine ();
      dashes->SetLineStyle (2);
      dashes->SetLineColor (kGray+1);
      dashes->DrawLine ((useTrkPt ? pTchBins : xhZBins)[nPtZBins-1][0], 1+yerr, (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1][(useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1]], 1+yerr);
      dashes->DrawLine ((useTrkPt ? pTchBins : xhZBins)[nPtZBins-1][0], 1, (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1][(useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1]], 1);
      dashes->DrawLine ((useTrkPt ? pTchBins : xhZBins)[nPtZBins-1][0], 1-yerr, (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1][(useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1]], 1-yerr);
      dashes->SetLineColor (kBlack);
      dashes->DrawLine ((useTrkPt ? pTchBins : xhZBins)[nPtZBins-1][0], 0, (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1][(useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1]], 0);
    }

    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      TH1D* h = (useTrkPt ? (iCent == 0 ? h_trk_pt_ptz : h_trk_pt_ptz_sub) : (iCent == 0 ? h_trk_xhz_ptz : h_trk_xhz_ptz_sub))[iSpc][iPtZ][iCent];
      TH1D* hden = (useTrkPt ? truthComp->h_trk_pt_ptz : truthComp->h_trk_xhz_ptz)[2][iPtZ][iCent];

      TGAE* g = (TGAE*) GetTGAE (h)->Clone ();

      double _x, _y;
      for (int ix = 0; ix < g->GetN (); ix++) {
        g->GetPoint (ix, _x, _y);
        int _ix = ix+1;
        while (_ix <= hden->GetNbinsX () && hden->GetBinCenter (_ix) != _x) _ix++;
        g->SetPoint (ix, _x, (_y / hden->GetBinContent (_ix)));
        //g->SetPointEYhigh (ix, sqrt (pow (g->GetErrorYhigh (ix), 2) + pow (hden->GetBinError (_ix), 2)) * hden->GetBinWidth (_ix));
        //g->SetPointEYlow (ix, sqrt (pow (g->GetErrorYlow (ix), 2) + pow (hden->GetBinError (_ix), 2)) * hden->GetBinWidth (_ix));
        g->SetPointEYhigh (ix, fabs (_y / hden->GetBinContent (_ix)) * sqrt (pow (g->GetErrorYhigh (ix) / _y, 2) + pow (hden->GetBinError (_ix) / hden->GetBinContent (_ix), 2)));
        g->SetPointEYlow (ix, fabs (_y / hden->GetBinContent (_ix)) * sqrt (pow (g->GetErrorYlow (ix) / _y, 2) + pow (hden->GetBinError (_ix) / hden->GetBinContent (_ix), 2)));
      }

      RecenterGraph (g);
      ResetXErrors (g);
      deltaize (g, 0.08*(iPtZ-2 - 0.5*(nPtZBins-3)), true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
      ResetXErrors (g);

      if (plotAsSyst) {
        SetConstantXErrors (g, 0.032, true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
        g->SetMarkerSize (0); 
        g->SetLineWidth (1);
        g->SetLineColor (finalColors[iPtZ-1]);
        g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

        if (iCent != 0) ((TGAE*) g->Clone ())->Draw ("5P");
      }
      else {
        Style_t markerStyle = markerStyles[iPtZ-2];
        double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.0 : 1.5);

        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (markerSize);
        g->SetMarkerColor (finalColors[iPtZ-1]);
        g->SetLineColor (finalColors[iPtZ-1]);
        g->SetLineWidth (2);

        ((TGAE*) g->Clone ())->Draw ("P");

        if (IsFullMarker (markerStyle)) {
          markerStyle = FullToOpenMarker (markerStyle);

          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (markerSize);
          g->SetLineWidth (0);
          g->SetMarkerColor (kBlack);
  
          ((TGAE*) g->Clone ())->Draw ("P");
        }
      }

      SaferDelete (&g);
    } // end loop over iPtZ
  } // end loop over iCent
  
  c->SaveAs (Form ("%s/TrkYields/closureCheck_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in Z Pt compared with truth yields as a function of S/B ratio
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotClosureCheck_SigToBkg (PhysicsAnalysis* truthComp, const bool useTrkPt, const short iSpc, const char* saveFileName) {
  CalculateSigToBkg ();
  const int axisTextSize = 36;
  const int axisLabelSize = 36;

  const char* canvasName = Form ("c_ClosureCheck_SigToBkg_%s_dPtZ", useTrkPt ? "pTch" : "xhZ");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 800);
    gDirectory->Add (c);
  }

  c->cd ();

  gPad->SetLogx ();


  if (!canvasExists) {
    TH1D* htemp = new TH1D ("htemp", "", 1, 5e-3, 1e2);

    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    yax->SetRangeUser (0.6, 1.4);

    //xax->SetMoreLogLabels ();

    xax->SetTitle (Form ("Y_{sig} (%s) / Y_{bkg} (%s)", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}"));
    yax->SetTitle ("Reco. level / Pythia Truth");

    xax->SetTitleFont (43);
    xax->SetTitleSize (axisTextSize);
    xax->SetLabelFont (43);
    xax->SetLabelSize (axisLabelSize);

    yax->SetTitleFont (43);
    yax->SetTitleSize (axisTextSize);
    yax->SetLabelFont (43);
    yax->SetLabelSize (axisLabelSize);

    //xax->SetTitleOffset (1.0 * xax->GetTitleOffset ());
    //yax->SetTitleOffset (1.0 * yax->GetTitleOffset ());

    htemp->SetLineWidth (0);
    htemp->SetMarkerSize (0);
  
    htemp->DrawCopy ("hist ][");
    SaferDelete (&htemp);

    const float yerr = 0.05;

    TGAE* g = new TGAE ();
    g->SetPoint (g->GetN (), 5e-3, 1.);
    g->SetPointEYhigh (g->GetN () - 1, yerr);
    g->SetPointEYlow (g->GetN () - 1, yerr);
    g->SetPoint (g->GetN (), 1e2, 1.);
    g->SetPointEYhigh (g->GetN () - 1, yerr);
    g->SetPointEYlow (g->GetN () - 1, yerr);
    g->SetFillColorAlpha (kGray+1, 0.4);
    g->Draw ("3");

    TLine* dashes = new TLine ();
    dashes->SetLineStyle (2);
    dashes->SetLineColor (kGray+1);
    dashes->DrawLine (5e-3, 1+yerr, 1e2, 1+yerr);
    dashes->DrawLine (5e-3, 1, 1e2, 1);
    dashes->DrawLine (5e-3, 1-yerr, 1e2, 1-yerr);
    dashes->SetLineColor (kBlack);

    myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.038);
    myText (0.22, 0.83, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.79, kBlack, "PowhegPythia + Hijing Overlay", 0.032);
    myText (0.22, 0.75, kBlack, "All centralities", 0.032);

    myMarkerAndBoxAndLineText (0.67, 0.320, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "15 < #it{p}_{T}^{Z} < 30 GeV", 0.036);
    myMarkerAndBoxAndLineText (0.67, 0.270, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "30 < #it{p}_{T}^{Z} < 60 GeV", 0.036);
    myMarkerAndBoxAndLineText (0.67, 0.220, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "#it{p}_{T}^{Z} > 60 GeV", 0.036);

  }

  TGAE* g_allPoints = new TGAE ();

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      TH1D* h = (useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[iSpc][iPtZ][iCent];
      TH1D* hden = (useTrkPt ? truthComp->h_trk_pt_ptz : truthComp->h_trk_xhz_ptz)[2][iPtZ][iCent];

      TGAE* g = (TGAE*) GetTGAE (h)->Clone ();

      double _x, _y;
      for (int ix = 0; ix < g->GetN (); ix++) {
        g->GetPoint (ix, _x, _y);
        int _ix = ix+1;
        while (_ix <= hden->GetNbinsX () && hden->GetBinCenter (_ix) != _x) _ix++;
        g->SetPoint (ix, _x, (_y / hden->GetBinContent (_ix)));
        //g->SetPointEYhigh (ix, sqrt (pow (g->GetErrorYhigh (ix), 2) + pow (hden->GetBinError (_ix), 2)) * hden->GetBinWidth (_ix));
        //g->SetPointEYlow (ix, sqrt (pow (g->GetErrorYlow (ix), 2) + pow (hden->GetBinError (_ix), 2)) * hden->GetBinWidth (_ix));
        g->SetPointEYhigh (ix, fabs (_y / hden->GetBinContent (_ix)) * sqrt (pow (g->GetErrorYhigh (ix) / _y, 2) + pow (hden->GetBinError (_ix) / hden->GetBinContent (_ix), 2)));
        g->SetPointEYlow (ix, fabs (_y / hden->GetBinContent (_ix)) * sqrt (pow (g->GetErrorYlow (ix) / _y, 2) + pow (hden->GetBinError (_ix) / hden->GetBinContent (_ix), 2)));
      }

      TH1D* h_s2b = (useTrkPt ? h_trk_pt_ptz_sig_to_bkg : h_trk_xhz_ptz_sig_to_bkg)[iSpc][iPtZ][iCent];
      double x, y;
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y); 
        double s2b = -1, s2b_err = -1;
        for (int iX = 1; s2b == -1 && iX <= h_s2b->GetNbinsX (); iX++) {
          if (h_s2b->GetBinLowEdge (iX) <= x && x <= h_s2b->GetBinLowEdge (iX) + h_s2b->GetBinWidth (iX)) {
            s2b = h_s2b->GetBinContent (iX);
            s2b_err = h_s2b->GetBinError (iX);
          } 
        }
        if (isnan (s2b) || s2b <= 0) continue;
        g->SetPoint (i, s2b, y);
        g->SetPointEXhigh (i, s2b_err);
        g->SetPointEXlow (i, s2b_err);
      }

      if (plotAsSyst) {
        SetConstantXErrors (g, 0.032, true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
        g->SetMarkerSize (0); 
        g->SetLineWidth (1);
        g->SetLineColor (finalColors[iPtZ-1]);
        g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

        if (iCent != 0) ((TGAE*) g->Clone ())->Draw ("5P");
      }
      else {
        Style_t markerStyle = markerStyles[iPtZ-2];
        double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.0 : 1.5);

        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (markerSize);
        g->SetMarkerColor (finalColors[iPtZ-1]);
        g->SetLineColor (finalColors[iPtZ-1]);
        g->SetLineWidth (2);

        ((TGAE*) g->Clone ())->Draw ("P");

        if (IsFullMarker (markerStyle)) {
          markerStyle = FullToOpenMarker (markerStyle);

          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (markerSize);
          g->SetLineWidth (0);
          g->SetMarkerColor (kBlack);
  
          ((TGAE*) g->Clone ())->Draw ("P");
        }

        double x, y;
        for (int i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);

          g_allPoints->SetPoint (g_allPoints->GetN (), x, y);
          g_allPoints->SetPointEXlow (g_allPoints->GetN ()-1, g->GetErrorXlow (i));
          g_allPoints->SetPointEXhigh (g_allPoints->GetN ()-1, g->GetErrorXhigh (i));
          g_allPoints->SetPointEYlow (g_allPoints->GetN ()-1, g->GetErrorYlow (i));
          g_allPoints->SetPointEYhigh (g_allPoints->GetN ()-1, g->GetErrorYhigh (i));
        }

      } 
      SaferDelete (&g);
    } // end loop over iPtZ
  } // end loop over iCent

  
  if (!plotAsSyst) {
    //TF1* fit1 = new TF1 ("f_closure_sigToBkg_1", "[0]+[1]/x", 5e-3, 2e2);
    //fit1->SetParameter (0, 1);
    //fit1->SetParameter (1, 0);
    //g_allPoints->Fit (fit1, "RN0Q");
    //fit1->SetLineColor (kBlack);
    //fit1->Draw ("same");

    //TF1* fit2 = new TF1 ("f_closure_sigToBkg_2", "[0]+[1]/([2]+x)", 5e-3, 2e2);
    //fit2->SetParameter (0, 1);
    //fit2->SetParameter (1, 0);
    //fit2->SetParameter (2, 0);
    //g_allPoints->Fit (fit2, "RN0Q");
    //fit2->SetLineColor (kBlue);
    //fit2->Draw ("same");

    //TF1* fit3 = new TF1 ("f_closure_sigToBkg_3", "[0]+[1]/x + [2]/(x^2)", 5e-3, 2e2);
    //fit3->SetParameter (0, 1);
    //fit3->SetParameter (1, 0);
    //fit3->SetParameter (2, 0);
    //g_allPoints->Fit (fit3, "RN0Q");
    //fit3->SetLineColor (kMagenta);
    //fit3->Draw ("same");

    TF1* fit = new TF1 (Form ("f_closure_sigToBkg_%s", useTrkPt ? "pTch" : "xhZ"), "[0]+[1]*(x^[2])+[3]*(x^[4])", 5e-3, 2e2);
    fit->SetParameter (0, 1);
    fit->SetParameter (1, 0);
    fit->SetParameter (2, -1);
    fit->SetParameter (3, 0);
    fit->SetParameter (4, -1);
    g_allPoints->Fit (fit, "RN0");
    fit->SetLineColor (kBlack);
    fit->Draw ("same");

    if (!canvasExists) {
      myLineText (0.25, 0.71, kBlack, 1, "Fit to [0]+[1]/x^[2]+[3]/x^[4]", 1.4, 0.028);
      myText (0.25, 0.67,  kBlack, Form ("#chi^{2}/dof = %.2f/%i", fit->GetChisquare (), fit->GetNDF ()), 0.028);
    }

    if (string (saveFileName) != "") {
      TFile* saveFile = new TFile (saveFileName, "recreate");
      fit->Write ();
      saveFile->Close ();
      SaferDelete (&saveFile);
    }
  }
  delete g_allPoints;
  
  c->SaveAs (Form ("%s/TrkYields/closureCheck_SigToBkg_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in Z Pt compared with truth yields as a function of S/B ratio
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotClosureCheck_SigToBkg_BkgNormalization (PhysicsAnalysis* truthComp, const bool useTrkPt, const short iSpc, const char* saveFileName) {
  CalculateSigToBkg ();
  const int axisTextSize = 36;
  const int axisLabelSize = 36;

  const char* canvasName = Form ("c_ClosureCheck_SigToBkg_BkgNormalization_%s_dPtZ", useTrkPt ? "pTch" : "xhZ");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 800);
    gDirectory->Add (c);
  }

  c->cd ();

  gPad->SetLogx ();

  if (!canvasExists) {
    TH1D* htemp = new TH1D ("htemp", "", 1, 5e-3, 1e2);

    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    yax->SetRangeUser (-10, 10);

    xax->SetTitle (Form ("Y_{sig} (%s) / Y_{bkg} (%s)", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}"));
    yax->SetTitle ("Bkg. normalization uncertainty [%]");

    xax->SetTitleFont (43);
    xax->SetTitleSize (axisTextSize);
    xax->SetLabelFont (43);
    xax->SetLabelSize (axisLabelSize);

    yax->SetTitleFont (43);
    yax->SetTitleSize (axisTextSize);
    yax->SetLabelFont (43);
    yax->SetLabelSize (axisLabelSize);

    //xax->SetTitleOffset (1.0 * xax->GetTitleOffset ());
    //yax->SetTitleOffset (1.0 * yax->GetTitleOffset ());

    htemp->SetLineWidth (0);
    htemp->SetMarkerSize (0);
  
    htemp->DrawCopy ("hist ][");
    SaferDelete (&htemp);

    const float yerr = 0.3;

    TGAE* g = new TGAE ();
    g->SetPoint (g->GetN (), 5e-3, 0.);
    g->SetPointEYhigh (g->GetN () - 1, yerr);
    g->SetPointEYlow (g->GetN () - 1, yerr);
    g->SetPoint (g->GetN (), 1e2, 0.);
    g->SetPointEYhigh (g->GetN () - 1, yerr);
    g->SetPointEYlow (g->GetN () - 1, yerr);
    g->SetFillColorAlpha (kGray+1, 0.4);
    g->Draw ("3");

    TLine* dashes = new TLine ();
    dashes->SetLineStyle (2);
    dashes->SetLineColor (kGray+1);
    dashes->DrawLine (5e-3, yerr, 1e2, yerr);
    dashes->DrawLine (5e-3, 0, 1e2, 0);
    dashes->DrawLine (5e-3, -yerr, 1e2, -yerr);
    dashes->SetLineColor (kBlack);

    myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.038);
    myText (0.22, 0.83, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.22, 0.79, kBlack, "PowhegPythia + Hijing Overlay", 0.032);
    myText (0.22, 0.75, kBlack, "All centralities", 0.032);

    myMarkerAndBoxAndLineText (0.25, 0.320, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "15 < #it{p}_{T}^{Z} < 30 GeV", 0.036);
    myMarkerAndBoxAndLineText (0.25, 0.270, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "30 < #it{p}_{T}^{Z} < 60 GeV", 0.036);
    myMarkerAndBoxAndLineText (0.25, 0.220, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "#it{p}_{T}^{Z} > 60 GeV", 0.036);

  }

  TGAE* g_allPoints = new TGAE ();

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      TH1D* h = (useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[iSpc][iPtZ][iCent];
      TH1D* hden = (useTrkPt ? truthComp->h_trk_pt_ptz : truthComp->h_trk_xhz_ptz)[2][iPtZ][iCent];

      TGAE* g = (TGAE*) GetTGAE (h)->Clone ();

      double _x, _y;
      for (int ix = 0; ix < g->GetN (); ix++) {
        g->GetPoint (ix, _x, _y);
        int _ix = ix+1;
        while (_ix <= hden->GetNbinsX () && hden->GetBinCenter (_ix) != _x) _ix++;
        g->SetPoint (ix, _x, (_y / hden->GetBinContent (_ix)));
        //g->SetPointEYhigh (ix, sqrt (pow (g->GetErrorYhigh (ix), 2) + pow (hden->GetBinError (_ix), 2)) * hden->GetBinWidth (_ix));
        //g->SetPointEYlow (ix, sqrt (pow (g->GetErrorYlow (ix), 2) + pow (hden->GetBinError (_ix), 2)) * hden->GetBinWidth (_ix));
        g->SetPointEYhigh (ix, fabs (_y / hden->GetBinContent (_ix)) * sqrt (pow (g->GetErrorYhigh (ix) / _y, 2) + pow (hden->GetBinError (_ix) / hden->GetBinContent (_ix), 2)));
        g->SetPointEYlow (ix, fabs (_y / hden->GetBinContent (_ix)) * sqrt (pow (g->GetErrorYlow (ix) / _y, 2) + pow (hden->GetBinError (_ix) / hden->GetBinContent (_ix), 2)));
      }

      TH1D* h_s2b = (useTrkPt ? h_trk_pt_ptz_sig_to_bkg : h_trk_xhz_ptz_sig_to_bkg)[iSpc][iPtZ][iCent];
      double x, y;
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y); 
        double s2b = -1, s2b_err = -1;
        for (int iX = 1; s2b == -1 && iX <= h_s2b->GetNbinsX (); iX++) {
          if (h_s2b->GetBinLowEdge (iX) <= x && x <= h_s2b->GetBinLowEdge (iX) + h_s2b->GetBinWidth (iX)) {
            s2b = h_s2b->GetBinContent (iX);
            s2b_err = h_s2b->GetBinError (iX);
          } 
        }
        if (isnan (s2b) || s2b <= 0) continue;
        g->SetPoint (i, s2b, y);
        g->SetPointEXhigh (i, s2b_err);
        g->SetPointEXlow (i, s2b_err);
      }

      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y); 
        if (x <= 0) continue;
        g->SetPoint (i, x, 100*(1.-y)*x);
        g->SetPointEYhigh (i, fabs (100*(1.-y)*x)*sqrt (pow (g->GetErrorYhigh (i) / (1.-y), 2) + pow (g->GetErrorX (i) / x, 2)));
        g->SetPointEYlow (i, fabs (100*(1.-y)*x)*sqrt (pow (g->GetErrorYlow (i) / (1.-y), 2) + pow (g->GetErrorX (i) / x, 2)));
      }

      if (plotAsSyst) {
        SetConstantXErrors (g, 0.032, true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
        g->SetMarkerSize (0); 
        g->SetLineWidth (1);
        g->SetLineColor (finalColors[iPtZ-1]);
        g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

        if (iCent != 0) ((TGAE*) g->Clone ())->Draw ("5P");
      }
      else {
        Style_t markerStyle = markerStyles[iPtZ-2];
        double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.0 : 1.5);

        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (markerSize);
        g->SetMarkerColor (finalColors[iPtZ-1]);
        g->SetLineColor (finalColors[iPtZ-1]);
        g->SetLineWidth (2);

        ((TGAE*) g->Clone ())->Draw ("P");

        if (IsFullMarker (markerStyle)) {
          markerStyle = FullToOpenMarker (markerStyle);

          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (markerSize);
          g->SetLineWidth (0);
          g->SetMarkerColor (kBlack);
  
          ((TGAE*) g->Clone ())->Draw ("P");
        }

        double x, y;
        for (int i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);

          g_allPoints->SetPoint (g_allPoints->GetN (), x, y);
          g_allPoints->SetPointEXlow (g_allPoints->GetN ()-1, g->GetErrorXlow (i));
          g_allPoints->SetPointEXhigh (g_allPoints->GetN ()-1, g->GetErrorXhigh (i));
          g_allPoints->SetPointEYlow (g_allPoints->GetN ()-1, g->GetErrorYlow (i));
          g_allPoints->SetPointEYhigh (g_allPoints->GetN ()-1, g->GetErrorYhigh (i));
        }

      } 
      SaferDelete (&g);
    } // end loop over iPtZ
  } // end loop over iCent

  
  if (!plotAsSyst) {
    TF1* fit = new TF1 (Form ("f_closure_sigToBkg_BkgNorm_%s", useTrkPt ? "pTch" : "xhZ"), "[0]+[1]*log(x)+[2]*log(x)*log(x)", 5e-3, 2e1);
    fit->SetParameter (0, 0);
    fit->SetParameter (1, 0);
    g_allPoints->Fit (fit, "RN0");
    fit->SetLineColor (kBlack);
    fit->Draw ("same");

    if (!canvasExists) {
      myLineText (0.25, 0.71, kBlack, 1, "Fit to [0]+[1]*log(x)+[2]*log^2(x)", 1.4, 0.028);
      myText (0.25, 0.67,  kBlack, Form ("#chi^{2}/dof = %.2f/%i", fit->GetChisquare (), fit->GetNDF ()), 0.028);
    }

    if (string (saveFileName) != "") {
      TFile* saveFile = new TFile (saveFileName, "recreate");
      fit->Write ();
      saveFile->Close ();
      SaferDelete (&saveFile);
    }
  }

  delete g_allPoints;
  
  c->SaveAs (Form ("%s/TrkYields/closureCheck_SigToBkg_BkgNormalization_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Y(pT or xZh) binned in Z Pt
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotAllYields_dPtZ_SpcComp (const bool useTrkPt) {
  const double padRatio = 0.9; // ratio of size of upper pad to middle & lower pads. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (padRatio+1.0);
  const int axisTextSize = 23;

  const char* canvasName = Form ("c_AllYields_%s_dPtZ", useTrkPt ? "pTch" : "xhZ");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 395*numCentBins, 685);
    gDirectory->Add (c);
  }

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    const char* topPadName = Form ("p_AllYields_top_%s_dPtZ_iCent%i", useTrkPt ? "pTch" : "xhZ", iCent);
    const char* bottomPadName = Form ("p_AllYields_bottom_%s_dPtZ_iCent%i", useTrkPt ? "pTch" : "xhZ", iCent);

    TPad* topPad = nullptr, *bottomPad = nullptr;
    if (!canvasExists) {
      topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
      bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

      gDirectory->Add (topPad);
      gDirectory->Add (bottomPad);

      topPad->SetTopMargin (0.04);
      topPad->SetBottomMargin (0);
      topPad->SetLeftMargin (0.20);
      topPad->SetRightMargin (0.06);
      bottomPad->SetTopMargin (0);
      bottomPad->SetBottomMargin (0.20);
      bottomPad->SetLeftMargin (0.20);
      bottomPad->SetRightMargin (0.06);
      topPad->Draw ();
      bottomPad->Draw ();
    }
    else {
      topPad = dynamic_cast<TPad*> (gDirectory->Get (topPadName));
      bottomPad = dynamic_cast<TPad*> (gDirectory->Get (bottomPadName));
    }

    topPad->cd ();
    GetDrawnObjects ();
    gPad->SetLogx ();
    gPad->SetLogy ();

    double min = 1e30, max = 0;
    if (useTrkPt) {
      if (iCent == 0)       { min = 1e-5; max = 4e1; }
      else if (iCent == 1)  { min = 1e-5; max = 4e2; }
      else if (iCent == 2)  { min = 5e-5; max = 2e3; }
      else if (iCent == 3)  { min = 8e-5; max = 6e3; }
    }
    else {
      if (iCent == 0)       { min = 2e-5; max = 2e3; }
      else if (iCent == 1)  { min = 1e-3; max = 5e4; }
      else if (iCent == 2)  { min = 1e-3; max = 9e4; }
      else if (iCent == 3)  { min = 1e-3; max = 9e5; }
    }

    if (!canvasExists) {
      TH1D* htemp = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);

      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      useTrkPt ? xax->SetLimits (trk_min_pt, trk_max_pt) : xax->SetLimits (allXhZBins[0], allXhZBins[maxNXhZBins]);
      yax->SetRangeUser (min, max);

      xax->SetMoreLogLabels ();

      xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
      yax->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T} d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x} d#Delta#phi");

      xax->SetTitleFont (43);
      xax->SetTitleSize (axisTextSize);
      xax->SetLabelFont (43);
      xax->SetLabelSize (axisTextSize);

      yax->SetTitleFont (43);
      yax->SetTitleSize (axisTextSize);
      yax->SetLabelFont (43);
      yax->SetLabelSize (axisTextSize);

      yax->SetTitleOffset (1.5 * yax->GetTitleOffset ());

      htemp->SetLineWidth (0);
      htemp->SetMarkerSize (0);
  
      htemp->DrawCopy ("hist ][");
      SaferDelete (&htemp);

      if (iCent == 0) {
        myText (0.25, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.06);
        myText (0.25, 0.14, kBlack, "Before subtraction", 0.060);
      }
      else myText (0.25, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

      if (iCent == 0)       myText (0.485, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      else if (iCent == 1)  myText (0.485, 0.87, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.06);
      else if (iCent == numCentBins - 1) {
        for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
          if (iPtZ == 2) {
            myText (0.37, 0.87, kBlack, "#it{ee}", 0.06);
            myText (0.47, 0.87, kBlack, "#it{#mu#mu}", 0.06);
          }
          myMarkerTextNoLine (0.42, 0.852-0.06*(iPtZ-2), colors[iPtZ-1], kFullCircle, "", 1.5, 0.054); // for plotting MC reco vs truth
          myMarkerTextNoLine (0.52, 0.852-0.06*(iPtZ-2), colors[iPtZ-1], kOpenCircle, "", 1.5, 0.054);

          if (iPtZ == nPtZBins-1) myText (0.56, 0.85-0.06*(iPtZ-2), kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.054);
          else myText (0.56, 0.85-0.06*(iPtZ-2), kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.054);
        }
      }
    }

    for (short iSpc = 0; iSpc < 2; iSpc++) {
      for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        const Style_t markerStyle = (iSpc == 1 ? kOpenCircle : kFullCircle);

        TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[iSpc][iPtZ][iCent]);
        RecenterGraph (g);

        if (!plotAsSyst) {
          ResetXErrors (g);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerColor (colors[iPtZ-1]);
          g->SetLineColor (colors[iPtZ-1]);
          g->SetMarkerSize (1);
          g->SetLineWidth (2);
        } else {
          g->SetMarkerSize (0); 
          g->SetLineWidth (1);
          g->SetLineColor (colors[iPtZ-1]);
          g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
        }

        if (!plotAsSyst) g->Draw ("P");
        else {
          ((TGAE*)g->Clone ())->Draw ("5P");
          g->Draw ("2P");
        }
      } // end loop over iPtZ
    } // end loop over iSpc

    bottomPad->cd ();
    GetDrawnObjects ();
    //plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
    gPad->SetLogx ();

    if (!canvasExists) {
      TH1D* htemp = new TH1D ("htemp", "", (useTrkPt ? nPtchBins : nXhZBins)[nPtZBins-1], (useTrkPt ? pTchBins : xhZBins)[nPtZBins-1]);

      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      useTrkPt ? xax->SetLimits (trk_min_pt, trk_max_pt) : xax->SetLimits (allXhZBins[0], allXhZBins[maxNXhZBins]);
      yax->SetRangeUser (0, 2);

      xax->SetMoreLogLabels ();

      xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
      yax->SetTitle ("Electrons / Muons");

      xax->SetTitleFont (43);
      xax->SetTitleSize (axisTextSize);
      xax->SetLabelFont (43);
      xax->SetLabelSize (axisTextSize);

      yax->SetTitleFont (43);
      yax->SetTitleSize (axisTextSize);
      yax->SetLabelFont (43);
      yax->SetLabelSize (axisTextSize);

      xax->SetTitleOffset (1.8 * xax->GetTitleOffset ());
      yax->SetTitleOffset (1.5 * yax->GetTitleOffset ());

      htemp->SetLineWidth (0);
      htemp->SetMarkerSize (0);
  
      htemp->DrawCopy ("hist ][");
      SaferDelete (&htemp);
    }

    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

      TH1D* h = (useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[0][iPtZ][iCent];
      h->Divide ((useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[1][iPtZ][iCent]);

      TGAE* g = GetTGAE (h);
      SaferDelete (&h);
      RecenterGraph (g);

      if (!plotAsSyst) {
        ResetXErrors (g);
        g->SetMarkerColor (colors[iPtZ-1]);
        g->SetLineColor (colors[iPtZ-1]);
        g->SetMarkerSize (1);
        g->SetLineWidth (2);
      } else {
        g->SetMarkerSize (0); 
        g->SetLineWidth (1);
        g->SetLineColor (colors[iPtZ-1]);
        g->SetFillColorAlpha (fillColors[iPtZ-1], 0.3);
      }

      if (!plotAsSyst) g->Draw ("P");
      else {
        ((TGAE*)g->Clone ())->Draw ("5P");
        g->Draw ("2P");
      }
    } // end loop over iPtZ
  } // end loop over iCent
  
  c->SaveAs (Form ("%s/TrkYields/allYields_%s_SpcComp.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots hadron yield in pp for ee and mumu separately
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotPPYields_dPtZ_SpcComp (const bool useTrkPt, const short iPtZ) {

  const char* canvasName = Form ("c_pp_yields_channelCompare_%s_iPtZ%i", useTrkPt ? "pTch" : "xhZ", iPtZ);
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  TPad* uPad = nullptr, *dPad = nullptr;
  if (canvasExists) {
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
    dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
  }
  else {
    c = new TCanvas (canvasName, "", 800, 800);
    c->cd ();
    uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
    dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
    uPad->SetBottomMargin (0);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    uPad->Draw ();
    dPad->Draw ();
    gDirectory->Add (c);
    gDirectory->Add (uPad);
    gDirectory->Add (dPad);
  }
  c->cd ();

  uPad->cd ();
  gPad->SetLogx ();
  gPad->SetLogy ();

  if (!canvasExists) {
    TH1D* htemp = new TH1D ("htemp", "", 1, (useTrkPt ? pTchBins : xhZBins)[iPtZ][0], (useTrkPt ? pTchBins : xhZBins)[iPtZ][(useTrkPt ? nPtchBins : nXhZBins)[iPtZ]]);
    htemp->Reset ();
    htemp->GetYaxis ()->SetRangeUser (useTrkPt ? 1e-2 : 8e-2, useTrkPt ? 1e1 : 8e1);
    htemp->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
    htemp->GetYaxis ()->SetTitle (Form ("(1/N_{Z}) (d^{2}N_{ch}/d%sd#Delta#phi)%s", useTrkPt ? "#it{p}_{T}" : "#it{x}_{hZ}", useTrkPt ? " [GeV^{-1}]" : ""));
    htemp->GetXaxis ()->SetTitleSize (0);
    htemp->GetYaxis ()->SetTitleSize (0.04/0.6);
    htemp->GetXaxis ()->SetLabelSize (0);
    htemp->GetYaxis ()->SetLabelSize (0.04/0.6);
    htemp->GetYaxis ()->SetTitleOffset (1.1);
    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);
  }

  for (short iSpc : {0, 1}) {
    TH1D* h = (useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[iSpc][iPtZ][0];

    TGraphAsymmErrors* g = GetTGAE (h);
    ResetXErrors (g);

    Style_t markerStyle = kFullCircle;
    const float markerSize = 1.25;
    Color_t markerColor = (iSpc == 0 ? finalColors[3] : finalColors[1]);

    g->SetMarkerStyle (markerStyle);
    g->SetMarkerColor (markerColor);
    g->SetMarkerSize (markerSize);
    g->SetLineWidth (2);
    g->SetLineColor (markerColor);
    ((TGAE*) g->Clone ())->Draw ("P");

    markerStyle = FullToOpenMarker (markerStyle);

    g->SetMarkerStyle (markerStyle);
    g->SetMarkerSize (markerSize);
    g->SetMarkerColor (kBlack);

    ((TGAE*) g->Clone ())->Draw ("P");

    SaferDelete (&g);
  }

  dPad->cd ();
  gPad->SetLogx ();

  {
    TH1D* htemp = new TH1D ("htemp", "", 1, (useTrkPt ? pTchBins : xhZBins)[iPtZ][0], (useTrkPt ? pTchBins : xhZBins)[iPtZ][(useTrkPt ? nPtchBins : nXhZBins)[iPtZ]]);
    htemp->GetXaxis ()->SetMoreLogLabels ();
    htemp->GetYaxis ()->SetRangeUser (-0.05, 0.05);
    htemp->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
    htemp->GetYaxis ()->SetTitle ("(Y_{#it{ll}} #minus Y_{comb.}) / Y_{comb.}");
    htemp->GetXaxis ()->SetTitleSize (0.04/0.4);
    htemp->GetYaxis ()->SetTitleSize (0.04/0.4);
    htemp->GetXaxis ()->SetLabelSize (0.04/0.4);
    htemp->GetYaxis ()->SetLabelSize (0.04/0.4);
    htemp->GetXaxis ()->SetTitleOffset (1.2);
    htemp->GetYaxis ()->SetTitleOffset (1.1*0.4/0.6);
    htemp->GetYaxis ()->CenterTitle ();
    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);

    TLine* l = new TLine (0, 1, 300, 1);
    l->SetLineColor (46);
    l->SetLineWidth (2);
    l->SetLineStyle (5);
    l->Draw ("same");
  }

  for (short iSpc : {0, 1}) {
    TH1D* h = (TH1D*) (useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[iSpc][iPtZ][0]->Clone ("h_ratio");
    h->Add ((useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[2][iPtZ][0], -1);
    h->Divide ((useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[2][iPtZ][0]);

    TGraphAsymmErrors* g = GetTGAE (h);
    ResetXErrors (g);

    Style_t markerStyle = kFullCircle;
    const float markerSize = 1.25;
    Color_t markerColor = (iSpc == 0 ? finalColors[3] : finalColors[1]);

    g->SetMarkerStyle (markerStyle);
    g->SetMarkerColor (markerColor);
    g->SetMarkerSize (markerSize);
    g->SetLineWidth (2);
    g->SetLineColor (markerColor);
    ((TGAE*) g->Clone ())->Draw ("P");

    markerStyle = FullToOpenMarker (markerStyle);

    g->SetMarkerStyle (markerStyle);
    g->SetMarkerSize (markerSize);
    g->SetMarkerColor (kBlack);

    ((TGAE*) g->Clone ())->Draw ("P");

    SaferDelete (&g);
  }

  uPad->cd ();

  if (isMC) myText (0.44, 0.85, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.045/0.6);
  else      myText (0.62, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);
  myText (0.62, 0.75, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
  if (iPtZ == nPtZBins-1) myText (0.62, 0.67, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04/0.6);
  else                    myText (0.62, 0.67, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04/0.6);
  myMarkerText (0.75, 0.60, finalColors[3], kFullCircle, "#it{Z} #rightarrow #it{ee}", 1.25, 0.04/0.6, true);
  myMarkerText (0.75, 0.53, finalColors[1], kFullCircle, "#it{Z} #rightarrow #it{#mu#mu}", 1.25, 0.04/0.6, true);


  c->SaveAs (Form ("%s/TrkYields/pp_yields_%s_iPtZ%i_channelCompare.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iPtZ));
  

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots all track yield covariance matrices.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotAllCovMatrices () {
  for (short iSpc : {0, 1})
    for (bool iPt : {true, false})
      for (short iCent : {0, 1, 2, 3})
        for (short iPtZ : {2, 3, 4})
          PlotCovMatrix (iPt, iSpc, iPtZ, iCent);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the track yield covariance matrix.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotCovMatrix (const bool useTrkPt, const short pSpc, const short pPtZ, const short pCent) {
  const char* canvasName = Form ("c_CovMatrix_%s_iSpc%i_iPtZ%i_iCent%i", useTrkPt ? "pTch" : "xhZ", pSpc, pPtZ, pCent);
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    FormatTH2Canvas (c, true);
    gDirectory->Add (c);
  }

  c->cd ();

  gPad->SetLogx ();
  gPad->SetLogy ();
  //gPad->SetLogz ();

  TH2D* h2 = (useTrkPt ? h2_trk_pt_ptz_cov : h2_trk_xhz_ptz_cov)[pSpc][pPtZ][pCent];
  TH2D* h2temp = (TH2D*) h2->Clone ("h2temp");
  //TH1D* h1 = (useTrkPt ? h_trk_pt_ptz : h_trk_xhz_ptz)[pSpc][pPtZ][pCent];

  for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
    for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
      //h2->SetBinContent (iX, iY, sqrt (h2->GetBinContent (iX, iY) / ((h1->GetBinContent (iX)) * (h1->GetBinContent (iY)))));
      h2->SetBinContent (iX, iY, h2temp->GetBinContent (iX, iY) / (sqrt (h2temp->GetBinContent (iX, iX)) * sqrt (h2temp->GetBinContent (iY, iY))));
    }
  }
  SaferDelete (&h2temp);

  h2->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
  h2->GetYaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");

  h2->GetXaxis ()->SetMoreLogLabels ();
  h2->GetYaxis ()->SetMoreLogLabels ();

  //h2->GetZaxis ()->SetTitle (useTrkPt ? "[cov (Y(x), Y(y)) / N_{evt} #LTY(x)#GT#LTY(y)#GT]^{1/2}" : "[cov (Y(x), Y(y)) / N_{evt} #LTY(x)#GT#LTY(y)#GT]^{1/2}");
  h2->GetZaxis ()->SetTitle ("#rho #equiv #sigma_{xy} / #sigma_{x}#sigma_{y}");
  //h2->GetZaxis ()->SetTitle (useTrkPt ? "cov (Y(x), Y(y)) / N_{evt} [GeV^{-2}]" : "cov (Y(x), Y(y)) / N_{evt}");

  h2->GetZaxis ()->SetTitleOffset (1.2 * h2->GetZaxis ()->GetTitleOffset ());

  h2->GetZaxis ()->SetRangeUser (-1, 1); // correlation coefficient is between -1 and 1
  //h2->GetZaxis ()->SetRangeUser (0.5 * h2->GetMinimum (0), 2 * h2->GetMaximum ()); // log z scale
  //h2->GetZaxis ()->SetRangeUser (1.2 * h2->GetMinimum (), 1.2 * h2->GetMaximum ()); // linear z scale
  h2->Draw ("colz");

  myText (0.23, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  if (pCent == 0) myText (0.23, 0.80, kBlack, "#it{pp}, 5.02 TeV", 0.04);
  else            myText (0.23, 0.80, kBlack, Form ("Pb+Pb, 5.02 TeV, %i-%i%%", (int)centCuts[pCent], (int)centCuts[pCent-1]), 0.04);
  if (pPtZ == nPtZBins-1) myText (0.23, 0.75, kBlack, Form ("#it{Z} #rightarrow %s, #it{p}_{T}^{Z} > %g GeV", (pSpc == 0 ? "#it{e}#it{e}" : (pSpc == 1 ? "#it{#mu}#it{#mu}" : "#it{l}#it{l}")), zPtBins[pPtZ]), 0.04);
  else                    myText (0.23, 0.75, kBlack, Form ("#it{Z} #rightarrow %s, %g < #it{p}_{T}^{Z} < %g GeV", (pSpc == 0 ? "#it{e}#it{e}" : (pSpc == 1 ? "#it{#mu}#it{#mu}" : "#it{l}#it{l}")), zPtBins[pPtZ], zPtBins[pPtZ+1]), 0.04);

  c->SaveAs (Form ("%s/CovMatrices/covMatrix_%s_iSpc%i_iPtZ%i_iCent%i.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", pSpc, pPtZ, pCent));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: CalculateIAA (const bool overwrite) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (iaaCalculated && !overwrite)
    return;
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      TH1D* ppHist = nullptr, *PbPbHist = nullptr;

      ppHist = h_trk_pt_ptz_sub[iSpc][iPtZ][0];
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        SaferDelete (&(h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]));
        PbPbHist = (TH1D*)(h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_pt_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())));
        PbPbHist->Divide (ppHist);
        h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent] = PbPbHist;
      } // end loop over iCent

      ppHist = h_trk_xhz_ptz_sub[iSpc][iPtZ][0];
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        SaferDelete (&(h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]));
        PbPbHist = (TH1D*)(h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_trk_xhz_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())));
        PbPbHist->Divide (ppHist);
        h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent] = PbPbHist;
      } // end loop over iCent

      //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
      //  TH1D* ppHist = h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][0];

      //  for (short iCent = 1; iCent < numCentBins; iCent++) {
      //    SaferDelete (&(h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent]));
      //    TH1D* PbPbHist = (TH1D*)(h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_pt_dphi_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())));
      //    PbPbHist->Divide (ppHist);
      //    h_trk_pt_dphi_iaa[iSpc][iPtZ][iPhi][iCent] = PbPbHist;
      //  } // end loop over iCent
      //  ppHist = h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][0];

      //  for (short iCent = 1; iCent < numCentBins; iCent++) {
      //    SaferDelete (&(h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent]));
      //    TH1D* PbPbHist = (TH1D*)(h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_trk_xhz_dphi_iaa_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ())));
      //    PbPbHist->Divide (ppHist);
      //    h_trk_xhz_dphi_iaa[iSpc][iPtZ][iPhi][iCent] = PbPbHist;
      //  } // end loop over iCent
      //} // end loop over iPhi
    } // end loop over iPtZ
  } // end loop over iSpc
  iaaCalculated = true;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots IAA for each deltaPhi bin in all centralities
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAA_dPhi (const bool useTrkPt, const short pSpc, const short pPtZ) {
  if (!iaaCalculated)
    CalculateIAA ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_iaa_%s_dPhi_iPtZ%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numCentBins-1), 502);
        c->Divide (numCentBins-1, 1);
        gDirectory->Add (c);
      }
      c->cd ();

      double xmin = 0, xmax = 0;
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);

        if (!canvasExists) {
          gPad->SetLogx ();
  
          TH1D* htemp = new TH1D ("htemp", "", useTrkPt ? nPtchBins[nPtZBins-1] : nXhZBins[nPtZBins-1], useTrkPt ? pTchBins[nPtZBins-1] : xhZBins[nPtZBins-1]);

          TAxis* xax = htemp->GetXaxis ();
          TAxis* yax = htemp->GetYaxis ();
  
          xax->SetMoreLogLabels ();
  
          xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
          yax->SetTitle ("I_{AA}");
          xmin = xax->GetXmin ();
          xmax = xax->GetXmax ();
  
          yax->SetRangeUser (0, max_iaa);
  
          xax->SetTitleFont (43);
          xax->SetTitleSize (axisTextSize);
          xax->SetLabelFont (43);
          xax->SetLabelSize (axisTextSize);
  
          yax->SetTitleFont (43);
          yax->SetTitleSize (axisTextSize);
          yax->SetLabelFont (43);
          yax->SetLabelSize (axisTextSize);
  
          xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
          //yax->SetTitleOffset (0.9 * yax->GetTitleOffset ());

          htemp->SetLineWidth (0);
          htemp->SetMarkerSize (0);
  
          htemp->DrawCopy ("hist ][");
          SaferDelete (&htemp);

          if (iCent == 1)
            myText (0.50, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
          else if (iCent == 2) {
            if (iPtZ == nPtZBins-1) myText (0.50, 0.85, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
            else myText (0.50, 0.85, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
          }
          else if (iCent == numCentBins-1) {
            myText (0.625, 0.88, kBlack, "#Delta#phi", 0.05);
            for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
              const char* lo = GetPiString (phiLowBins[iPhi]);
              const char* hi = GetPiString (phiHighBins[iPhi]);
              //myMarkerTextNoLine (0.60, 0.94-0.06*(iPhi+1), colors[iPhi], kFullCircle, Form ("(%s, %s)", lo, hi), 2.0, 0.05);
              //myMarkerTextNoLine (0.60, 0.94-0.06*(iPhi+1), kBlack, kOpenCircle, "", 2.0, 0.05);
              myMarkerAndBoxAndLineText (0.60, 0.94-0.06*(iPhi+1), 2.0, 1001, fillColors[iPhi], 0.30, colors[iPhi], kFullCircle, 1.8, Form ("(%s, %s)", lo, hi), 0.05);
              //myText (0.62, 0.91-0.06*(iPhi+1), kBlack, Form ("(%s, %s)", lo, hi), 0.05);
            }
          }
          myText (0.22, 0.20, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

          TLine* l = new TLine (xmin, 1, xmax, 1);
          l->SetLineStyle (2);
          l->SetLineWidth (2);
          l->SetLineColor (kPink-8);
          l->Draw ("same");
        }

        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi_iaa : h_trk_xhz_dphi_iaa)[iSpc][iPtZ][iPhi][iCent]);
          RecenterGraph (g);
          ResetXErrors (g);
          deltaize (g, 0.09*(iPhi-1 - 0.5*(numPhiBins-2)), true, useTrkPt ? pTchBins[nPtZBins-1][0] : xhZBins[nPtZBins-1][0], useTrkPt ? pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]] : xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]]);
          ResetXErrors (g);

          if (plotAsSyst) {
            SetConstantXErrors (g, 0.040, true, useTrkPt ? pTchBins[nPtZBins-1][0] : xhZBins[nPtZBins-1][0], useTrkPt ? pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]] : xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]]);
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (colors[iPhi]);
            g->SetFillColorAlpha (fillColors[iPhi], 0.3);

            ((TGAE*) g->Clone ())->Draw ("5P");
          }
          else {
            Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);
            double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

            g->SetMarkerStyle (markerStyle);
            g->SetMarkerSize (markerSize);
            g->SetMarkerColor (colors[iPhi]);
            g->SetLineColor (colors[iPhi]);
            g->SetLineWidth (2);

            ((TGAE*) g->Clone ())->Draw ("P");

            markerStyle = FullToOpenMarker (markerStyle);

            g->SetMarkerStyle (markerStyle);
            g->SetMarkerSize (markerSize);
            g->SetLineWidth (0);
            g->SetMarkerColor (kBlack);
  
            ((TGAE*) g->Clone ())->Draw ("P");
          }

          SaferDelete (&g);

        } // end loop over iPhi
      } // end loop over iCent
      c->cd (1);

      c->SaveAs (Form ("%s/IAA/iaa_%s_dPhi_iPtZ%i_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots IAA for each centrality in all deltaPhi bins
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAA_dCent (const bool useTrkPt, const short pSpc, const short pPtZ) {
  if (!iaaCalculated)
    CalculateIAA ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_iaa_%s_dCent_iPtZ%i_%s", useTrkPt ? "pTch" : "xhZ", iPtZ, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numPhiBins-1), 500);
        c->Divide (numPhiBins-1, 1);
        gDirectory->Add (c);
      }
      c->cd ();

      double xmin = 0, xmax = 0;
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        c->cd (iPhi);
        gPad->SetLogx ();

        for (int iCent = 1; iCent < numCentBins; iCent++) {
          const Style_t markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);

          TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_dphi_iaa : h_trk_xhz_dphi_iaa)[iSpc][iPtZ][iPhi][iCent]);
          RecenterGraph (g);

          if (!plotAsSyst) {
            ResetXErrors (g);
            deltaize (g, 1+(iCent-2)*0.02, true); // 2.5 = 0.5*(numPhiBins-1)
            g->SetLineColor (colors[iCent]);
            g->SetMarkerColor (colors[iCent]);
            g->SetMarkerStyle (markerStyle);
            g->SetMarkerSize (1.2);
            g->SetLineWidth (2);
          } else {
            g->SetMarkerSize (0); 
            g->SetLineWidth (1);
            g->SetLineColor (colors[iCent]);
            g->SetFillColorAlpha (fillColors[iCent], 0.3);
          }

          useTrkPt ? g->GetXaxis ()->SetLimits (trk_min_pt, pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[nPtZBins-1]]);
          g->GetYaxis ()->SetRangeUser (0, max_iaa);
          xmin = g->GetXaxis ()->GetXmin ();
          xmax = g->GetXaxis ()->GetXmax ();

          g->GetXaxis ()->SetMoreLogLabels ();

          g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
          g->GetYaxis ()->SetTitle ("I_{AA}");

          g->GetXaxis ()->SetTitleFont (43);
          g->GetXaxis ()->SetTitleSize (axisTextSize);
          g->GetXaxis ()->SetLabelFont (43);
          g->GetXaxis ()->SetLabelSize (axisTextSize);

          g->GetYaxis ()->SetTitleFont (43);
          g->GetYaxis ()->SetTitleSize (axisTextSize);
          g->GetYaxis ()->SetLabelFont (43);
          g->GetYaxis ()->SetLabelSize (axisTextSize);

          g->GetXaxis ()->SetTitleOffset (0.9 * g->GetXaxis ()->GetTitleOffset ());
          //g->GetYaxis ()->SetTitleOffset (0.9 * g->GetYaxis ()->GetTitleOffset ());

          if (!plotAsSyst) {
            string drawString = string (!canvasExists && iCent == 1 ? "AP" : "P");
            g->Draw (drawString.c_str ());
          } else {
            string drawString = string (!canvasExists && iCent == 1 ? "A5P" : "5P");
            ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
            g->Draw ("2P");
          }

          if (!canvasExists) LabelIAA_dCent (iCent, iPhi, iPtZ);
        } // end loop over iPhi
      } // end loop over iCent
      c->cd (1);

      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        c->cd (iPhi);
        const char* lo = GetPiString (phiLowBins[iPhi]);
        const char* hi = GetPiString (phiHighBins[iPhi]);
        myText (0.22, 0.24, kBlack, Form ("%s < #Delta#phi < %s", lo, hi), 0.06);

        TLine* l = new TLine (xmin, 1, xmax, 1);
        l->SetLineStyle (2);
        l->SetLineWidth (2);
        l->SetLineColor (kPink-8);
        l->Draw ("same");
      } // end loop over iCent

      c->SaveAs (Form ("%s/IAA/iaa_%s_dCent_iPtZ%i_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iPtZ, spc));
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary (non-virtual) plot labelling for IAA measurements
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: LabelIAA_dCent (const short iCent, const short iPhi, const short iPtZ) {

  if (iPhi == 1 && iCent == 1) {
    myText (0.50, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
    if (iPtZ == nPtZBins-1) myText (0.50, 0.78, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
    else myText (0.50, 0.78, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
  }
  else if (iPhi == numPhiBins-1) {
    //TVirtualPad* cPad = gPad; // store current pad
    //TBox* b = TBoxNDC (0.61-0.024, 0.91-0.06*iPhi-0.016, 0.61+0.024, 0.91-0.06*iPhi+0.016);
    //b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
    //b->Draw ("l");
    //cPad->cd ();
    //myMarkerTextNoLine (0.61, 0.912-0.06*iPhi, colors[iPhi], kOpenCircle, "", 1.4, 0.05);
    //myMarkerTextNoLine (0.512, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
    myMarkerTextNoLine (0.63, 0.925-0.06*iCent, colors[iCent], kFullCircle, "", 1.4, 0.05);
    myText (0.65, 0.91-0.06*iCent, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.05);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots IAA for each pT^Z bin in all centralities
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAA_dPtZ (const bool useTrkPt, const short pSpc) {
  if (!iaaCalculated)
    CalculateIAA ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_iaa_%s_dPtZ_%s", useTrkPt ? "pTch" : "xhZ", spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 500*(numCentBins-1), 502);
      c->Divide (numCentBins-1, 1);
      gDirectory->Add (c);
    }
    c->cd ();

    double xmin = 0, xmax = 0;
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      c->cd (iCent);

      if (!canvasExists) {
        gPad->SetLogx ();
        gPad->SetLogy ();

        TH1D* htemp = new TH1D ("htemp", "", useTrkPt ? nPtchBins[nPtZBins-1] : nXhZBins[nPtZBins-1], useTrkPt ? pTchBins[nPtZBins-1] : xhZBins[nPtZBins-1]);

        TAxis* xax = htemp->GetXaxis ();
        TAxis* yax = htemp->GetYaxis ();

        xax->SetMoreLogLabels ();

        xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
        yax->SetTitle ("I_{AA}");
        xmin = xax->GetXmin ();
        xmax = xax->GetXmax ();

        yax->SetRangeUser (0.05, 10);
        //yax->SetRangeUser (0, max_iaa);

        xax->SetTitleFont (43);
        xax->SetTitleSize (axisTextSize);
        xax->SetLabelFont (43);
        xax->SetLabelSize (axisTextSize);

        yax->SetTitleFont (43);
        yax->SetTitleSize (axisTextSize);
        yax->SetLabelFont (43);
        yax->SetLabelSize (axisTextSize);

        xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
        //yax->SetTitleOffset (0.9 * yax->GetTitleOffset ());

        htemp->SetLineWidth (0);
        htemp->SetMarkerSize (0);
  
        htemp->DrawCopy ("hist ][");
        SaferDelete (&htemp);

        myText (0.22, 0.20, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);

        if (iCent == 1) {
          myText (0.26, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.065);
          myText (0.26, 0.345, kBlack, "Pb+Pb, 5.02 TeV, 1.4-1.7 nb^{-1}", 0.05);
          myText (0.26, 0.28, kBlack, "#it{pp}, 5.02 TeV, 260 pb^{-1}", 0.05);
        }
        else if (iCent == 2) {
          const char* lo = GetPiString (phiLowBins[1]);
          const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
          myText (0.50, 0.86, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.05);
        }
        else if (iCent == numCentBins-1) {
          for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
            const Style_t markerStyle = markerStyles[iPtZ-2];
            const double markerSize = (markerStyle == kOpenDiamond || markerStyle == kFullDiamond ? 2.0 : 1.3);
            if (iPtZ == nPtZBins-1)
              myMarkerAndBoxAndLineText (0.56, 0.88-0.07*(iPtZ-2), 2.0, 1001, finalFillColors[iPtZ-1], 0.30, finalColors[iPtZ-1], markerStyle, markerSize, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.05);
            else
              myMarkerAndBoxAndLineText (0.56, 0.88-0.07*(iPtZ-2), 2.0, 1001, finalFillColors[iPtZ-1], 0.30, finalColors[iPtZ-1], markerStyle, markerSize, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
          } // end loop over iPtZ
        }

        TLine* l = new TLine (xmin, 1, xmax, 1);
        l->SetLineStyle (2);
        l->SetLineWidth (2);
        l->SetLineColor (kPink-8);
        l->Draw ("same");
      }

      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_iaa : h_trk_xhz_ptz_iaa)[iSpc][iPtZ][iCent]);
        RecenterGraph (g);
        ResetXErrors (g);
        deltaize (g, 0.09*(iPtZ-2 - 0.5*(nPtZBins-3)), true, useTrkPt ? pTchBins[nPtZBins-1][0] : xhZBins[nPtZBins-1][0], useTrkPt ? pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]] : xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]]);
        ResetXErrors (g);

        if (plotAsSyst) {
          SetConstantXErrors (g, 0.040, true, useTrkPt ? pTchBins[nPtZBins-1][0] : xhZBins[nPtZBins-1][0], useTrkPt ? pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]] : xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]]);
          g->SetMarkerSize (0); 
          g->SetLineWidth (1);
          g->SetLineColor (finalColors[iPtZ-1]);
          g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

          ((TGAE*) g->Clone ())->Draw ("5P");
        }
        else {
          Style_t markerStyle = markerStyles[iPtZ-2];
          double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (markerSize);
          g->SetMarkerColor (finalColors[iPtZ-1]);
          g->SetLineColor (finalColors[iPtZ-1]);
          g->SetLineWidth (2);

          ((TGAE*) g->Clone ())->Draw ("P");

          markerStyle = FullToOpenMarker (markerStyle);

          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (markerSize);
          g->SetLineWidth (0);
          g->SetMarkerColor (kBlack);
  
          ((TGAE*) g->Clone ())->Draw ("P");
        }
      } // end loop over iPtZ
    } // end loop over iCent

    c->SaveAs (Form ("%s/IAA/iaa_%s_dPtZ_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", spc));
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots IAA for each pT^Z bin in a single centrality
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotSingleIAA_dPtZ (const bool useTrkPt, const short pPtZ, const short iCent, const short pSpc) {
  if (!iaaCalculated)
    CalculateIAA ();

  short iPtZLo = 2;
  short iPtZHi = nPtZBins;
  if (pPtZ != -1) {
    iPtZLo = pPtZ;
    iPtZHi = pPtZ+1;
  }

  const int axisTextSize = 35;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_singlePlot_iaa_%s_iPtZ%i_iCent%i_%s", useTrkPt ? "pTch" : "xhZ", pPtZ, iCent, spc);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 800, 800);
      gDirectory->Add (c);
    }
    c->cd ();

    double xmin = 0, xmax = 0;
    gPad->SetLogx ();
    gPad->SetLogy ();

    if (!canvasExists) {
      TH1D* htemp = new TH1D ("htemp", "", useTrkPt ? nPtchBins[nPtZBins-1] : nXhZBins[nPtZBins-1], useTrkPt ? pTchBins[nPtZBins-1] : xhZBins[nPtZBins-1]);

      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetMoreLogLabels ();
      yax->SetMoreLogLabels ();

      xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
      yax->SetTitle ("I_{AA}");
      xmin = xax->GetXmin ();
      xmax = xax->GetXmax ();

      yax->SetRangeUser (0.05, 10);

      xax->SetTitleFont (43);
      xax->SetTitleSize (axisTextSize);
      xax->SetLabelFont (43);
      xax->SetLabelSize (axisTextSize);

      yax->SetTitleFont (43);
      yax->SetTitleSize (axisTextSize);
      yax->SetLabelFont (43);
      yax->SetLabelSize (axisTextSize);

      xax->SetTitleOffset (1.2 * xax->GetTitleOffset ());
      //yax->SetTitleOffset (0.9 * yax->GetTitleOffset ());

      htemp->SetLineWidth (0);
      htemp->SetMarkerSize (0);
  
      htemp->DrawCopy ("hist ][");
      SaferDelete (&htemp);

      myText (0.20, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.044);
      myText (0.20, 0.26, kBlack, "Pb+Pb, 5.02 TeV, 1.4-1.7 nb^{-1}", 0.036);
      myText (0.20, 0.215, kBlack, "#it{pp}, 5.02 TeV, 260 pb^{-1}", 0.036);
      myText (0.76, 0.215, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.044);
      //myText (0.20, 0.77, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.036);

      myMarkerAndBoxAndLineText (0.68, 0.880, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "#it{p}_{T}^{Z} = 15-30 GeV", 0.036);
      myMarkerAndBoxAndLineText (0.68, 0.820, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "#it{p}_{T}^{Z} = 30-60 GeV", 0.036);
      myMarkerAndBoxAndLineText (0.68, 0.760, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "#it{p}_{T}^{Z} = 60+ GeV", 0.036);

      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    }

    for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
      Style_t markerStyle = markerStyles[iPtZ-iPtZLo];
      double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);

      TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_iaa : h_trk_xhz_ptz_iaa)[iSpc][iPtZ][iCent]);
      RecenterGraph (g);
      ResetXErrors (g);
      deltaize (g, 0.09*(iPtZ-iPtZLo - 0.5*(iPtZHi-iPtZLo)), true, useTrkPt ? pTchBins[nPtZBins-1][0] : xhZBins[nPtZBins-1][0], useTrkPt ? pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]] : xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]]);
      ResetXErrors (g);

      if (plotAsSyst) {
        SetConstantXErrors (g, 0.04, true, useTrkPt ? pTchBins[nPtZBins-1][0] : xhZBins[nPtZBins-1][0], useTrkPt ? pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]] : xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]]);
        g->SetMarkerSize (0); 
        g->SetLineWidth (1);
        g->SetLineColor (finalColors[iPtZ-1]);
        g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

        ((TGAE*) g->Clone ())->Draw ("5P");
      }
      else {
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (markerSize);
        g->SetMarkerColor (finalColors[iPtZ-1]);
        g->SetLineColor (finalColors[iPtZ-1]);
        g->SetLineWidth (2);

        ((TGAE*) g->Clone ())->Draw ("P");

        markerStyle = FullToOpenMarker (markerStyle);

        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (markerSize);
        g->SetLineWidth (0);
        g->SetMarkerColor (kBlack);
  
        ((TGAE*) g->Clone ())->Draw ("P");
      }

      SaferDelete (&g);
    } // end loop over iPtZ
    c->SaveAs (Form ("%s/IAA/iaa_%s_iCent%i_dPtZ_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iCent, spc));
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between some Pb+Pb centrality and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotIAA_dPtZ_SpcComp (const bool useTrkPt) {
  if (!iaaCalculated)
    CalculateIAA ();

  short iPtZLo = 2;
  short iPtZHi = nPtZBins;
  short iCentLo = 1;
  short iCentHi = numCentBins;

  const int axisTextSize = 22;

  const char* canvasName = Form ("c_iaa_%s_SpcComp", useTrkPt ? "pTch" : "xhZ");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 400*(numCentBins-1), 1200);
    c->Divide (numCentBins-1, iPtZHi - iPtZLo);
    gDirectory->Add (c);
  }
  c->cd ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {

    for (short iCent = iCentLo; iCent < iCentHi; iCent++) {
      for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
        const Style_t markerStyle = kFullCircle;
        c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo) + iCent-iCentLo + 1);

        gPad->SetLogx ();

        TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_iaa : h_trk_xhz_ptz_iaa)[iSpc][iPtZ][iCent]);
        RecenterGraph (g);

        if (!plotAsSyst) {
          ResetXErrors (g);
          deltaize (g, 1+iSpc*0.05-0.025, true); // 2.5 = 0.5*(numPhiBins-1)
          g->SetLineColor (colors[iSpc+1]);
          g->SetMarkerColor (colors[iSpc+1]);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (3);
        } else {
          g->SetMarkerSize (0); 
          g->SetLineWidth (1);
          g->SetLineColor (colors[iSpc+1]);
          g->SetFillColorAlpha (fillColors[iSpc+1], 0.3);
        }

        useTrkPt ? g->GetXaxis ()->SetLimits (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]) : g->GetXaxis ()->SetLimits (xhZBins[nPtZBins-1][0], xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]]);
        g->GetYaxis ()->SetRangeUser (0, max_iaa);

        g->GetXaxis ()->SetMoreLogLabels ();
        g->GetYaxis ()->SetMoreLogLabels ();

        g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
        g->GetYaxis ()->SetTitle ("I_{AA}");

        g->GetXaxis ()->SetTitleFont (43);
        g->GetXaxis ()->SetTitleSize (axisTextSize);
        g->GetXaxis ()->SetLabelFont (43);
        g->GetXaxis ()->SetLabelSize (axisTextSize);

        g->GetYaxis ()->SetTitleFont (43);
        g->GetYaxis ()->SetTitleSize (axisTextSize);
        g->GetYaxis ()->SetLabelFont (43);
        g->GetYaxis ()->SetLabelSize (axisTextSize);

        g->GetXaxis ()->SetTitleOffset (3.0 * g->GetXaxis ()->GetTitleOffset ());
        g->GetYaxis ()->SetTitleOffset (3.0 * g->GetYaxis ()->GetTitleOffset ());

        if (!plotAsSyst) {
          string drawString = string (!canvasExists && iSpc == 0 ? "AP" : "P");
          g->Draw (drawString.c_str ());
        } else {
          string drawString = string (!canvasExists && iSpc == 0 ? "A5P" : "5P");
          ((TGAE*)g->Clone ())->Draw (drawString.c_str ());
          g->Draw ("2P");
        }
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc

  if (!canvasExists) {
    for (short iCent = iCentLo; iCent < iCentHi; iCent++) {
      for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
        c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo) + iCent-iCentLo + 1);
        if (iPtZ == iPtZHi-1)
          myText (0.65, 0.77, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.08);
        if (iCent == iCentLo && iPtZ == iPtZLo) {
          myText (0.25, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.08);
          myText (0.25, 0.75, kBlack, "Pb+Pb, 5.02 TeV", 0.07);
        }
        //if (iCent == iCentLo+1 && iPtZ == iPtZLo) {
        //  const char* lo = "3#pi/4";
        //  const char* hi = "#pi";
        //  myText (0.25, 0.85, kBlack, Form ("%s < |#Delta#phi| < %s", lo, hi), 0.06);
        //}
      }
      for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
        if (iCent == iCentLo+1) {
          c->cd (iCent-iCentLo + 1);
          if (iPtZ == iPtZLo) {
            myText (0.656-0.08, 0.875, kBlack, "#it{ee}", 0.06);
            myText (0.736-0.08, 0.875, kBlack, "#it{#mu#mu}", 0.06);
            myText (0.816-0.08, 0.875, kBlack, "Comb.", 0.06);
            myMarkerTextNoLine (0.720-0.08, 0.810-0.10*(iPtZ-iPtZLo), colors[2*(iPtZ-iPtZLo)+1], kFullCircle, "", 1.3, 0.06);
            myMarkerTextNoLine (0.800-0.08, 0.810-0.10*(iPtZ-iPtZLo), colors[2*(iPtZ-iPtZLo)+2], kFullCircle, "", 1.3, 0.06);
            myMarkerTextNoLine (0.880-0.08, 0.810-0.10*(iPtZ-iPtZLo), colors[2*(iPtZ-iPtZLo)+3], kFullCircle, "", 1.3, 0.06);
          }
        }
        else if (iCent == iCentLo+2) {
          c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo)+iCent);
          if (iPtZ == nPtZBins-1) myText (0.510, 0.85, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.07);
          else                    myText (0.430, 0.85, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.07);
        }
      } // end loop over iPtZ
    } // end loop over iCent
  }

  for (short iCent = iCentLo; iCent < iCentHi; iCent++) {
    for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
      c->cd ((iPtZ-iPtZLo)*(iCentHi-iCentLo) + iCent-iCentLo + 1);
      const float xmin = useTrkPt ? trk_min_pt : xhZBins[nPtZBins-1][0];
      const float xmax = useTrkPt ? pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]] : xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]];
      TLine* l = new TLine (xmin, 1, xmax, 1);
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kPink-8);
      l->Draw ("same");
    }
  } // end loop over iCent

  c->SaveAs (Form ("%s/IAA/iaa_%s_dPtZ_SpcComp.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Writes cubic polynomial fits of IAA to a file.
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: WriteIAAs () {

  TDirectory* _gDirectory = gDirectory;

  const char* outFileName = "DataAnalysis/Nominal/data18hi_iaa_fits.root"; 
  TFile* outFile = new TFile (Form ("%s/%s", rootPath.Data (), outFileName), "recreate");

  TF1**** f_trk_pt_ptz_iaa = Get3DArray <TF1*> (3, nPtZBins, numCentBins);
  TF1**** f_trk_xhz_ptz_iaa = Get3DArray <TF1*> (3, nPtZBins, numCentBins);

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent] = new TF1 (Form ("f_trk_pt_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "[0]+[1]*log(x)+[2]*(log(x))^2+[3]*(log(x))^3", allPtchBins[0], allPtchBins[maxNPtchBins]);
        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (0, 1);
        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (1, 0);
        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (2, 0);
        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (3, 0);
        h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Fit (f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent], "RN0Q");

        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent] = new TF1 (Form ("f_trk_xhz_ptz_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "[0]+[1]*log(x)+[2]*(log(x))^2+[3]*(log(x))^3", allXhZBins[0], allXhZBins[maxNXhZBins]);
        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (0, 1);
        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (1, 0);
        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (2, 0);
        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->SetParameter (3, 0);
        h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Fit (f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent], "RN0Q");

        f_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Write ();
        h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Write ();
        f_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Write ();
        h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Write ();
      }
    }
  }

  outFile->Close ();

  Delete3DArray (&f_trk_pt_ptz_iaa, 3, nPtZBins, numCentBins);
  Delete3DArray (&f_trk_xhz_ptz_iaa, 3, nPtZBins, numCentBins);

  _gDirectory->cd ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots signal-to-background panels for data
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotSignalToBkg (const bool useTrkPt, const short iSpc) {
  CalculateSigToBkg ();

  TCanvas* c = new TCanvas (Form ("c_stb_%s", useTrkPt ? "pTch" : "xhZ"), "", 1200, 400);
  c->Divide (3, 1);

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    c->cd (iPtZ-2+1);
    gPad->SetLogx ();
    gPad->SetLogy ();

    TH1D* htemp = new TH1D (Form ("htemp_iPtZ%i", iPtZ), "", (useTrkPt ? nPtchBins : nXhZBins)[iPtZ], (useTrkPt ? pTchBins : xhZBins)[iPtZ]);

    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
    yax->SetTitle ("Y / Y_{bkg}");

    htemp->SetMarkerSize (0);
    htemp->SetLineWidth (0);

    xax->SetMoreLogLabels ();

    yax->SetRangeUser (8e-4, 2e4);

    htemp->SetLineWidth (0);
    htemp->SetMarkerSize (0);
  
    htemp->DrawCopy ("hist ][");
    SaferDelete (&htemp);

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
      double markerSize = (markerStyle == kOpenDiamond || markerStyle == kFullDiamond ? 2.0 : 1.2);

      TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_sig_to_bkg : h_trk_xhz_ptz_sig_to_bkg)[iSpc][iPtZ][iCent]);
      RecenterGraph (g);
      ResetXErrors (g);
      deltaize (g, 0.08*(iCent - 0.5*numCentBins), true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
      ResetXErrors (g);
      
      g->SetMarkerStyle (markerStyle);
      g->SetMarkerSize (markerSize);
      g->SetMarkerColor (finalColors[iCent]);
      g->SetLineColor (finalColors[iCent]);
      g->SetLineWidth (2);
      ((TGAE*) g->Clone ())->Draw ("P");

      if (IsFullMarker (markerStyle)) {
        markerStyle = FullToOpenMarker (markerStyle);

        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (markerSize);
        g->SetLineWidth (0);
        g->SetMarkerColor (kBlack);
  
        ((TGAE*) g->Clone ())->Draw ("P");
      }

      SaferDelete (&g);
    } // end loop over iCent

    if (iPtZ == 2) {
      myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.060);
      myText (0.22, 0.81, kBlack, "Pb+Pb, 5.02 TeV, 1.7 nb^{-1}", 0.050);
      myText (0.22, 0.74, kBlack, "#it{pp}, 5.02 TeV, 260 pb^{-1}", 0.050);
      myText (0.56, 0.24, kBlack, "15 < #it{p}_{T}^{Z} < 30 GeV", 0.050);
    }
    else if (iPtZ == 3) {
      myText (0.56, 0.24, kBlack, "30 < #it{p}_{T}^{Z} < 60 GeV", 0.050);
      myMarkerTextNoLine (0.30, 0.88, finalColors[0], kOpenCircle, "#it{pp}", 1.2, 0.05);
      myMarkerTextNoLine (0.30, 0.815, finalColors[1], markerStyles[0], "30-80%", 1.2, 0.05);
      myMarkerTextNoLine (0.30, 0.815, kBlack, FullToOpenMarker (markerStyles[0]), "", 1.2, 0.05);
    }
    else if (iPtZ == 4) {
      myMarkerTextNoLine (0.30, 0.88, finalColors[2], markerStyles[1], "10-30%", 1.2, 0.05);
      myMarkerTextNoLine (0.30, 0.88, kBlack, FullToOpenMarker (markerStyles[1]), "", 1.2, 0.05);
      myMarkerTextNoLine (0.30, 0.815, finalColors[3], markerStyles[2], "0-10%", 2.0, 0.05);
      myMarkerTextNoLine (0.30, 0.815, kBlack, FullToOpenMarker (markerStyles[2]), "", 2.0, 0.05);
      myText (0.62, 0.24, kBlack, "#it{p}_{T}^{Z} > 60 GeV", 0.05);
    }
  } // end loop over iPtZ

  c->SaveAs (Form ("%s/TrkYields/sigToBkg_%s_dPtZ.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots mean track distributions vs. pT^Z
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotTrackMeans (const bool useTrkPt, const short iSpc) {
  const char* canvasName = Form ("c_mean_%s_%s", useTrkPt ? "pTch" : "xhZ", iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 800);
    gDirectory->Add (c);
  }
  c->cd ();

  if (!canvasExists) {
    TH1D* htemp = new TH1D ("htemp", "", 1, 5, 120);

    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();
  
    xax->SetRangeUser (5, 120);
    yax->SetRangeUser (0, useTrkPt ? 15 : 0.25);

    xax->SetMoreLogLabels ();

    xax->SetTitle ("#LT#it{p}_{T}^{Z}#GT [GeV]");
    yax->SetTitle (useTrkPt ? "#LT#it{p}_{T}^{ch}#GT [GeV]" : "#LT#it{x}_{hZ}#GT");

    htemp->SetLineWidth (0);
    htemp->SetMarkerSize (0);
  
    htemp->DrawCopy ("hist ][");
    SaferDelete (&htemp);

    if (useTrkPt) {
      myText (0.22, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.055);
      if (iSpc == 0)      myText (0.72, 0.87, kBlack, "#it{Z} #rightarrow #it{ee}", 0.045);
      else if (iSpc == 1) myText (0.72, 0.87, kBlack, "#it{Z} #rightarrow #mu#mu", 0.045);
      myOnlyBoxText (0.28, 0.815, 1.4, finalFillColors[0], finalColors[0], 1, "", 0.045, 1001, 0.3);
      myOnlyBoxText (0.28, 0.755, 1.4, finalFillColors[1], finalColors[1], 1, "", 0.045, 1001, 0.3);
      myOnlyBoxText (0.50, 0.815, 1.4, finalFillColors[2], finalColors[2], 1, "", 0.045, 1001, 0.3);
      myOnlyBoxText (0.50, 0.755, 1.4, finalFillColors[3], finalColors[3], 1, "", 0.045, 1001, 0.3);
      myMarkerText (0.28, 0.815, finalColors[0], markerStyles[0], "#it{pp}", 1.5, 0.045);
      myMarkerText (0.28, 0.755, finalColors[1], markerStyles[1], "30-80%", 1.5, 0.045);
      myMarkerText (0.50, 0.815, finalColors[2], markerStyles[2], "10-30%", 2.4, 0.045);
      myMarkerText (0.50, 0.755, finalColors[3], markerStyles[3], "0-10%", 1.5, 0.045);
      myText (0.46, 0.27, kBlack, "Pb+Pb, 5.02 TeV, 1.7 nb^{-1}", 0.045);
      myText (0.46, 0.21, kBlack, "#it{pp}, 5.02 TeV, 260 pb^{-1}", 0.045);
      myText (0.22, 0.68, kBlack, "2 < #it{p}_{T}^{ch} < 240 GeV", 0.045);
    }
    else {
      myText (0.22, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.055);
      if (iSpc == 0)      myText (0.72, 0.87, kBlack, "#it{Z} #rightarrow #it{ee}", 0.045);
      else if (iSpc == 1) myText (0.72, 0.87, kBlack, "#it{Z} #rightarrow #mu#mu", 0.045);
      myOnlyBoxText (0.28, 0.815, 1.4, finalFillColors[0], finalColors[0], 1, "", 0.045, 1001, 0.3);
      myOnlyBoxText (0.28, 0.755, 1.4, finalFillColors[1], finalColors[1], 1, "", 0.045, 1001, 0.3);
      myOnlyBoxText (0.50, 0.815, 1.4, finalFillColors[2], finalColors[2], 1, "", 0.045, 1001, 0.3);
      myOnlyBoxText (0.50, 0.755, 1.4, finalFillColors[3], finalColors[3], 1, "", 0.045, 1001, 0.3);
      myMarkerText (0.28, 0.815, finalColors[0], markerStyles[0], "#it{pp}", 1.5, 0.045);
      myMarkerText (0.28, 0.755, finalColors[1], markerStyles[1], "30-80%", 1.5, 0.045);
      myMarkerText (0.50, 0.815, finalColors[2], markerStyles[2], "10-30%", 2.4, 0.045);
      myMarkerText (0.50, 0.755, finalColors[3], markerStyles[3], "0-10%", 1.5, 0.045);
      myText (0.22, 0.27, kBlack, "Pb+Pb, 5.02 TeV, 1.7 nb^{-1}", 0.045);
      myText (0.22, 0.21, kBlack, "#it{pp}, 5.02 TeV, 260 pb^{-1}", 0.045);
      myText (0.22, 0.68, kBlack, "1/15 < #it{x}_{hZ} < 2", 0.045);
    }
  }

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    TGAE* g = (useTrkPt ? g_trk_avg_pt_ptz : g_trk_avg_xhz_ptz)[iSpc][iCent];

    if (!plotAsSyst) {
      g->SetMarkerStyle (markerStyles[iCent]);
      g->SetMarkerSize ((markerStyles[iCent] == kOpenDiamond || markerStyles[iCent] == kFullDiamond) ? 2.4 : 1.5);
      g->SetMarkerColor (finalColors[iCent]);
      g->SetLineColor (finalColors[iCent]);
      g->SetLineWidth (2);
      g->Draw ("P");
    } else {
      g->SetMarkerSize (0); 
      //g->SetLineWidth (0);
      g->SetLineWidth (1);
      g->SetLineColor (finalColors[iCent]);
      g->SetFillColorAlpha (finalFillColors[iCent], 0.3);
      ((TGAE*)g->Clone ())->Draw ("5P");
      g->Draw ("2P");
    }

  } // end loop over iCent

  c->SaveAs (Form ("%s/TrkYields/Avg_%s_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb")));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot UE subtracted Y(pT or xZh) binned in Z Pt with power law fits overlaid
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PlotSubYields_dPtZ_Fits (short pSpc) {
  //const double padRatio = 0.9; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  //const double dPadY = padRatio / (padRatio+1.0);

  if (pSpc != 0 && pSpc != 1)
    pSpc = 2;
 
  const char* spc = (pSpc == 0 ? "ee" : (pSpc == 1 ? "mumu" : "comb")); 

  const char* canvasName = Form ("c_subYields_dPtZ_%s", spc);
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1500, 1000);
    gDirectory->Add (c);
  }
  c->cd ();

  const double llMargin = 0.17;
  const double lrMargin = 0.024;
  const double clMargin = 0.032;
  const double crMargin = 0.024;
  const double rlMargin = 0.032;
  const double rrMargin = 0.040;

  const double deltaL = (1. - llMargin - lrMargin);
  const double deltaC = (1. - clMargin - crMargin);
  const double deltaR = (1. - rlMargin - rrMargin);

  const double a = (double) (deltaR * deltaC / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));
  const double b = (double) (deltaR * deltaL / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));

  const double xPadLCMiddle = a;
  const double xPadCRMiddle = a+b;

  double yPadMiddle = 0.5;

  TPad* luPad = new TPad (Form ("p_subYield_dPtZ_%s_luPad", spc), "", 0, yPadMiddle, xPadLCMiddle, 1);
  TPad* ldPad = new TPad (Form ("p_subYield_dPtZ_%s_ldPad", spc), "", 0, 0, xPadLCMiddle, yPadMiddle);
  TPad* cuPad = new TPad (Form ("p_subYield_dPtZ_%s_cuPad", spc), "", xPadLCMiddle, yPadMiddle, xPadCRMiddle, 1);
  TPad* cdPad = new TPad (Form ("p_subYield_dPtZ_%s_cdPad", spc), "", xPadLCMiddle, 0, xPadCRMiddle, yPadMiddle);
  TPad* ruPad = new TPad (Form ("p_subYield_dPtZ_%s_ruPad", spc), "", xPadCRMiddle, yPadMiddle, 1, 1);
  TPad* rdPad = new TPad (Form ("p_subYield_dPtZ_%s_rdPad", spc), "", xPadCRMiddle, 0, 1, yPadMiddle);

  luPad->SetLeftMargin (llMargin);
  luPad->SetRightMargin (lrMargin);
  cuPad->SetLeftMargin (clMargin);
  cuPad->SetRightMargin (crMargin);
  ruPad->SetLeftMargin (rlMargin);
  ruPad->SetRightMargin (rrMargin);
  ldPad->SetLeftMargin (llMargin);
  ldPad->SetRightMargin (lrMargin);
  cdPad->SetLeftMargin (clMargin);
  cdPad->SetRightMargin (crMargin);
  rdPad->SetLeftMargin (rlMargin);
  rdPad->SetRightMargin (rrMargin);

  luPad->Draw ();
  cuPad->Draw ();
  ruPad->Draw ();
  ldPad->Draw ();
  cdPad->Draw ();
  rdPad->Draw ();

  TPad* pads[6] = {luPad, cuPad, ruPad, ldPad, cdPad, rdPad};

  for (bool useTrkPt : {true, false}) {
    for (short iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
      pads[(useTrkPt ? 0 : 3) + iPtZ-2]->cd ();
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* htemp = new TH1D ("h", "", useTrkPt ? nPtchBins[iPtZ] : nXhZBins[iPtZ], useTrkPt ? pTchBins[iPtZ] : xhZBins[iPtZ]);

      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      useTrkPt ? htemp->GetXaxis ()->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : htemp->GetXaxis ()->SetRangeUser (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
      useTrkPt ? htemp->GetYaxis ()->SetRangeUser (2e-3, 80) : yax->SetRangeUser (5e-3, 2e3);

      xax->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
      yax->SetTitle (useTrkPt ? "(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]" : "(1/N_{Z}) (d^{2}N_{ch} / d#it{x} d#Delta#phi)");

      xax->SetTitleFont (43);
      yax->SetTitleFont (43);
      //xax->SetLabelFont (43);
      yax->SetLabelFont (43);

      xax->SetTitleSize (30);
      yax->SetTitleSize (30);
      xax->SetTitleOffset (2.35);
      yax->SetTitleOffset (2.5);
      xax->SetLabelSize (0);
      yax->SetLabelSize (gPad != luPad && gPad != ldPad ? 0 : 27);

      htemp->SetLineWidth (0);
      htemp->SetMarkerSize (0);
  
      htemp->DrawCopy ("hist ][");
      SaferDelete (&htemp);

      if (useTrkPt) {
        TLatex* tl = new TLatex ();
        tl->SetTextFont (43);
        tl->SetTextSize (25);
        tl->SetTextAlign (21);
        double yoff = 9.33e-4;
        tl->DrawLatex (1,  yoff, "1");
        tl->DrawLatex (2,  yoff, "2");
        tl->DrawLatex (3,  yoff, "3");
        tl->DrawLatex (4,  yoff, "4");
        tl->DrawLatex (5,  yoff, "5");
        tl->DrawLatex (6,  yoff, "6");
        tl->DrawLatex (7,  yoff, "7");
        tl->DrawLatex (10, yoff, "10");
        tl->DrawLatex (20, yoff, "20");
        tl->DrawLatex (30, yoff, "30");
        tl->DrawLatex (40, yoff, "40");
        tl->DrawLatex (60, yoff, "60");
      }
      else {
        TLatex* tl = new TLatex ();
        tl->SetTextFont (43);
        tl->SetTextSize (25);
        tl->SetTextAlign (21);
        double yoff = 1.95e-3;
        if (iPtZ > 2) {
          if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
          else tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
        }
        tl->DrawLatex (1e-1,  yoff, "10^{-1}");
        tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
        tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
        tl->DrawLatex (1,     yoff, "1");
      }
    } // end loop over iPtZ
  } // end loop over useTrkPt

  for (bool useTrkPt : {true, false}) {
    for (short iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
      pads[(useTrkPt ? 0 : 3) + iPtZ-2]->cd ();
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        TGAE* g = GetTGAE ((useTrkPt ? h_trk_pt_ptz_sub : h_trk_xhz_ptz_sub)[pSpc][iPtZ][iCent]);
        RecenterGraph (g);
        ResetXErrors (g);
        deltaize (g, 0.08*(iCent - 0.5*numCentBins), true, useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
        ResetXErrors (g);

        if (plotAsSyst) {
          SetConstantXErrors (g, 0.032, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);

          g->SetMarkerSize (0);
          g->SetLineWidth (1);
          g->SetLineColor (finalColors[iCent]);
          g->SetFillColorAlpha (finalFillColors[iCent], 0.3);

          ((TGAE*) g->Clone ())->Draw ("5P");

          SaferDelete (&g);
        }
        else {
          Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
          double markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (markerSize);
          g->SetLineWidth (2);
          g->SetMarkerColor (finalColors[iCent]);
          g->SetLineColor (finalColors[iCent]);

          ((TGAE*) g->Clone ())->Draw ("P");

          if (iCent > 0) {
            markerStyle = FullToOpenMarker (markerStyle);
            markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  
            g->SetMarkerStyle (markerStyle);
            g->SetMarkerSize (markerSize);
            g->SetLineWidth (0);
            g->SetMarkerColor (kBlack);
  
            ((TGAE*) g->Clone ())->Draw ("P");
          }

          SaferDelete (&g);
        }
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over useTrkPt

  if (!plotAsSyst) {
    for (bool useTrkPt : {true, false}) {
      for (short iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
        pads[(useTrkPt ? 0 : 3) + iPtZ-2]->cd ();
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          TF1* f = (TF1*) (useTrkPt ? f_trk_fit_pt_ptz : f_trk_fit_xhz_ptz)[pSpc][iPtZ][iCent]->Clone ();
          if (useTrkPt) f->SetParameter (0, f->GetParameter (0) * pow (pTchBins[iPtZ][0]/pTchBins[iPtZ][nPtchBins[iPtZ]], 0.5 * (0.08*(iCent - 0.5*numCentBins)) * f->GetParameter (1)));
          else          f->SetParameter (0, f->GetParameter (0) * pow (xhZBins[iPtZ][0]/xhZBins[iPtZ][nXhZBins[iPtZ]],    0.5 * (0.08*(iCent - 0.5*numCentBins)) * f->GetParameter (1)));
          f->SetLineColor (finalColors[iCent]);
          f->Draw ("same");
        } // end loop over iCent
      } // end loop over iPtZ
    } // end loop over useTrkPt
  }

  luPad->cd ();
  myText (0.21, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.06); myText (0.32, 0.79, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}", 0.050);
  myText (0.32, 0.72, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02, 1.4-1.7 nb^{-1}", 0.050);

  ldPad->cd ();
  myText (0.21, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.06);
  myText (0.32, 0.79, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}", 0.050);
  myText (0.32, 0.72, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02, 1.4-1.7 nb^{-1}", 0.050);

  luPad->cd ();
  myText (0.60, 0.86, kBlack, "#it{p}_{T}^{Z} = 15-30 GeV", 0.055);
  myText (0.60, 0.65, kBlack, Form ("Z #rightarrow #it{%s}", pSpc == 0 ? "ee" : (pSpc == 1 ? "#mu#mu" : "ll")), 0.050);
  cuPad->cd ();
  myText (0.56, 0.86, kBlack, "#it{p}_{T}^{Z} = 30-60 GeV", 0.055);
  ruPad->cd ();
  myText (0.61, 0.86, kBlack, "#it{p}_{T}^{Z} = 60+ GeV", 0.055);

  ldPad->cd ();
  myText (0.60, 0.86, kBlack, "#it{p}_{T}^{Z} = 15-30 GeV", 0.055);
  myText (0.60, 0.65, kBlack, Form ("Z #rightarrow #it{%s}", pSpc == 0 ? "ee" : (pSpc == 1 ? "#mu#mu" : "ll")), 0.050);
  cdPad->cd ();
  myText (0.56, 0.86, kBlack, "#it{p}_{T}^{Z} = 30-60 GeV", 0.055);
  rdPad->cd ();
  myText (0.61, 0.86, kBlack, "#it{p}_{T}^{Z} = 60+ GeV", 0.055);

  cuPad->cd ();
  myMarkerAndBoxAndLineText (0.22, 0.280, 2.0, 1001, finalFillColors[0], 0.40, finalColors[0], kOpenCircle,     1.8, "#it{pp}", 0.052);
  myMarkerAndBoxAndLineText (0.22, 0.220, 2.0, 1001, finalFillColors[1], 0.40, finalColors[1], markerStyles[0], 1.8, "Pb+Pb 30-80\%", 0.052);
  ruPad->cd ();
  myMarkerAndBoxAndLineText (0.22, 0.280, 2.0, 1001, finalFillColors[2], 0.40, finalColors[2], markerStyles[1], 1.8, "Pb+Pb 10-30\%", 0.052);
  myMarkerAndBoxAndLineText (0.22, 0.220, 2.0, 1001, finalFillColors[3], 0.40, finalColors[3], markerStyles[2], 2.4, "Pb+Pb 0-10\%", 0.052);

  cdPad->cd ();
  myMarkerAndBoxAndLineText (0.22, 0.280, 2.0, 1001, finalFillColors[0], 0.40, finalColors[0], kOpenCircle,     1.8, "#it{pp}", 0.052);
  myMarkerAndBoxAndLineText (0.22, 0.220, 2.0, 1001, finalFillColors[1], 0.40, finalColors[1], markerStyles[0], 1.8, "Pb+Pb 30-80\%", 0.052);
  rdPad->cd ();
  myMarkerAndBoxAndLineText (0.22, 0.280, 2.0, 1001, finalFillColors[2], 0.40, finalColors[2], markerStyles[1], 1.8, "Pb+Pb 10-30\%", 0.052);
  myMarkerAndBoxAndLineText (0.22, 0.220, 2.0, 1001, finalFillColors[3], 0.40, finalColors[3], markerStyles[2], 2.4, "Pb+Pb 0-10\%", 0.052);
    
  c->SaveAs (Form ("%s/TrkYields/subYields_withFits_%s.pdf", plotPath.Data (), spc));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Prints IAA values for the bin specified
////////////////////////////////////////////////////////////////////////////////////////////////
void PhysicsAnalysis :: PrintIAA (const bool printErrs, const bool useTrkPt, const short iCent, const short iPtZ, const short iSpc) {
  TH1D* h = nullptr;
  if (useTrkPt)
    h = h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent];
  else
    h = h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent];
  cout << name << "\t";
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    cout << h->GetBinContent (ix) << "\t";
    if (printErrs)
      cout << h->GetBinError (ix) << "\t";
    cout << endl;
  }
}

#endif

