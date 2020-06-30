#ifndef __TreeVariables_h__
#define __TreeVariables_h__

#include <TTree.h>
#include <vector>

/** 
 * This file defines the TreeVariables struct used by this analysis.
 * Author: Jeff Ouellette
 * Dated: 8/20/2018
 */

using namespace std;

namespace ZHadronAnalysis {

struct TreeVariables {
  private:
    TTree* tree = nullptr;

    int currEntry = -1;
    bool isMC = false;
    bool getCollisionRateInfo = false;
    bool getVertices = false;
    bool getFCals = false;
    bool getZdc = false;
    bool getTracks = false;
    bool getElectrons = false;
    bool getClusters = false;
    bool getMuons = false;
    bool getJets = false;
    bool getTruthTracks = false;
    bool getTruthElectrons = false;
    bool getTruthMuons = false;
    bool getTruthJets = false;

  public:

    // public functions
    TreeVariables (TTree* t, const bool _isMC = false);
    ~TreeVariables () {}
    void SetBranchAddresses ();
    void PrintAll (const long long entry);
    void GetEntry (const long long entry) { tree->GetEntry (entry); currEntry = entry; }
    int GetCurrEntry () { return currEntry; }

    // setter functions
    void SetGetCollisionRateInfo  (const bool _getCollisionRateInfo = true);
    void SetGetVertices           (const bool _getVertices = true);
    void SetGetFCals              (const bool _getFCals = true);
    void SetGetZdc                (const bool _getZdc = true);
    void SetGetTracks             (const bool _getTracks = true);
    void SetGetTruthTracks        (const bool _getTruthTracks = true);
    void SetGetElectrons          (const bool _getElectrons = true);
    void SetGetClusters           (const bool _getClusters = true);
    void SetGetTruthElectrons     (const bool _getTruthElectrons = true);
    void SetGetMuons              (const bool _getMuons = true);
    void SetGetTruthMuons         (const bool _getTruthMuons = true);
    void SetGetJets               (const bool _getJets = true);
    void SetGetTruthJets          (const bool _getTruthJets = true);

    // public variables
    unsigned int run_number = 0;
    unsigned int event_number = 0;
    unsigned int lumi_block = 0;
    bool passes_toroid = false;
    bool isOOTPU = false;
    bool BlayerDesyn = false;

    vector<float>* mcEventWeights = nullptr;
    // truth event info
    int nTruthEvt = 0;
    int nPart1[5];
    int nPart2[5];
    float impactParameter[5];
    int nColl[5];
    int nSpectatorNeutrons[5];
    int nSpectatorProtons[5];
    float eccentricity[5];
    float eventPlaneAngle[5];

    float actualInteractionsPerCrossing = 0;
    float averageInteractionsPerCrossing = 0;

    int nvert = 0;
    float   vert_x[30];
    float   vert_y[30];
    float   vert_z[30];
    int     vert_ntrk[30];
    int     vert_type[30];

    // fcal energies
    float fcalA_et = 0;
    float fcalC_et = 0;
    float fcalA_et_Cos2 = 0;
    float fcalC_et_Cos2 = 0;
    float fcalA_et_Sin2 = 0;
    float fcalC_et_Sin2 = 0;
    float fcalA_et_Cos3 = 0;
    float fcalC_et_Cos3 = 0;
    float fcalA_et_Sin3 = 0;
    float fcalC_et_Sin3 = 0;
    float fcalA_et_Cos4 = 0;
    float fcalC_et_Cos4 = 0;
    float fcalA_et_Sin4 = 0;
    float fcalC_et_Sin4 = 0;

    // ZDC
    float  ZdcCalibEnergy_A   = 0;
    float  ZdcCalibEnergy_C   = 0;
    bool   L1_ZDC_A           = false;
    bool   L1_ZDC_A_tbp       = false;
    bool   L1_ZDC_A_tap       = false;
    bool   L1_ZDC_A_tav       = false;
    float  L1_ZDC_A_prescale  = 0;
    bool   L1_ZDC_C           = false;
    bool   L1_ZDC_C_tbp       = false;
    bool   L1_ZDC_C_tap       = false;
    bool   L1_ZDC_C_tav       = false;
    float  L1_ZDC_C_prescale  = 0;
 
    // tracking info (0th or primary vertex only)
    int     ntrk = 0;
    float   trk_pt[10000];
    float   trk_eta[10000];
    float   trk_phi[10000];
    float   trk_charge[10000];
    bool    trk_HItight[10000];
    bool    trk_HIloose[10000];
    bool    trk_TightPrimary[10000];
    float   trk_d0[10000];
    float   trk_d0sig[10000];
    float   trk_z0[10000];
    float   trk_z0sig[10000];
    float   trk_theta[10000];
    float   trk_vz[10000];
    int     trk_nBLayerHits[10000];
    int     trk_nBLayerSharedHits[10000];
    int     trk_nPixelHits[10000];
    int     trk_nPixelDeadSensors[10000];
    int     trk_nPixelSharedHits[10000];
    int     trk_nSCTHits[10000];
    int     trk_nSCTDeadSensors[10000];
    int     trk_nSCTSharedHits[10000];
    int     trk_nTRTHits[10000];
    int     trk_nTRTSharedHits[10000];
    float   trk_prob_truth[10000];
    float   trk_truth_pt[10000];
    float   trk_truth_eta[10000];
    float   trk_truth_phi[10000];
    float   trk_truth_charge[10000];
    int     trk_truth_type[10000];
    int     trk_truth_orig[10000];
    int     trk_truth_barcode[10000];
    int     trk_truth_pdgid[10000];
    float   trk_truth_vz[10000];
    int     trk_truth_nIn[10000];
    bool    trk_truth_isHadron[10000];

    int     electron_n = 0;
    float   electron_pt_precalib[40];
    float   electron_pt[40];
    float   electron_eta[40];
    float   electron_phi[40];
    int     electron_charge[40];
    bool    electron_lhtight[40];
    bool    electron_lhmedium[40];
    bool    electron_lhloose[40];
    bool    electron_lhmedium_hi[40];
    bool    electron_lhloose_hi[40];
    bool    electron_matched[40];
    float   electron_etcone20[40];
    float   electron_etcone30[40];
    float   electron_etcone40[40];
    float   electron_topoetcone20[40];
    float   electron_topoetcone30[40];
    float   electron_topoetcone40[40];
    float   electron_id_track_pt[40];
    float   electron_id_track_eta[40];
    float   electron_id_track_phi[40];
    float   electron_id_track_charge[40];
    float   electron_id_track_d0[40];
    float   electron_id_track_d0sig[40];
    float   electron_id_track_z0[40];
    float   electron_id_track_z0sig[40];
    float   electron_id_track_theta[40];
    float   electron_id_track_vz[40];
    bool    electron_id_track_tightprimary[40];
    bool    electron_id_track_hiloose[40];
    bool    electron_id_track_hitight[40];
    float   electron_pt_sys[40];
    float   electron_eta_sys[40];
    float   electron_phi_sys[40];

    // egamma calo cluster info, see https://ucatlas.github.io/RootCoreDocumentation/2.4.28/dc/d4b/CaloCluster__v1_8h_source.html
    int             Cluster_n         = 0;
    vector<float>*  Cluster_pt        = nullptr;  
    vector<float>*  Cluster_et        = nullptr;
    vector<float>*  Cluster_eta       = nullptr;
    vector<float>*  Cluster_phi       = nullptr;
    vector<float>*  Cluster_energyBE  = nullptr;
    vector<float>*  Cluster_etaBE     = nullptr;
    vector<float>*  Cluster_phiBE     = nullptr;
    vector<float>*  Cluster_calE      = nullptr;
    vector<float>*  Cluster_calEta    = nullptr;
    vector<float>*  Cluster_calPhi    = nullptr;
    vector<int>*    Cluster_size      = nullptr;
    vector<int>*    Cluster_status    = nullptr;

    int     muon_n = 0;
    float   muon_pt_precalib[40];
    float   muon_pt[40];
    float   muon_ms_pt_precalib[40];
    float   muon_ms_pt[40];
    float   muon_eta[40];
    float   muon_phi[40];
    int     muon_charge[40];
    bool    muon_tight[40];
    bool    muon_medium[40];
    bool    muon_loose[40];
    bool    muon_matched[40];
    float   muon_etcone20[40];
    float   muon_etcone30[40];
    float   muon_etcone40[40];
    float   muon_topoetcone20[40];
    float   muon_topoetcone30[40];
    float   muon_topoetcone40[40];
    float   muon_id_track_pt[40]; //!
    float   muon_id_track_eta[40]; //!
    float   muon_id_track_phi[40]; //!
    float   muon_id_track_charge[40]; //!
    float   muon_id_track_d0[40]; //!
    float   muon_id_track_d0sig[40]; //!
    float   muon_id_track_z0[40]; //!
    float   muon_id_track_z0sig[40]; //!
    float   muon_id_track_theta[40]; //!
    float   muon_id_track_vz[40]; //!
    bool    muon_id_track_tightprimary[40]; //!
    bool    muon_id_track_hiloose[40]; //!
    bool    muon_id_track_hitight[40]; //!
    float   muon_ms_track_pt[40]; //!
    float   muon_ms_track_eta[40]; //!
    float   muon_ms_track_phi[40]; //!
    float   muon_ms_track_charge[40]; //!
    float   muon_ms_track_d0[40]; //!
    float   muon_ms_track_d0sig[40]; //!
    float   muon_ms_track_z0[40]; //!
    float   muon_ms_track_z0sig[40]; //!
    float   muon_ms_track_theta[40]; //!
    float   muon_ms_track_vz[40]; //!
  
    vector <vector <double>>* muon_pt_sys;
    vector <vector <double>>* muon_eta_sys;
    vector <vector <double>>* muon_phi_sys;

    int     akt4emtopo_jet_n = 0;
    float   akt4emtopo_jet_pt[40];
    float   akt4emtopo_jet_eta[40];
    float   akt4emtopo_jet_phi[40];
    float   akt4emtopo_jet_e[40];

    int     truth_electron_n = 0;
    float   truth_electron_pt[1000];
    float   truth_electron_eta[1000];
    float   truth_electron_phi[1000];
    int     truth_electron_charge[1000];
    int     truth_electron_barcode[1000];

    int     truth_muon_n = 0;
    float   truth_muon_pt[1000];
    float   truth_muon_eta[1000];
    float   truth_muon_phi[1000];
    int     truth_muon_charge[1000];
    int     truth_muon_barcode[1000];

    int     truth_trk_n = 0;
    float   truth_trk_pt[10000];
    float   truth_trk_eta[10000];
    float   truth_trk_phi[10000];
    float   truth_trk_charge[10000];
    int     truth_trk_pdgid[10000];
    int     truth_trk_barcode[10000];
    bool    truth_trk_isHadron[10000];

    int     truth_jet_n = 0;
    float   truth_jet_pt[40];
    float   truth_jet_eta[40];
    float   truth_jet_phi[40];
    float   truth_jet_e[40];
   
};

} // end namespace

#endif
