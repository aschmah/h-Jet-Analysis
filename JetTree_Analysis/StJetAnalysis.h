#ifndef StJetAnalysis_hh
#define StJetAnalysis_hh

#define USEEVE

//#include "StarClassLibrary/SystemOfUnits.h" // ALEX
#include "TString.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCutG.h"
#include "TCanvas.h"
#include "TTree.h"
//#include "StMessMgr.h" // ALEX
//#include <iostream.h>
#include "TChain.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TFile.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TGaxis.h"
#include <fstream>
#include "TMath.h"
#include "TColor.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH3D.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLCamera.h"
#include "TGLPerspectiveCamera.h"
#include "TGFrame.h"
#include "TGLUtil.h"
#include "TGLLightSet.h"
#include "TGLCameraOverlay.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"

//#include "StJetAnalysisLinkDef.h"



using namespace std;

//#include "/home/ceres/schmah/ALICE/Jet_Analysis/StJetTrackEventLinkDef.h"
//#include "/home/ceres/schmah/ALICE/Jet_Analysis/StJetTrackEvent.h"
#include "../StJetTrackEvent.h"
#include "../StJetTrackEventLinkDef.h"

ClassImp(StEMCal)
ClassImp(StJetTrackParticle)
ClassImp(StJetTrackEvent)

// #include "StarClassLibrary/StThreeVectorF.hh" // ALEX
//#include "StRoot/StRefMultCorr/StRefMultCorr.h" // ALEX
// #include "StRoot/StRefMultCorr/CentralityMaker.h" // ALEX

#include "fastjet/config.h"             // will allow a test for FJ3
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"


//#include "StPhysicalHelixD.hh" // ALEX

using namespace fastjet;
using namespace std;


#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#if defined(USEEVE)
#include "TEveBox.h"
#include "TEveArrow.h"
#include <TEveManager.h>
#include "TEveLine.h"
#include "TEvePointSet.h"
#endif

typedef std::vector< std::vector<Int_t> >  Two_dim_Int_vector;
typedef std::vector< std::vector< std::vector<Int_t> > >  Three_dim_Int_vector;
typedef std::vector< std::vector< std::vector<StJetTrackEvent> > >  Three_dim_StJetTrackEvent_vector;

static const Double_t Pi = TMath::Pi();


#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"



class StJetAnalysis {
public:
    StJetAnalysis();
    ~StJetAnalysis();



    vector<Double_t> Get_Helix_params_from_kine(TLorentzVector TLV_particle, TVector3 TV3_vertex, Double_t charge);
    void  setGeomDir(TString iGeomDir) {GeomDir = iGeomDir;};
    void  setInListDir(TString iInListDir) {InListDir = iInListDir;};
    void  setSEList(TString iSEList) {SEList = iSEList;};
    void  setInputDir(TString ipinputdir) {pinputdir = ipinputdir;};
    void  setInputDirPYTHIA(TString ipinputdirPYTHIA) {pinputdirPYTHIA = ipinputdirPYTHIA;};
    void  setStartEvent(Long64_t istart_event) {start_event_use = istart_event;};
    void  setStopEvent(Long64_t istop_event) {stop_event_use = istop_event;};
    void  setOutputfile(TString ioutputfile_name) {outputfile_name = ioutputfile_name;};
    void  setEnergy(Int_t iBeamTimeNum) {eBeamTimeNum = iBeamTimeNum;};
    void  setRandom(Int_t iRandom) {eRandom = iRandom;};
    void  setMode(Int_t iMode) {eMode = iMode;};
    void  setIn_Mode(Int_t iIn_Mode) {eIn_Mode = iIn_Mode;};
    void  setOutdir(TString iOutdir) {eOutdir = iOutdir;};
    void  setSampleHisto(TString iSampleHist) {eSampleHist = iSampleHist;};
    void  setSE_Et_Histo(TString iSE_Et_Hist) {eSE_Et_Hist = iSE_Et_Hist;};
    void  setME_Et_Histo(TString iME_Et_Hist) {eME_Et_Hist = iME_Et_Hist;};
    void  setSuffix(TString iSuffix) {eSuffix = iSuffix;};
    void  setReCenteringFile(TString ire_centering_name) {re_centering_name = ire_centering_name;};
    void  setPYTHIA_eff_factor(Double_t iPYTHIA_eff_factor) {ePYTHIA_eff_factor = iPYTHIA_eff_factor;};
    void  setBins(Int_t iz_bin, Int_t imult_bin, Int_t iPsi_bin) {ez_bin = iz_bin; emult_bin = imult_bin; ePsi_bin = iPsi_bin;};
    void  setCentrality(Int_t iCentrality) {eCentrality = iCentrality;};
    void  setPrim_Glob(Int_t iflab_prim_glob) {eflab_prim_glob = iflab_prim_glob;};
    void  setJet_R(Double_t iJet_R) {eJet_R = iJet_R;};
    void  setBkg_R(Double_t iBkg_R) {eBkg_R = iBkg_R;};
    void  setRem_n_hardest(Int_t iRem_n_hardest) {eRem_n_hardest = iRem_n_hardest;};
    void  set_Max_pt_down_scale(Double_t imax_pt_down_scale) {emax_pt_down_scale = imax_pt_down_scale;};
    void  set_N_vertex_mult_Psi_bins(Int_t N_z_vertex_bins_set, Int_t N_mult_bins_set, Int_t N_Psi_bins_set) {N_z_vertex_bins = N_z_vertex_bins_set; N_mult_bins = N_mult_bins_set; N_Psi_bins = N_Psi_bins_set;};
    void  SplitTrees();
    void  MakeJets(Int_t graphics, Int_t ME_flag);
    void  MakeME(Int_t i_ME_mult_bin, Int_t i_ME_EP_bin);
    Double_t MakeEventPlane();
    Int_t LoopEvent(Int_t LoopEvent, Int_t graphics);
    void  setMultHist(TString iMultHistName) {eMultHistName = iMultHistName;}

    //----------------------------------------------------------------------
    // Function for track reconstruction efficiency
    void      Read_ABC_params_eff_track_rec_function();             //read root file  for Global_ABC
    Double_t* Parameter_eff_track_rec_function(Int_t cent, Int_t runID,  Int_t PID, Int_t use_all);  //return A,B,C,runID for different centraltiy
    void      setEff_track_rec_parameter_file(TString iEff_file) {eEff_file = iEff_file;};
    Int_t     Get_runID_range_Index(Int_t runid);
    //----------------------------------------------------------------------


    void Init3DGraphics();
    void InitSplit();
    void InitJet();
    Long64_t ReadData();
    void WriteJet();
    Int_t Make();
    Int_t Finish();


private:

    TRandom ran;
    TRandom3 r3;
    TString HistName;
    char NoP[50];
    TRandom ran_gen;
    Int_t N_z_vertex_bins, N_mult_bins, N_Psi_bins;

    vector< vector< vector<StJetTrackEvent*> > > JetTrackEvent_Fill;
    vector< vector< vector<TTree*> > >           Tree_JetTrackEvent_Fill;
    vector< vector< vector<TFile*> > >           vec_Outputfiles;
    TFile* outputfile_jet;

    TVector3 Qvec_eta_pos, Qvec_eta_neg;
    vector< vector<TVector3> > vec_TV3_Qvec_eta;
    Double_t Psi_full_global;

    vector<PseudoJet> vec_PJ_particles;
    Int_t eflab_prim_glob, eRem_n_hardest;
    Double_t eJet_R, eBkg_R, emax_pt_down_scale;
    TString eSuffix, InListDir, GeomDir, eMultHistName;
    Double_t ePYTHIA_eff_factor;
    Int_t ez_bin, emult_bin, ePsi_bin, eCentrality;
    TString eOutdir, eSampleHist, eSE_Et_Hist, eME_Et_Hist;
    Int_t eMode, eIn_Mode;
    Int_t eBeamTimeNum;
    Int_t eRandom;
    Long64_t file_entries_SE_ME[2];
    TString SEList, pinputdir, pinputdirPYTHIA, outputfile_name, re_centering_name;
    TChain* input_SE_ME[2];
    TChain* input_PYTHIA[11]; // PYTHIA hard bins, only for mode 312 -> embedding
    StJetTrackEvent     *JetTrackEvent_PYTHIA[11]; // PYTHIA hard bins, only for mode 312 -> embedding
    StJetTrackParticle  *JetTrackParticle_PYTHIA[11]; // PYTHIA hard bins, only for mode 312 -> embedding
    Long64_t start_event_use, stop_event_use;
    TFile* Outputfile;
    TFile* Inputfile;
    TFile* Inputfile_Et[2];
    StJetTrackEvent     *JetTrackEvent;
    StJetTrackParticle  *JetTrackParticle;
    StEMCal             *JetEMCalParticle;
    StJetTrackParticle  *JetTrackParticle_Fill;

    vector<Int_t> vec_ME_mult_bins; // multiplicity
    vector<Int_t> vec_ME_EP_bins; // event plane
    vector< vector< vector<Int_t> > > vec_ME_event_index_mult;

    TH2D* h2D_dEdx_vs_p;
    TH1D* h_jet_sub_SE;
    TH1D* h_jet_sub_ME;
    TH1D* h_jet_SE;
    TH1D* h_jet_ME;
    TH1D* h_rho_sub_SE;
    TH1D* h_rho_sub_ME;
    TH1D* h_area_sub_SE;
    TH1D* h_area_sub_ME;
    TH1D* h_mult_particles;
    TH1D* h_input_mult_hist;
    vector<TH1D*> vec_h_input_mult_hist;
    TH2D* h2D_Psi_pos_vs_Psi_neg;
    TH2D* h2D_jet_SE_vs_dPsi_phi;
    TH1D* h_Psi_full;

    //----------------------------------------------------------------------
    // Parameters for track reconstruction efficiency
    // Centrality: 0-5,5-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80
    // runID: Index 0 is for all run ranges, then starting from 1: 12138025,12145021,12152017,12154022,12165032,12171017
    // PID: p, anti-p", pi+, pi-, K+, K-
    Double_t Global_ABC_Parameters_A[9][7][6];  //centrality(9),runID(7),PID(6)
    Double_t Global_ABC_Parameters_B[9][7][6];  //centrality(9),runID(7),PID(6)
    Double_t Global_ABC_Parameters_C[9][7][6];  //centrality(9),runID(7),PID(6)
    TString eEff_file;
    //----------------------------------------------------------------------


    //----------------------------------------------------------------------
    vector< vector<TH1D*> > vec_TH1D_TRD_geometry; // store for all 540 chambers the 8 corner vertices per detector
    TEveLine* TEveLine_beam_axis = NULL;
    TEvePointSet* TEveP_primary_vertex;
    TEveLine* TPL3D_helix = NULL;
    vector<TEveLine*> vec_TPL3D_helix;
    vector<TEveLine*> vec_TPL3D_helix_inner;
    vector<TEveLine*> vec_TPL3D_helix_hull;
    vector<TEveBox*> vec_eve_TRD_detector_box;
    vector<TEveBox*> vec_eve_EMCal_cluster;
    vector<TEveBox*> vec_eve_Jets;
    vector<TEveBox*> vec_eve_tracks_jets;
    //----------------------------------------------------------------------


    ClassDef(StJetAnalysis,1)

};

#endif
