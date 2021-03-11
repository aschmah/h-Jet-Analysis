#include "StJetAnalysis.h"

ClassImp(StJetAnalysis);



//----------------------------------------------------------------------------------------
static TRandom ran;
static TRandom3 r3;
static TString HistName;
static char NoP[50];
static TRandom ran_gen;
//----------------------------------------------------------------------------------------

//static StRefMultCorr* refmultCorrUtil; // ALEX
static Int_t erefMult_bin;
static Int_t erefMult_bin16;


static TString ALEX_EVENT_TREE   = "JetTrackEvent";
static TString ALEX_EVENT_BRANCH = "Events";

static Double_t jet_R;
static Double_t jet_R_background;
static const Int_t N_Beamtime            = 8; // 0 == 7.7, 1 == 11.5, 2 == 39, 3 == 62.4, 4 == 19.6, 5 == 27, 6 == 200, 7 == 14.5
static const Int_t N_jet_areas           = 16; //
static const Int_t N_z_vertex_bins       = 20; // 8
static const Int_t N_mult_bins           = 8; // 4
static const Int_t N_Psi_bins            = 4; // 5
static const Int_t N_global_bin          = N_Psi_bins*N_mult_bins*N_z_vertex_bins;
static const Double_t max_pt_val_embed   = 50.0; // maximum pt value embeded for Delta pt calculation
static const Double_t jet_pt_cut_closure = 3.0;

static const  Double_t track_eff_scaling_factor = 81499.0/73525.0; // ratio of reconstructed tracks (1000 events) for Nhitsfit > 14 and Nhitsfit > 19

static Int_t Remove_N_hardest; // 2 for same event, 0 for mixed event (if high pT tracks were removed or scaled down)
static const Int_t N_Et_bins               = 5;
static const Int_t N_leading_pt_bins       = 20; // 15
static const Double_t leading_pt_split_cut = 4.0; // cut at which the tracks are split into tracks with smaller pT // 2.5 // 1.5
static const Double_t leading_pt_split_val = 0.5;
static const Double_t max_pt_threshold     = 30.0; // 30
static Double_t max_pt_downscale_threshold; // Used only for eMode == 32, all particles with this momentum will be downscaled by downscale_factor
static const Double_t downscale_factor     = 0.1;
static const Int_t N_track_pt_bins_eta_phi = 3;
static const Double_t array_pt_bins_eta_phi[N_track_pt_bins_eta_phi] = {0.5,2.0,max_pt_threshold}; // binning in pT (upper values) for eta vs. phi 2D histograms
static const Double_t jet_delta_eta_cut    = 100.0;
static const Double_t jet_delta_phi_cut    = 45.0*(2.0*Pi/360.0);
static const Double_t track_pT_assoc_threshold = 2.0;
static const Double_t array_areas[N_jet_areas] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75};
static const Int_t N_Delta_pt_bins = 5;
static const Double_t array_pt_bins_Delta_pt[N_Delta_pt_bins] = {1,3,5,10,15};

// R = 0.2 --> area = 0.13
// R = 0.3 --> area = 0.28
// R = 0.4 --> area = 0.5
// R = 0.5 --> area = 0.79

static Double_t PYTHIA_hard_bin_weigh_factors[11]    = {1.6137644,0.136423,0.022972,0.005513,0.001650,0.000574,0.000390,0.000010200809,0.000000501270,0.000000028611,0.000000001446}; // in mb
static Double_t PYTHIA_hard_bin_high_pT_N_events[11] = {5,61,1275,9844,39541,103797,338451,1655573,3084417,3997019,4471887}; // Number of PYTHIA events with at least one 9 GeV/c track
static TH1D* h_PYTHIA_hard_bin_weigh_factors;
static TH1D* h_PYTHIA_hard_bin_high_pT_N_events;

static Long64_t N_PYTHIA_events[11];
static Long64_t counter_PYTHIA[11];

static TH1D* h_pT_eta_phi;
static TH2D* h_Delta_pt_vs_embed_pt[N_jet_areas][2][N_Psi_bins];
static TH2D* h_Delta_pt_vs_embed_pt_weight[N_jet_areas][2][N_Psi_bins];
static TH1D* h_jet_pt[N_jet_areas][N_leading_pt_bins+1][4];
static TH1D* h_jet_pt_sub[N_jet_areas][N_leading_pt_bins+1][4];
static TH2D* h2D_dijet_pt[N_jet_areas][N_leading_pt_bins+1][4];
static TH2D* h2D_dijet_pt_sub[N_jet_areas][N_leading_pt_bins+1][4];
static TH2D* h2D_dijet_pt_sub_test;
static TH2D* h2D_dijet_pt_sub_weight_test;
static TH1D* h_N_tracks_dijet[N_leading_pt_bins+1][2];
static TH2D* h2D_mult_vs_global_bin[N_leading_pt_bins+1][2];
static TH1D* h_trigger_track[4];
static TH2D* h_trigger_track_vs_global_bin[4];
static TH1D* h_jet_area[N_jet_areas][2];
static TH1D* h_jet_rho[N_jet_areas][2];
static TH1D* h_jet_per_event[N_jet_areas][2];
static TH1D* h_dijet_per_event[N_jet_areas][2];
static TProfile* p_jet_area_values;


static TH1D* h_jet_area_array[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH1D* h_jet_rho_array[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH1D* h_jet_rhoarea_array[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH2D* h2D_jet_rho_vs_mult_array[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH2D* h2D_track_eta_vs_phi[N_z_vertex_bins][N_mult_bins][N_Psi_bins][N_track_pt_bins_eta_phi];
static TH2D* h2D_jet_rho_vs_mult_array_In[N_z_vertex_bins][N_mult_bins][N_Psi_bins][2];
static TH1D* h_jet_Et_array[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH1D* h_jet_Et_array_weight[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH1D* h_jet_Et_array_In[N_z_vertex_bins][N_mult_bins][N_Psi_bins][2];
static TH1D* h_jet_per_event_array[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH1D* h_tracks_per_event_array[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH1D* h_tracks_per_event_array_In[N_z_vertex_bins][N_mult_bins][N_Psi_bins][2];
static TH2D* h_tracks_vs_z_vertex_array[N_z_vertex_bins][N_mult_bins];
static TH2D* h_Psi_vs_z_vertex_array[N_z_vertex_bins][N_Psi_bins];
static TH2D* h_tracks_vs_z_vertex_sample;
static TH1D* h_Psi2_sample;
static TH2D* h_jet_rho_vs_Et;
static TH2D* h_Phi_vs_eta[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH2D* h_Phi_vs_eta_random_phi[N_z_vertex_bins];
static TH1D* h_Momentum[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH2D* h_tracks_vs_z_vertex[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH1D* h_Psi2[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH1D* h_Et[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static TH2D* h_ratio_sub_lead_to_lead_pt_vs_lead_pt;
static TH2D* h_sub_lead_vs_lead_pt;
static TH2D* h_PsiA_vs_PsiB;
static TH1D* h_PsiA;
static TH1D* h_PsiB;
static TH1D* h_Psi_Full;
static TH1D* h_Psi_etapos;
static TH1D* h_Psi_etaneg;
static TH1D* h_phi;
static TH2D* h_area_vs_jet_pt;
static TH2D* h_area_vs_recoil_jet_pt[N_leading_pt_bins+1];
static TH1D* h_rho[N_leading_pt_bins+1];
static TH1D* h_area;
static TH1D* h_track_pT;
static TH1D* h_track_pT_piplus;
static TH1D* h_track_pT_cut;
static TH1D* h_tracks_above_threshold_per_event;

static TH1D* h_eta_recoil[4];
static TH1D* h_phi_recoil[4];
static TH1D* h_eta_trigger[4];
static TH1D* h_phi_trigger[4];
static TH1D* h_delta_phi_trigger_recoil[4];
static TH1D* h_Delta_phi_trig_jet;

static TProfile* p_pt_Ach[2];
static TH1F*     h_Ach;

static const Int_t N2D_tracks_above_threshold = 10;
static TH2D* h2D_tracks_above_threshold[N2D_tracks_above_threshold];

static TH1D* h_1D_det_diff_part_lvl_area;
static TH1D* h_1D_det_diff_part_lvl_area_lt5GeV;
static TH2D* h_2D_det_vs_part_lvl_area;
static TH2D* h_2D_det_vs_part_lvl_area_lt5GeV;
static TH1D* h_Sim_original_pT_cut[N_jet_areas];
static TH1D* h_Sim_original_pT_all[N_jet_areas];
static TH1D* h_Sim_smeared_pT_not_matched[N_jet_areas];
static TH2D* h2D_Sim_matched_pT_vs_original_pT[N_jet_areas];
static TH2D* h2D_Sim_original_pT_vs_matched_pT[N_jet_areas];
static TH1D* h_matched_tracks_fraction[N_jet_areas];
static TH2D* h2D_matched_tracks_fraction_vs_original_pT[N_jet_areas];
static TH1D* h_N_accepted_recoil_jets[N_jet_areas][N_leading_pt_bins+1];


static const Int_t N_m2_radial_dist = 10;
static TH2D* h2D_m2_vs_pT[N_m2_radial_dist][2]; // [radial distance][+,- charge]
static Double_t Array_m2_radial_dist[N_m2_radial_dist] = {1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1};

static const Int_t   Array_runID_eff_track_rec[7] = {12126100,12138025,12145021,12152017,12154022,12165032,12171017};
static const TString Array_PID_eff_track_rec[6]   = {"pplus","pminus","piplus","piminus","kplus","kminus"};

static TF1* f_EfficiencyVsPt[9][7][6];  //centrality(9),runID(7),PID(6)


static TProfile* p_v2_vs_pt;
static TProfile* p_v2_vs_pt_jet;
static TProfile* p_parameters;
static TProfile* p_Array_leading_pt_bins[2];
Float_t vertex_z_start_stop_delta[N_Beamtime][3] =
{
    {-40.0,40.0,1.0},
    {-40.0,40.0,1.0},
    {-40.0,40.0,1.0},
    {-40.0,40.0,1.0},
    {-40.0,40.0,1.0},
    {-40.0,40.0,1.0},
    {-40.0,40.0,1.0},
    {-40.0,40.0,1.0}
};

static const Double_t z_acceptance[N_Beamtime]  = {70.0,50.0,40.0,40.0,70.0,70.0,30.0,50.0};

Double_t Psi_start_stop_delta[N_Beamtime][3] =
{
    {-TMath::Pi()/2.0,TMath::Pi()/2.0,1.0},
    {-TMath::Pi()/2.0,TMath::Pi()/2.0,1.0},
    {-TMath::Pi()/2.0,TMath::Pi()/2.0,1.0},
    {-TMath::Pi()/2.0,TMath::Pi()/2.0,1.0},
    {-TMath::Pi()/2.0,TMath::Pi()/2.0,1.0},
    {-TMath::Pi()/2.0,TMath::Pi()/2.0,1.0},
    {-TMath::Pi()/2.0,TMath::Pi()/2.0,1.0},
    {-TMath::Pi()/2.0,TMath::Pi()/2.0,1.0}
};

//static const Float_t Array_leading_pt_bins[2][N_leading_pt_bins+1] =
//{
//    {0.0,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,10.0,15.0,0.0,5.0,10.0,7.0,0.0},
//    {0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,10.0,15.0,20.0,5.0,10.0,20.0,20.0,500.0}
//};

static const Float_t Array_leading_pt_bins[2][N_leading_pt_bins+1] =
{
    {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,0.0},
    {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,500.0}
};

static Double_t Integrated_events_SE_ME[2];
static Double_t SE_ME_integral_scale_factor;

static const Int_t n_z_vertex_bins = 10;
static const Int_t n_params = 3;
static TH1D* h_rc_QxQy_etapm_z_vs_index_[6][n_z_vertex_bins][n_params]; // stores the fit parameters as a function of the index
static TH1D* h_runId_index_rc;
static TH1D* h_runId_index_rc_inv;
static TH1D* h_runId_index_rc_clone;
static TFile* re_centering;
static Float_t start_z;
static Float_t delta_z;
static Int_t gBeamTimeNum;
static Int_t random_trigger_index;


//static TF1* LevyFit_pT;


//------------------------------------------------------------------------------------------------------------------
// For mixed event
Float_t mult_start_stop_delta[N_Beamtime][3] =
{
    {700.0,1000.0,1.0},
    {700.0,1000.0,1.0},
    {700.0,1000.0,1.0},
    {700.0,1000.0,1.0},
    {700.0,1000.0,1.0},
    {700.0,1000.0,1.0},
    {700.0,1000.0,1.0},
    {700.0,1000.0,1.0}
};

//OK


TFile* F_mixed_event[N_z_vertex_bins][N_mult_bins][N_Psi_bins];


//static StJetTrackEvent blubb;
//static StJetTrackEvent JetTrackEvent_Fill[N_z_vertex_bins][N_mult_bins][N_Psi_bins]; // -> Problem

#if 0

static StJetTrackEvent *JetTrackEvent_ptr_Fill[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static StJetTrackParticle  *JetTrackParticle_Fill;
static StJetTrackParticle  *JetTrackParticle_ME;
static TTree        *Tree_JetTrackEvent_Fill[N_z_vertex_bins][N_mult_bins][N_Psi_bins];
static const char JETTRACK_EVENT_TREE[]   = "JetTrackEvent";
static const char JETTRACK_EVENT_BRANCH[] = "Events";

static std::vector<StJetTrackEvent*> JetTrackEvent_ME;

static Int_t N_max_events;


//------------------------------------------------------------------------------------------------------------------


//static TFile* f_track_efficiencies;
//static TF1* func_effLow;
//static TF1* func_effHigh;


//TCanvas* c_3D;
//TH3D* h_3D_dummy;
static TNtuple *NT_ReCoil_Jet;
static Float_t ReCoil_Jet_NTDataArray[14];

static Int_t N_orig_smear = 1; // used only for PYTHIA



//------------------------------------------------------------------------------------------------------------------------------------
// bad run lists
static const Int_t n_bad_run_numbers[N_Beamtime] = {328,27,38,105,35,34,179,1};
static const Int_t bad_run_list_7GeV[328]   = {11114084,11114085,11114086,11114088,11114089,11114094,11114095,11114100,11114109,11115005,11115007,11115013,11115019,11115025,11115027,11115028,11115030,11115032,11115051,11115062,11115064,11115069,11115072,11115078,11115079,11115080,11115086,11115088,11115094,11116001,11116002,11116005,11116006,11116010,11116014,11116020,11116023,11116028,11116060,11116061,11116062,11116064,11116068,11116070,11116072,11116073,11116075,11117002,11117006,11117031,11117033,11117034,11117036,11117039,11117044,11117045,11117046,11117052,11117055,11117063,11117064,11117071,11117075,11117085,11117088,11117089,11117090,11117093,11117094,11117095,11117098,11117100,11117103,11117104,11117107,11118007,11118008,11118016,11118024,11118025,11118026,11118039,11118044,11119001,11119003,11119006,11119007,11119009,11119012,11119013,11119015,11119016,11119017,11119022,11119024,11119026,11119029,11119030,11119056,11119057,11119060,11119062,11119067,11119069,11119070,11119071,11119074,11119075,11119077,11119079,11119081,11119090,11119091,11119100,11119101,11120003,11120006,11120008,11120011,11120014,11120019,11120023,11120030,11120034,11120037,11120039,11120040,11120045,11120052,11120057,11120062,11120063,11120069,11120070,11120071,11120074,11120077,11120078,11120084,11120092,11121006,11121014,11121015,11121019,11121029,11121030,11121034,11121035,11121043,11121044,11121054,11121058,11121066,11121067,11121070,11121075,11121082,11122001,11122007,11122008,11122010,11122017,11122024,11122037,11122038,11122047,11122048,11122049,11122050,11122053,11122058,11122062,11122069,11122073,11122078,11122085,11122097,11123003,11123004,11123015,11123026,11123028,11123040,11123044,11123055,11123057,11123058,11123059,11123067,11123075,11123076,11123077,11123079,11123081,11123084,11123086,11123088,11123089,11123093,11123094,11123095,11123100,11123101,11123102,11123104,11124001,11124005,11124007,11124008,11124015,11124016,11124018,11124041,11124046,11124050,11124051,11124052,11124053,11124058,11124060,11124061,11124062,11124063,11124064,11124065,11124066,11124069,11125002,11125003,11125004,11125005,11125006,11125008,11125012,11125013,11125014,11125015,11125016,11125017,11125020,11125021,11125022,11125023,11125073,11125081,11125089,11125090,11125096,11125097,11126005,11126006,11126007,11126016,11126018,11126022,11126023,11127001,11127002,11127043,11128005,11128012,11128018,11128050,11128056,11128072,11129018,11129022,11129028,11129051,11130027,11130034,11130057,11131038,11131062,11132013,11132070,11133006,11133019,11134053,11134060,11134067,11134076,11135068,11136003,11136005,11136006,11136007,11136008,11136012,11136013,11136014,11136061,11136076,11136101,11136130,11136160,11136163,11137019,11138027,11138049,11138086,11138124,11139014,11140076,11140086,11141063,11142117,11143026,11143028,11144001,11144009,11144031,11144033,11144040,11144043,11144052,11145008,11145028,11145035,11146061,11146076,11146079,11147004,11147006,11147014,11147017,11147021,11147023};
static const Int_t bad_run_list_11GeV[27]   = {11148039,11148045,11149001,11149008,11149010,11149011,11149015,11149047,11150016,11150025,11150028,11151036,11151040,11151050,11152016,11152036,11152078,11153032,11153042,11155001,11155009,11156003,11156009,11157012,11158006,11158022,11158024};
static const Int_t bad_run_list_19GeV[35]   = {12113091,12114007,12114035,12114078,12114092,12114116,12115009,12115014,12115015,12115016,12115018,12115019,12115020,12115022,12115023,12115062,12115073,12115093,12115094,12116012,12116054,12117010,12117016,12117020,12117065,12119040,12119042,12120017,12120026,12121017,12121022,12121034,12121050,12121067,12122019};
static const Int_t bad_run_list_27GeV[34]   = {12172050,12172051,12172055,12173030,12173031,12173032,12173033,12173034,12174067,12174085,12175062,12175087,12175113,12175114,12175115,12176001,12176044,12176054,12176071,12177015,12177061,12177092,12177099,12177101,12177106,12177107,12177108,12178003,12178004,12178005,12178006,12178013,12178099,12178120};
static const Int_t bad_run_list_39GeV[38]   = {11199124,11100002,11100045,11101046,11102012,11102051,11102052,11102053,11102054,11102055,11102058,11103035,11103056,11103058,11103092,11103093,11105052,11105053,11105054,11105055,11107007,11107042,11107057,11107061,11107065,11107074,11108101,11109013,11109077,11109088,11109090,11109127,11110013,11110034,11110073,11110076,11111084,11111085};
static const Int_t bad_run_list_62GeV[105]  = {11080072,11081023,11081025,11082012,11082013,11082046,11082056,11082057,11084009,11084011,11084012,11084013,11084020,11084021,11084035,11084044,11084064,11085015,11085025,11085030,11085046,11085055,11085056,11085057,11086005,11086007,11087001,11087002,11087003,11087004,11088013,11089026,11089028,11089029,11089055,11089068,11089072,11091007,11091015,11091021,11091078,11092010,11092011,11092012,11092032,11092033,11092034,11092067,11092096,11093001,11094016,11094017,11094018,11094019,11094020,11094021,11094022,11094023,11094024,11094027,11094028,11094042,11094044,11094045,11094046,11094047,11094048,11094050,11094051,11094052,11094053,11094054,11094055,11094074,11094075,11094077,11095001,11095002,11095003,11095004,11095005,11095006,11095009,11095010,11095011,11095012,11095013,11095014,11095015,11095022,11095040,11095048,11095050,11095051,11095061,11095062,11095063,11095064,11095082,11095087,11096024,11096039,11096043,11096044,11097093};
static const Int_t bad_run_list_200GeV[179] = {12126101,12127003,12127017,12127018,12127032,12128025,12132043,12133018,12134023,12136005,12136006,12136014,12136017,12136022,12136023,12136024,12136025,12136027,12136028,12136029,12136030,12136031,12136034,12136054,12138017,12138021,12138081,12138082,12139006,12139007,12139015,12139016,12139028,12139059,12139075,12139076,12139077,12139078,12139079,12139080,12140010,12140011,12140012,12140013,12140014,12140015,12140016,12140018,12140019,12140020,12140021,12140025,12140026,12140027,12140028,12140029,12140042,12140051,12140052,12140053,12140054,12140055,12140056,12140064,12140066,12140067,12141001,12141002,12141003,12141004,12141005,12141006,12141009,12141014,12141015,12141016,12141017,12141018,12141019,12141026,12141027,12141028,12141029,12141030,12141032,12141033,12141034,12141035,12141036,12141041,12141042,12141043,12141044,12141045,12141046,12141048,12141050,12141051,12141052,12141056,12141059,12141060,12141061,12141062,12141063,12141064,12141065,12141066,12141067,12141071,12141072,12142001,12142002,12142003,12142006,12142013,12142014,12142015,12142016,12142017,12142018,12142019,12142020,12142021,12142022,12142023,12142026,12142027,12142032,12142033,12142034,12142046,12142047,12142048,12142049,12142050,12142051,12142061,12142062,12142063,12142076,12142077,12143016,12143018,12143054,12143075,12144001,12144002,12144013,12144014,12144027,12144028,12157038,12157051,12158040,12158041,12158054,12158056,12158057,12162055,12162056,12162057,12162058,12164037,12164078,12164079,12166002,12166003,12167015,12167024,12167052,12168002,12168009,12168022,12168077,12170044,12170045,12170054,12170056};
static const Int_t bad_run_list_15GeV[1]   = {1};
//------------------------------------------------------------------------------------------------------------------------------------
/*
#include "StJetAnalysis_Func.h"
*/


//------------------------------------------------------------------------------------------------------------------
StJetAnalysis::StJetAnalysis()
{

}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
StJetAnalysis::~StJetAnalysis()
{

}
//------------------------------------------------------------------------------------------------------------------



/*
//------------------------------------------------------------------------------------------------------------------
Int_t StJetAnalysis::Get_runID_range_Index(Int_t runID)
{
    for(Int_t i_range = 0; i_range < 6; i_range++)
    {
        if(runID > Array_runID_eff_track_rec[i_range] && runID <= Array_runID_eff_track_rec[i_range+1])
        {
            return i_range;
        }
    }

    return -1;
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
Double_t* StJetAnalysis::Parameter_eff_track_rec_function(Int_t Centrality, Int_t runID, Int_t PID, Int_t All)
{
    Int_t runID_range = -1;
    Double_t *Parameter = new Double_t[4];
    for(int i = 0; i < 4; i++)
    {
        Parameter[i] = 0;
    }

    runID_range = Get_runID_range_Index(runID);
    if(runID_range < 0)
    {
        cout << "WARNING: runID is out of range" << endl;
        return 0;
    }
    if( PID < 0 || PID > 5)
    {
        cout << "PID is out of range" << endl;
        return 0;
    }

    runID_range += 1; // index 0 is for all runIDs together
    if(All == 1){runID_range = 0;}

    Parameter[0] = Global_ABC_Parameters_A[Centrality][runID_range][PID];
    Parameter[1] = Global_ABC_Parameters_B[Centrality][runID_range][PID];
    Parameter[2] = Global_ABC_Parameters_C[Centrality][runID_range][PID];
    Parameter[3] = runID_range;

    return Parameter;
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void StJetAnalysis::Read_ABC_params_eff_track_rec_function()  //read ABC Parameter from a root file
{
    Int_t Switch_Array[6]={2,3,4,5,0,1}; //change PID in turn

    cout << "Read track reconstruction efficiency parameter file: " << eEff_file.Data() << endl;
    TFile *Eff_file=new TFile(eEff_file.Data());
    TH3D* hA  =(TH3D*)Eff_file->Get("hA");
    TH3D* hB  =(TH3D*)Eff_file->Get("hB");
    TH3D* hC  =(TH3D*)Eff_file->Get("hC");
    TH2D* hsA =(TH2D*)Eff_file->Get("hsA");
    TH2D* hsB =(TH2D*)Eff_file->Get("hsB");
    TH2D* hsC =(TH2D*)Eff_file->Get("hsC");

    for(Int_t i = 0; i < 6; i++) // PID
    {
        for(Int_t j = 0; j < 9; j++) // Centrality
        {
            for(Int_t k = 0; k < 7; k++) // runID
            {
                Global_ABC_Parameters_A[j][k][Switch_Array[i]] = hA->GetBinContent(hA->FindBin(i,j,k));
                Global_ABC_Parameters_B[j][k][Switch_Array[i]] = hB->GetBinContent(hB->FindBin(i,j,k));
                Global_ABC_Parameters_C[j][k][Switch_Array[i]] = hC->GetBinContent(hC->FindBin(i,j,k));
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------------
*/


//------------------------------------------------------------------------------------------------------------------
Int_t StJetAnalysis::Init()
{
    cout << "Init started" << endl;
    r3.SetSeed(0);

    //refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr(); // ALEX

    jet_R                      = eJet_R;
    jet_R_background           = eBkg_R;
    Remove_N_hardest           = eRem_n_hardest;
    max_pt_downscale_threshold = emax_pt_down_scale;


    for(Int_t i = 0; i < 11; i++)
    {
        counter_PYTHIA[i] = 0;
    }


    h_pT_eta_phi = new TH1D("h_pT_eta_phi","h_pT_eta_phi",100,0,1);

    for(Int_t i_charge = 0; i_charge < 2; i_charge++)
    {
        HistName = "p_pt_Ach_";
        HistName += i_charge;
        p_pt_Ach[i_charge] = new TProfile(HistName.Data(),HistName.Data(),15,-0.25,0.25);
    }
    h_Ach           = new TH1F("h_Ach","h_Ach",60,-0.3,0.3);

    /*
    //------------------------------------------------------------
    cout << "Initialize track reconstruction efficiency functions" << endl;
    TString FuncName;
    for(Int_t loop_Centr = 0; loop_Centr < 9; loop_Centr++) // Centrality
    {
        for(Int_t loop_runID = 0; loop_runID < 7; loop_runID++) // runID
        {
            for(Int_t loop_PID = 0; loop_PID < 6; loop_PID++) // PID
            {
                Int_t loop_runID_use = loop_runID - 1; // runIDs start from loop_runID = 1, for loop_runID = 0 all runIDs are used
                Int_t use_all_runIDs = 0;
                if(loop_runID == 0) // for first index use all runIDs
                {
                    loop_runID_use = 0;
                    use_all_runIDs = 1;
                }
                Double_t* Parameter = Parameter_eff_track_rec_function(loop_Centr,Array_runID_eff_track_rec[loop_runID_use+1],loop_PID,use_all_runIDs);

                FuncName =  "f_EfficiencyVsPt_runID_Centr_";
                FuncName += loop_Centr;
                FuncName += "_runID_";
                FuncName += loop_runID;
                FuncName += "_PID_";
                FuncName += Array_PID_eff_track_rec[loop_PID];

                f_EfficiencyVsPt[loop_Centr][loop_runID][loop_PID] = new TF1(FuncName.Data(),Eff_track_rec_function,0,50.0,3);
                f_EfficiencyVsPt[loop_Centr][loop_runID][loop_PID]->SetParameter(0,Parameter[0]);
                f_EfficiencyVsPt[loop_Centr][loop_runID][loop_PID]->SetParameter(1,Parameter[1]);
                f_EfficiencyVsPt[loop_Centr][loop_runID][loop_PID]->SetParameter(2,Parameter[2]);

                delete[] Parameter;
            }
        }
    }

    cout << "" << endl;
    cout << "-----------------------------------------------------" << endl;
    // cent9_eff: 0,1 : 0-10%, 7,8: 60-80%
    // PID_eff: p, anti-p, pi+, pi-, K+, K-
    TString PID_array[6] = {"p","anti-p","pi+","pi-","K+","K-"};
    Double_t pt_val_array_eff[3] = {0.2,0.4,20.0};
    for(Int_t iPID = 0; iPID < 6; iPID++)
    {
        cout << "" << endl;
        cout << "------------------------" << endl;
        cout << "PID: " << PID_array[iPID] << endl;
        for(Int_t i = 0; i < 3; i++)
        {
            Double_t eff_centA = f_EfficiencyVsPt[0][0][iPID]->Eval(pt_val_array_eff[i]);
            Double_t eff_centB = f_EfficiencyVsPt[1][0][iPID]->Eval(pt_val_array_eff[i]);

            Double_t eff_perA  = f_EfficiencyVsPt[7][0][iPID]->Eval(pt_val_array_eff[i]);
            Double_t eff_perB  = f_EfficiencyVsPt[8][0][iPID]->Eval(pt_val_array_eff[i]);

            Double_t eff_central    = (eff_centA+eff_centB)/2.0;
            Double_t eff_peripheral = (eff_perA+eff_perB)/2.0;

            cout << "pT: " << pt_val_array_eff[i] << " GeV/c, eff(0-10%): " << eff_central
                << ", eff(60%-80%): " << eff_peripheral << endl;
        }
    }
    cout << "-----------------------------------------------------" << endl;
    cout << "" << endl;
    //------------------------------------------------------------



    //------------------------------------------------------------
    cout << "Open re-centering correction histograms" << endl;
    // Open the mean Qx, Qy histograms for re-centering correction

    start_z = -z_acceptance[eBeamTimeNum];
    delta_z = 2.0*z_acceptance[eBeamTimeNum]/((Float_t)n_z_vertex_bins);
    gBeamTimeNum = eBeamTimeNum;

    re_centering = TFile::Open(re_centering_name.Data());  // open the file
    cout << "Opened re centering file = " << re_centering_name.Data() << endl;

    for(Int_t n_qxy = 0; n_qxy < 6; n_qxy++)
    {
        for(Int_t n_z = 0; n_z < n_z_vertex_bins; n_z++)
        {
            for(Int_t n_par = 0; n_par < n_params; n_par++)
            {
                HistName = "h_rc_QxQy_etapm_z_vs_index_";
                HistName += n_qxy;
                HistName += "_";
                HistName += n_z;
                HistName += "_";
                HistName += n_par;
                h_rc_QxQy_etapm_z_vs_index_[n_qxy][n_z][n_par] = (TH1D*)re_centering->FindObjectAny(HistName.Data());
            }
        }
    }

    h_runId_index_rc = (TH1D*)re_centering->FindObjectAny("h_runId_index");
    h_runId_index_rc->SetName("h_runId_index_rc");
    Double_t h_min_rc   = h_runId_index_rc->GetMinimum();
    Double_t h_max_rc   = h_runId_index_rc->GetMaximum();
    Int_t    h_Nbins_rc = h_runId_index_rc->GetNbinsX();
    Int_t    h_Nbins_rc_inv = (Int_t)(h_max_rc-h_min_rc)+1;
    cout << "h_Nbins_rc = " << h_Nbins_rc << ", h_min_rc = " << h_min_rc << ", h_max_rc = " << h_max_rc << endl;
    h_runId_index_rc_inv = new TH1D("h_runId_index_rc_inv","h_runId_index_rc_inv",h_Nbins_rc_inv,h_min_rc,h_max_rc+1.0);

    cout << "h_Nbins_rc_inv = " << h_runId_index_rc_inv->GetNbinsX() << endl;
    for(Int_t i_rc = 1; i_rc < h_Nbins_rc+1; i_rc++)
    {
        h_runId_index_rc_inv->SetBinContent(h_runId_index_rc_inv->FindBin(h_runId_index_rc->GetBinContent(i_rc)),(Float_t)i_rc);
    }
    h_runId_index_rc_clone = (TH1D*)h_runId_index_rc->Clone("h_runId_index_rc_clone");
    h_runId_index_rc_clone->Reset();
    for(Int_t i_rc = 1; i_rc < h_runId_index_rc_inv->GetNbinsX()+1; i_rc++)
    {
        h_runId_index_rc_clone->SetBinContent(h_runId_index_rc_clone->FindBin(h_runId_index_rc_inv->GetBinContent(i_rc)),h_runId_index_rc_inv->GetBinCenter(i_rc));
    }
    h_runId_index_rc_clone->Add(h_runId_index_rc,-1);
    // End recentering files
    //------------------------------------------------------------



    if(eMode == 311 || eMode == 312) N_orig_smear = 2; // PYTHIA

    if(eCentrality == 1) // 60-80%
    {
        mult_start_stop_delta[0][0] = 5.0;    // 5.0
        mult_start_stop_delta[0][1] = 120.0;  // 135.0
        mult_start_stop_delta[0][2] = 1.0;

        mult_start_stop_delta[1][0] = 5.0;
        mult_start_stop_delta[1][1] = 120.0;
        mult_start_stop_delta[1][2] = 1.0;

        mult_start_stop_delta[2][0] = 5.0;
        mult_start_stop_delta[2][1] = 120.0;
        mult_start_stop_delta[2][2] = 1.0;

        mult_start_stop_delta[3][0] = 5.0;
        mult_start_stop_delta[3][1] = 120.0;
        mult_start_stop_delta[3][2] = 1.0;

        mult_start_stop_delta[4][0] = 5.0;    // 5.0
        mult_start_stop_delta[4][1] = 120.0;  // 135.0
        mult_start_stop_delta[4][2] = 1.0;

        mult_start_stop_delta[5][0] = 5.0;
        mult_start_stop_delta[5][1] = 120.0;
        mult_start_stop_delta[5][2] = 1.0;

        mult_start_stop_delta[6][0] = 5.0;
        mult_start_stop_delta[6][1] = 120.0;
        mult_start_stop_delta[6][2] = 1.0;

        mult_start_stop_delta[7][0] = 5.0;
        mult_start_stop_delta[7][1] = 120.0;
        mult_start_stop_delta[7][2] = 1.0;
    }

    //cout << "Open track efficiency file" << endl;
    //f_track_efficiencies = TFile::Open("/project/projectdirs/star/aschmah/Jet/Data/Efficiencies/eff_pp.root");  // open the file
    //func_effLow  = (TF1*)f_track_efficiencies->Get("effhL");
    //func_effHigh = (TF1*)f_track_efficiencies->Get("effhH");
    */

    //----------------------------------------------------------------------------------------
    //cout << "Define functions" << endl;
    //LevyFit_pT            = new TF1("LevyFit_pT",LevyFitFunc_pT,0.0,6.0,4);
    //----------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------
    if(eMode == 310)
    {
        for(Int_t i_hist = 0; i_hist < N_m2_radial_dist; i_hist++)
        {
            for(Int_t i_charge = 0; i_charge < 2; i_charge++)
            {
                HistName = "h2D_m2_vs_pT_Rdist_";
                HistName += i_hist;
                HistName += "_charge_";
                HistName += i_charge;
                h2D_m2_vs_pT[i_hist][i_charge] = new TH2D(HistName.Data(),HistName.Data(),100,0,10.0,2000,-0.5,1.5);
            }
        }
    }

    if(eMode == 0 || eMode == 1 || eMode == 311 || eMode == 312)
    {
        // Same event input
        for(Int_t i_SE_ME = 0; i_SE_ME < 1; i_SE_ME++)
        {
            TString SE_ME_List = SEList;
            //if(i_SE_ME == 1) SE_ME_List = MEList;
            if (!SE_ME_List.IsNull())   // if input file is ok
            {
                cout << "Open file list " << SE_ME_List << endl;
                ifstream in(SE_ME_List);  // input stream
                if(in)
                {
                    cout << "file list is ok" << endl;
                    if(eMode == 0 || eMode == 1 || eMode == 311)
                    {
                        //cout << "Test A" << endl;
                        input_SE_ME[i_SE_ME]  = new TChain( ALEX_EVENT_TREE, ALEX_EVENT_TREE );
                        //cout << "Test B" << endl;
                        char str[255];       // char array for each file name
                        Long64_t entries_save = 0;
                        while(in)
                        {
                            in.getline(str,255);  // take the lines of the file list
                            if(str[0] != 0)
                            {
                                TString addfile;
                                if(eMode == 0 || eMode == 1) addfile = "jet_trees_V2/";
                                addfile += str;
                                if(eMode == 1) addfile = pinputdir+addfile;
                                if(eMode == 311) addfile = pinputdirPYTHIA+addfile;
                                if(eMode == 312) addfile = pinputdir+addfile;
                                Long64_t file_entries;
                                input_SE_ME[i_SE_ME] ->AddFile(addfile.Data(),-1, ALEX_EVENT_TREE );
                                file_entries = input_SE_ME[i_SE_ME]->GetEntries();
                                cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                                entries_save = file_entries;
                            }
                        }
                    }
                    if(eMode == 312)
                    {
                        char str[255];       // char array for each file name
                        Int_t PYTHIA_file = 0;
                        while(in)
                        {
                            input_PYTHIA[PYTHIA_file] = new TChain( ALEX_EVENT_TREE, ALEX_EVENT_TREE );
                            in.getline(str,255);  // take the lines of the file list
                            if(str[0] != 0)
                            {
                                TString addfile;
                                TString inputdir_full = pinputdirPYTHIA;
                                inputdir_full += "Pythia/Trees_9GeV/";
                                addfile += str;
                                addfile = inputdir_full+addfile;
                                Long64_t file_entries;
                                input_PYTHIA[PYTHIA_file] ->AddFile(addfile.Data(),-1, ALEX_EVENT_TREE );
                                file_entries = input_PYTHIA[PYTHIA_file]->GetEntries();
                                N_PYTHIA_events[PYTHIA_file] = file_entries;
                                cout << "File added to data chain: " << addfile.Data() << " with " << file_entries << " entries" << endl;
                            }
                            PYTHIA_file++;
                        }
                    }
                }
                else
                {
                    cout << "WARNING: file input is problemtic" << endl;
                }

                if(eMode == 0 || eMode == 1 || eMode == 311)
                {
                    JetTrackEvent = new StJetTrackEvent();
                    input_SE_ME[i_SE_ME]  ->SetBranchAddress( ALEX_EVENT_BRANCH, &JetTrackEvent );
                }
                if(eMode == 312)
                {
                    for(Int_t PYTHIA_file = 0; PYTHIA_file < 11; PYTHIA_file++)
                    {
                        cout << "Set branch address for PYTHIA hard bin file: " << PYTHIA_file << endl;
                        JetTrackEvent_PYTHIA[PYTHIA_file] = new StJetTrackEvent();
                        input_PYTHIA[PYTHIA_file]  ->SetBranchAddress( ALEX_EVENT_BRANCH, &JetTrackEvent_PYTHIA[PYTHIA_file] );
                    }
                }
            }

            if(eMode == 0 || eMode == 1 || eMode == 311)
            {
                file_entries_SE_ME[i_SE_ME] = input_SE_ME[i_SE_ME]->GetEntries();
                cout << "Number of entries in chain: " << file_entries_SE_ME[i_SE_ME] << endl;
            }
        }
    }
    //----------------------------------------------------------------------------------------------------


    //----------------------------------------------------------------------------------------------------
    if(eMode == 2 || eMode == 3 || eMode == 4
       || eMode == 11 || eMode == 31 || eMode == 310 || eMode == 32 || eMode == 42 || eMode == 312
      )
    {
        N_max_events = mult_start_stop_delta[eBeamTimeNum][1]*1.1;

        JetTrackEvent_ME.resize(N_max_events);
        for(Int_t mix_loop = 0; mix_loop < N_max_events; mix_loop++)
        {
            JetTrackEvent_ME[mix_loop] = new StJetTrackEvent();
        }
        //Int_t size = JetTrackEvent_ME.size();
        //cout << "size = " << size << endl;
        //JetTrackEvent_ME.clear();
        //cout << "cleared" << endl;

        input_SE_ME[0]  = new TChain( ALEX_EVENT_TREE, ALEX_EVENT_TREE );
        Long64_t entries_save = 0;

        TString addfile;

        //addfile = "histo_out/mode_";
        addfile = "/histo_out_V2/mode_";
        if(eIn_Mode != 24)
        {
            addfile += eIn_Mode;
        }
        else
        {
            addfile += 2;
        }
        addfile += "/F_mixed_event_z_";
        addfile += ez_bin;
        addfile += "_mult_";
        addfile += emult_bin;
        addfile += "_Psi_";
        addfile += ePsi_bin;
        addfile += "_mode_";
        if(eIn_Mode != 24)
        {
            addfile += eIn_Mode;
        }
        if(eIn_Mode == 2)
        {
            addfile += "_In_mode_1";
        }
        if(eIn_Mode == 24)
        {
            addfile += "2_In_mode_4";
        }
        if(eCentrality == 0) addfile += "_0_10";
        if(eCentrality == 1) addfile += "_60_80";
        addfile += "_V2.root";
        addfile = pinputdir+addfile; // <<--
        //addfile = pinputdir+"test_single_out_V67.root"; // REMOVE

        Long64_t file_entries;
        input_SE_ME[0] ->AddFile(addfile.Data(),-1, ALEX_EVENT_TREE );
        file_entries = input_SE_ME[0]->GetEntries();
        cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
        entries_save = file_entries;

        JetTrackEvent = new StJetTrackEvent();
        input_SE_ME[0]  ->SetBranchAddress( ALEX_EVENT_BRANCH, &JetTrackEvent );


        file_entries_SE_ME[0] = input_SE_ME[0]->GetEntries();
        cout << "Number of entries in chain: " << file_entries_SE_ME[0] << endl;

        if(eMode != 11)
        {
            cout << "Load sample histogram file: " << eSampleHist.Data() << endl;
            Inputfile = TFile::Open(eSampleHist.Data());  // open the file
            for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
            {
                for(Int_t i_mult = 0; i_mult < N_mult_bins; i_mult++)
                {
                    HistName = "h_tracks_vs_z_vertex_array_z_";
                    HistName += i_z;
                    HistName += "_m_";
                    HistName += i_mult;
                    h_tracks_vs_z_vertex_array[i_z][i_mult] = (TH2D*)Inputfile->Get(HistName.Data());
                }
            }
            for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
            {
                for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
                {
                    HistName = "h_Psi_vs_z_vertex_array_z_";
                    HistName += i_z;
                    HistName += "_m_";
                    HistName += i_Psi;
                    h_Psi_vs_z_vertex_array[i_z][i_Psi] = (TH2D*)Inputfile->Get(HistName.Data());
                }
            }

            for(Int_t i_pT = 0; i_pT < N_track_pt_bins_eta_phi; i_pT++)
            {
                HistName = "h2D_track_eta_vs_phi_";
                HistName += ez_bin;
                HistName += "_m_";
                HistName += emult_bin;
                HistName += "_P_";
                HistName += ePsi_bin;
                HistName += "_pT_";
                HistName += i_pT;
                h2D_track_eta_vs_phi[ez_bin][emult_bin][ePsi_bin][i_pT] = (TH2D*)Inputfile->Get(HistName.Data());
            }
        }
        cout << "Load SE Et histogram file: " << eSE_Et_Hist.Data() << endl;
        Inputfile_Et[0] = TFile::Open(eSE_Et_Hist.Data());  // open the file
        cout << "Load ME Et histogram file: " << eME_Et_Hist.Data() << endl;
        Inputfile_Et[1] = TFile::Open(eME_Et_Hist.Data());  // open the file
        //h_tracks_vs_z_vertex_sample = (TH2D*)Inputfile->Get("h_tracks_vs_z_vertex");
        //h_tracks_vs_z_vertex_sample ->SetName("h_tracks_vs_z_vertex_sample");
        //h_Psi2_sample               = (TH1D*)Inputfile->Get("h_Psi2");
        //h_Psi2_sample               ->SetName("h_Psi2_sample");
        cout << "Histograms loaded" << endl;
    }
    //----------------------------------------------------------------------------------------------------


    //----------------------------------------------------------------------------------------------------
    p_jet_area_values = new TProfile("p_jet_area_values","p_jet_area_values",100,0,100);
    for(Int_t i_val = 0; i_val < N_jet_areas; i_val++)
    {
        p_jet_area_values->Fill(i_val,(Double_t)array_areas[i_val]);
    }

    p_parameters = new TProfile("p_parameters","p_parameters",100,0,100);
    p_parameters ->Fill(0.0,(Double_t)N_Beamtime);
    p_parameters ->Fill(1.0,(Double_t)N_jet_areas);
    p_parameters ->Fill(2.0,(Double_t)N_z_vertex_bins);
    p_parameters ->Fill(3.0,(Double_t)N_mult_bins);
    p_parameters ->Fill(4.0,(Double_t)N_Psi_bins);
    p_parameters ->Fill(5.0,(Double_t)N_Et_bins);
    p_parameters ->Fill(6.0,(Double_t)N_leading_pt_bins);
    p_parameters ->Fill(7.0,(Double_t)leading_pt_split_cut);
    p_parameters ->Fill(8.0,(Double_t)leading_pt_split_val);
    p_parameters ->Fill(9.0,(Double_t)max_pt_threshold);
    p_parameters ->Fill(10.0,(Double_t)jet_delta_eta_cut);
    p_parameters ->Fill(11.0,(Double_t)jet_delta_phi_cut);
    p_parameters ->Fill(12.0,(Double_t)track_pT_assoc_threshold);
    p_parameters ->Fill(13.0,(Double_t)jet_R);
    p_parameters ->Fill(14.0,(Double_t)array_areas[0]);
    p_parameters ->Fill(15.0,(Double_t)array_areas[1]);
    p_parameters ->Fill(16.0,(Double_t)array_areas[2]);
    p_parameters ->Fill(17.0,(Double_t)array_areas[3]);
    p_parameters ->Fill(18.0,(Double_t)array_areas[4]);
    p_parameters ->Fill(19.0,(Double_t)Remove_N_hardest);
    p_parameters ->Fill(20.0,(Double_t)ePYTHIA_eff_factor);
    p_parameters ->Fill(21.0,(Double_t)downscale_factor);
    p_parameters ->Fill(22.0,(Double_t)N_Delta_pt_bins);
    p_parameters ->Fill(23.0,(Double_t)array_pt_bins_Delta_pt[0]);
    p_parameters ->Fill(24.0,(Double_t)array_pt_bins_Delta_pt[1]);
    p_parameters ->Fill(25.0,(Double_t)array_pt_bins_Delta_pt[2]);
    p_parameters ->Fill(26.0,(Double_t)array_pt_bins_Delta_pt[3]);
    p_parameters ->Fill(27.0,(Double_t)array_pt_bins_Delta_pt[4]);
    p_parameters ->Fill(28.0,(Double_t)jet_R_background);
    p_parameters ->Fill(29.0,(Double_t)track_eff_scaling_factor);
    p_parameters ->Fill(30.0,(Double_t)max_pt_downscale_threshold);
    p_parameters ->Fill(31.0,(Double_t)max_pt_val_embed);

    cout << "Parameters set" << endl;

    for(Int_t ipt_lead = 0; ipt_lead < (N_leading_pt_bins+1); ipt_lead++)
    {
        HistName = "h_area_vs_recoil_jet_pt_";
        HistName += ipt_lead;
        h_area_vs_recoil_jet_pt[ipt_lead] = new TH2D(HistName.Data(),HistName.Data(),400,-20,80,400,0,2);
        h_area_vs_recoil_jet_pt[ipt_lead]->Sumw2();

        HistName = "h_rho_";
        HistName += ipt_lead;
        h_rho[ipt_lead] = new TH1D(HistName.Data(),HistName.Data(),400,0.0,60);
        h_rho[ipt_lead]->Sumw2();
    }

    h_area_vs_jet_pt = new TH2D("h_area_vs_jet_pt","h_area_vs_jet_pt",400,-20,80,400,0,2);
    h_area_vs_jet_pt->Sumw2();
    h_area           = new TH1D("h_area","h_area",400,0,2);
    h_area->Sumw2();
    h_track_pT       = new TH1D("h_track_pT","h_track_pT",500,0,200);
    h_track_pT->Sumw2();
    h_track_pT_piplus       = new TH1D("h_track_pT_piplus","h_track_pT_piplus",500,0,200);
    h_track_pT_piplus->Sumw2();
    h_track_pT_cut   = new TH1D("h_track_pT_cut","h_track_pT_cut",500,0,200);
    h_track_pT_cut->Sumw2();
    h_tracks_above_threshold_per_event = new TH1D("h_tracks_above_threshold_per_event","h_tracks_above_threshold_per_event",20,0,20);
    h_tracks_above_threshold_per_event->Sumw2();

    h_Delta_phi_trig_jet = new TH1D("h_Delta_phi_trig_jet","h_Delta_phi_trig_jet",200,0,2.0*TMath::Pi());

    for(Int_t i_hist = 0; i_hist < N2D_tracks_above_threshold; i_hist++)
    {
        HistName = "h2D_tracks_above_threshold_";
        HistName += i_hist;
        h2D_tracks_above_threshold[i_hist] = new TH2D(HistName.Data(),HistName.Data(),20,0,20,20,0,20);
        h2D_tracks_above_threshold[i_hist]->Sumw2();
    }

    if(eMode == 0 || eMode == 3 || eMode == 31 || eMode == 310 || eMode == 32 || eMode == 311 || eMode == 312 || eMode == 42)
    {
        if(eMode == 3 || eMode == 31 || eMode == 310 || eMode == 32 || eMode == 312 || eMode == 42)
        {
            outputfile_name = eOutdir;
            outputfile_name += "F_mixed_event_z_";
            outputfile_name += ez_bin;
            outputfile_name += "_mult_";
            outputfile_name += emult_bin;
            outputfile_name += "_Psi_";
            outputfile_name += ePsi_bin;
            outputfile_name += "_mode_";
            outputfile_name += eMode;
            outputfile_name += "_In_mode_";
            outputfile_name += eIn_Mode;
            if(eCentrality == 0) outputfile_name += "_0_10";
            if(eCentrality == 1) outputfile_name += "_60_80";
            outputfile_name += ".root";
        }
        cout << "Output file: " << outputfile_name.Data() << endl;
        Outputfile        = new TFile(outputfile_name.Data(),"RECREATE");

        NT_ReCoil_Jet      = new TNtuple("NT_ReCoil_Jet","NT_ReCoil_Jet Ntuple","EventId:JetId:rho:area:Jetphi:Jeteta:Jetpt:TrackId:eta:phi:pt:x:y:z");
        NT_ReCoil_Jet      ->SetAutoSave( 5000000 );
        //c_3D       = new TCanvas("c_3D","c_3D",10,10,800,800);
        //h_3D_dummy = new TH3D("h_3D_dummy","h_3D_dummy",200,-300,300,200,-300,300,200,-1000,1000);
        //Draw_STAR_3D();
    }
    if(eMode == 11)
    {
        outputfile_name  = eOutdir;
        outputfile_name += "Jet_sample_histos_out_z_";
        //outputfile_name += eSuffix.Data();
        outputfile_name += ez_bin;
        outputfile_name += "_mult_";
        outputfile_name += emult_bin;
        outputfile_name += "_Psi_";
        outputfile_name += ePsi_bin;
        outputfile_name += "_mode_";
        outputfile_name += eMode;
        outputfile_name += "_In_mode_";
        outputfile_name += eIn_Mode;
        if(eCentrality == 0) outputfile_name += "_0_10";
        if(eCentrality == 1) outputfile_name += "_60_80";
        outputfile_name += ".root";
        Outputfile       = new TFile(outputfile_name.Data(),"RECREATE");
    }

    if(eMode == 1 || eMode == 2 || eMode == 4) // create sub category trees
    {
        cout << "Create sub category output files" << endl;

        Int_t start_z    = 0;
        Int_t start_mult = 0;
        Int_t start_Psi  = 0;
        //Int_t start_Psi  = ePsi_bin;
        Int_t stop_z     = N_z_vertex_bins;
        Int_t stop_mult  = N_mult_bins;
        Int_t stop_Psi   = N_Psi_bins;
        //Int_t stop_Psi   = ePsi_bin+1;
        if(eMode == 2
           || eMode == 1
           || eMode == 4
          )
        {
            start_z    = ez_bin;
            start_mult = emult_bin;
            start_Psi  = ePsi_bin;
            stop_z     = ez_bin+1;
            stop_mult  = emult_bin+1;
            stop_Psi   = ePsi_bin+1;
        }

        for(Int_t i_z = start_z; i_z < stop_z; i_z++)
        {
            for(Int_t i_mult = start_mult; i_mult < stop_mult; i_mult++)
            {
                for(Int_t i_Psi = start_Psi; i_Psi < stop_Psi; i_Psi++)
                {
                    TString mixed_event_outfile_name = eOutdir;
                    mixed_event_outfile_name += "F_mixed_event_z_";
                    mixed_event_outfile_name += i_z;
                    mixed_event_outfile_name += "_mult_";
                    mixed_event_outfile_name += i_mult;
                    mixed_event_outfile_name += "_Psi_";
                    mixed_event_outfile_name += i_Psi;
                    mixed_event_outfile_name += "_mode_";
                    mixed_event_outfile_name += eMode;
                    if(eMode == 2)
                    {
                        mixed_event_outfile_name += "_In_mode_";
                        mixed_event_outfile_name += eIn_Mode;
                    }
                    mixed_event_outfile_name += eSuffix.Data();
                    if(eCentrality == 0) mixed_event_outfile_name += "_0_10";
                    if(eCentrality == 1) mixed_event_outfile_name += "_60_80";
                    mixed_event_outfile_name += ".root";
                    F_mixed_event[i_z][i_mult][i_Psi] = new TFile(mixed_event_outfile_name.Data(),"RECREATE");

                    cout << "Sub category output file" << mixed_event_outfile_name.Data() << " created" << endl;

                    Tree_JetTrackEvent_Fill[i_z][i_mult][i_Psi] = NULL;
                    JetTrackEvent_ptr_Fill[i_z][i_mult][i_Psi] = &JetTrackEvent_Fill[i_z][i_mult][i_Psi];
                    Tree_JetTrackEvent_Fill[i_z][i_mult][i_Psi] = new TTree(JETTRACK_EVENT_TREE ,  JETTRACK_EVENT_TREE);
                    Tree_JetTrackEvent_Fill[i_z][i_mult][i_Psi]->Branch(JETTRACK_EVENT_BRANCH , "StJetTrackEvent", &JetTrackEvent_ptr_Fill[i_z][i_mult][i_Psi]);
                    Long64_t maxtreesize = 2000000000;
                    Tree_JetTrackEvent_Fill[i_z][i_mult][i_Psi]->SetMaxTreeSize(5*Long64_t(maxtreesize));
                    Tree_JetTrackEvent_Fill[i_z][i_mult][i_Psi]->SetAutoSave( 10000000 );
                    //Tree_JetTrackEvent_Fill[i_z][i_mult][i_Psi]->SetBasketSize("*",128000);
                    Tree_JetTrackEvent_Fill[i_z][i_mult][i_Psi]->SetBasketSize("*",128000*10);

                    HistName = "h_Momentum_z_";
                    HistName += i_z;
                    HistName += "_mult_";
                    HistName += i_mult;
                    HistName += "_Psi_";
                    HistName += i_Psi;
                    h_Momentum[i_z][i_mult][i_Psi]           = new TH1D(HistName.Data(),HistName.Data(),500,-50,50);

                    HistName = "h_tracks_vs_z_vertex_z_";
                    HistName += i_z;
                    HistName += "_mult_";
                    HistName += i_mult;
                    HistName += "_Psi_";
                    HistName += i_Psi;
                    h_tracks_vs_z_vertex[i_z][i_mult][i_Psi] = new TH2D(HistName.Data(),HistName.Data(),100,-50,50,1300,0,1300);

                    HistName = "h_Psi2_z_";
                    HistName += i_z;
                    HistName += "_mult_";
                    HistName += i_mult;
                    HistName += "_Psi_";
                    HistName += i_Psi;
                    h_Psi2[i_z][i_mult][i_Psi]               = new TH1D(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());

                    HistName = "h_Et_z_";
                    HistName += i_z;
                    HistName += "_mult_";
                    HistName += i_mult;
                    HistName += "_Psi_";
                    HistName += i_Psi;
                    h_Et[i_z][i_mult][i_Psi]                 = new TH1D(HistName.Data(),HistName.Data(),1000,0,1000);

                    HistName = "h_Phi_vs_eta_z_";
                    HistName += i_z;
                    HistName += "_mult_";
                    HistName += i_mult;
                    HistName += "_Psi_";
                    HistName += i_Psi;
                    h_Phi_vs_eta[i_z][i_mult][i_Psi]         = new TH2D(HistName.Data(),HistName.Data(),200,-1.5,1.5,200,-1.0*TMath::Pi(),1.0*TMath::Pi());
                }
            }
        }


        h_PsiA_vs_PsiB = new TH2D("h_PsiA_vs_PsiB","h_PsiA_vs_PsiB",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
        h_PsiA         = new TH1D("h_PsiA","h_PsiA",100,-TMath::Pi(),TMath::Pi());
        h_PsiB         = new TH1D("h_PsiB","h_PsiB",100,-TMath::Pi(),TMath::Pi());
        h_Psi_Full     = new TH1D("h_Psi_Full","h_Psi_Full",100,-TMath::Pi(),TMath::Pi());
    }

    /*

    for(Int_t i = 0; i < 4; i++)
    {
        HistName = "h_delta_phi_trigger_recoil_smear_orig_";
        HistName += i;
        h_delta_phi_trigger_recoil[i]  = new TH1D(HistName.Data(),HistName.Data(),100,-2.0*TMath::Pi(),2.0*TMath::Pi());

        HistName = "h_eta_recoil_smear_orig_";
        HistName += i;
        h_eta_recoil[i]  = new TH1D(HistName.Data(),HistName.Data(),100,-1.2,1.2);

        HistName = "h_phi_recoil_smear_orig_";
        HistName += i;
        h_phi_recoil[i]  = new TH1D(HistName.Data(),HistName.Data(),100,-2.0*TMath::Pi(),2.0*TMath::Pi());

        HistName = "h_eta_trigger_smear_orig_";
        HistName += i;
        h_eta_trigger[i] = new TH1D(HistName.Data(),HistName.Data(),100,-1.2,1.2);

        HistName = "h_phi_trigger_smear_orig_";
        HistName += i;
        h_phi_trigger[i] = new TH1D(HistName.Data(),HistName.Data(),100,-2.0*TMath::Pi(),2.0*TMath::Pi());
    }


    cout << "QA histograms created" << endl;
    //----------------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------
    if(eMode == 11)
    {
        for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
        {
            for(Int_t i_mult = 0; i_mult < N_mult_bins; i_mult++)
            {
                HistName = "h_tracks_vs_z_vertex_array_z_";
                HistName += i_z;
                HistName += "_m_";
                HistName += i_mult;
                h_tracks_vs_z_vertex_array[i_z][i_mult] = new TH2D(HistName.Data(),HistName.Data(),100,-50,50,1300,0,1300);

                for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
                {
                    for(Int_t i_pT = 0; i_pT < N_track_pt_bins_eta_phi; i_pT++)
                    {
                        HistName = "h2D_track_eta_vs_phi_";
                        HistName += i_z;
                        HistName += "_m_";
                        HistName += i_mult;
                        HistName += "_P_";
                        HistName += i_Psi;
                        HistName += "_pT_";
                        HistName += i_pT;
                        h2D_track_eta_vs_phi[i_z][i_mult][i_Psi][i_pT] = new TH2D(HistName.Data(),HistName.Data(),80,-TMath::Pi(),TMath::Pi(),80,-1,1);
                    }
                }
            }
        }
        for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
        {
            for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
            {
                HistName = "h_Psi_vs_z_vertex_array_z_";
                HistName += i_z;
                HistName += "_m_";
                HistName += i_Psi;
                h_Psi_vs_z_vertex_array[i_z][i_Psi] = new TH2D(HistName.Data(),HistName.Data(),100,-50,50,100,-TMath::Pi()/2.0,TMath::Pi()/2.0);
            }
        }



    }
    cout << "Additional histograms created" << endl;
    //----------------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------
    if( (eMode == 3 || eMode == 31 || eMode == 310 || eMode == 32 || eMode == 312 || eMode == 42) && eRandom != -1)
    {
        for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
        {
            for(Int_t i_mult = 0; i_mult < N_mult_bins; i_mult++)
            {
                for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
                {
                    Float_t Int_SE_ME[2];
                    for(Int_t iSE_ME = 0; iSE_ME < 2; iSE_ME++)
                    {

                        HistName = "QA_hist_arrays/h2D_jet_rho_vs_mult_array_z";
                        HistName += i_z;
                        HistName += "_m";
                        HistName += i_mult;
                        HistName += "_P";
                        HistName += i_Psi;
                        //cout << "Old name" << HistName.Data() << endl;
                        h2D_jet_rho_vs_mult_array_In[i_z][i_mult][i_Psi][iSE_ME] = (TH2D*)Inputfile_Et[iSE_ME]->Get(HistName.Data());
                        if(!h2D_jet_rho_vs_mult_array_In[i_z][i_mult][i_Psi][iSE_ME])
                        {
                            cout << "WARNING: " << HistName.Data() << " does not exist!" << endl;
                            if(eCentrality != 1)
                            {
                                h2D_jet_rho_vs_mult_array_In[i_z][i_mult][i_Psi][iSE_ME] = new TH2D(HistName.Data(),HistName.Data(),100,720,1020,150,15,45);
                            }
                            else
                            {
                                h2D_jet_rho_vs_mult_array_In[i_z][i_mult][i_Psi][iSE_ME] = new TH2D(HistName.Data(),HistName.Data(),60,0,180,150,0,30);
                            }
                        }
                        HistName += "_iSE_ME";
                        HistName += iSE_ME;
                        //cout << "New name" << HistName.Data() << endl;
                        h2D_jet_rho_vs_mult_array_In[i_z][i_mult][i_Psi][iSE_ME]->SetName(HistName.Data());
                        Integrated_events_SE_ME[iSE_ME] += h2D_jet_rho_vs_mult_array_In[i_z][i_mult][i_Psi][iSE_ME]->Integral(1,h2D_jet_rho_vs_mult_array_In[i_z][i_mult][i_Psi][iSE_ME]->GetNbinsX(),1,h2D_jet_rho_vs_mult_array_In[i_z][i_mult][i_Psi][iSE_ME]->GetNbinsY());

                    }
                    if(Int_SE_ME[0] > 0.0 && Int_SE_ME[1] > 0.0)
                    {
                        // Normalization? I removed it to take into account the different number of entries in SE and ME for different bins
                        //h_jet_Et_array_In[i_z][i_mult][i_Psi][1]->Scale(Int_SE_ME[0]/Int_SE_ME[1]);
                    }
                }
            }
        }
        cout << "All normalization histograms loaded" << endl;

        SE_ME_integral_scale_factor = 1.0;
        if(Integrated_events_SE_ME[0] > 0 && Integrated_events_SE_ME[1] > 0)
        {
            SE_ME_integral_scale_factor = Integrated_events_SE_ME[0]/Integrated_events_SE_ME[1];
        }
        cout << "Integrated_events_SE_ME[0] = " << Integrated_events_SE_ME[0] << ", Integrated_events_SE_ME[1] = " << Integrated_events_SE_ME[1] << ", SE_ME_integral_scale_factor = " << SE_ME_integral_scale_factor << endl;

        // Apply global normalization factor to ME histograms -> same integrated number of counts
        for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
        {
            for(Int_t i_mult = 0; i_mult < N_mult_bins; i_mult++)
            {
                for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
                {
                    //h_tracks_per_event_array_In[i_z][i_mult][i_Psi][1]->Scale(SE_ME_integral_scale_factor);
                    h2D_jet_rho_vs_mult_array_In[i_z][i_mult][i_Psi][1]->Scale(SE_ME_integral_scale_factor);
                }
            }
        }
    }
    //----------------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------
    h_ratio_sub_lead_to_lead_pt_vs_lead_pt = new TH2D("h_ratio_sub_lead_to_lead_pt_vs_lead_pt","h_ratio_sub_lead_to_lead_pt_vs_lead_pt",300,0,50,300,0,2);
    h_sub_lead_vs_lead_pt = new TH2D("h_sub_lead_vs_lead_pt","h_sub_lead_vs_lead_pt",300,0,50,300,0,50);

    for(Int_t i_area_acc = 0; i_area_acc < N_jet_areas; i_area_acc++)
    {
        for(Int_t ipt_lead = 0; ipt_lead < (N_leading_pt_bins+1); ipt_lead++)
        {
            HistName = "h_N_accepted_recoil_jets_area_";
            HistName += i_area_acc;
            HistName += "_pt_";
            HistName += ipt_lead;
            h_N_accepted_recoil_jets[i_area_acc][ipt_lead] = new TH1D(HistName.Data(),HistName.Data(),20,0,20);
            h_N_accepted_recoil_jets[i_area_acc][ipt_lead]->Sumw2();
        }
    }

    for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
    {
        for(Int_t i_mult = 0; i_mult < N_mult_bins; i_mult++)
        {
            for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
            {
                HistName = "h_jet_area_array_z";
                HistName += i_z;
                HistName += "_m";
                HistName += i_mult;
                HistName += "_P";
                HistName += i_Psi;
                h_jet_area_array[i_z][i_mult][i_Psi] = new TH1D(HistName.Data(),HistName.Data(),1000,0.0,2.0);

                HistName = "h_jet_rho_array_z";
                HistName += i_z;
                HistName += "_m";
                HistName += i_mult;
                HistName += "_P";
                HistName += i_Psi;
                h_jet_rho_array[i_z][i_mult][i_Psi] = new TH1D(HistName.Data(),HistName.Data(),1000,0.0,100);

                HistName = "h_jet_rhoarea_array_z";
                HistName += i_z;
                HistName += "_m";
                HistName += i_mult;
                HistName += "_P";
                HistName += i_Psi;
                h_jet_rhoarea_array[i_z][i_mult][i_Psi] = new TH1D(HistName.Data(),HistName.Data(),1000,0.0,200);

                HistName = "h2D_jet_rho_vs_mult_array_z";
                HistName += i_z;
                HistName += "_m";
                HistName += i_mult;
                HistName += "_P";
                HistName += i_Psi;
                if(eCentrality != 1)
                {
                    h2D_jet_rho_vs_mult_array[i_z][i_mult][i_Psi] = new TH2D(HistName.Data(),HistName.Data(),100,720,1020,150,15,45);
                }
                else
                {
                    h2D_jet_rho_vs_mult_array[i_z][i_mult][i_Psi] = new TH2D(HistName.Data(),HistName.Data(),60,0,180,150,0,30);
                }

                HistName = "h_jet_Et_array_z";
                HistName += i_z;
                HistName += "_m";
                HistName += i_mult;
                HistName += "_P";
                HistName += i_Psi;
                h_jet_Et_array[i_z][i_mult][i_Psi] = new TH1D(HistName.Data(),HistName.Data(),2000,0.0,1000);

                HistName = "h_jet_Et_array_weight_z";
                HistName += i_z;
                HistName += "_m";
                HistName += i_mult;
                HistName += "_P";
                HistName += i_Psi;
                h_jet_Et_array_weight[i_z][i_mult][i_Psi] = new TH1D(HistName.Data(),HistName.Data(),2000,0.0,1000);


                HistName = "h_jet_per_event_array_z";
                HistName += i_z;
                HistName += "_m";
                HistName += i_mult;
                HistName += "_P";
                HistName += i_Psi;
                h_jet_per_event_array[i_z][i_mult][i_Psi] = new TH1D(HistName.Data(),HistName.Data(),100,0.0,100);

                HistName = "h_tracks_per_event_array_z";
                HistName += i_z;
                HistName += "_m";
                HistName += i_mult;
                HistName += "_P";
                HistName += i_Psi;
                h_tracks_per_event_array[i_z][i_mult][i_Psi] = new TH1D(HistName.Data(),HistName.Data(),1300,0.0,1300);
            }
        }
    }

    h_jet_rho_vs_Et = new TH2D("h_jet_rho_vs_Et","h_jet_rho_vs_ET",2000,0.0,1000,1000,0.0,100);
    h_trigger_track[0] = new TH1D("h_trigger_track","h_trigger_track",N_leading_pt_bins+1,0,N_leading_pt_bins+1);
    h_trigger_track[1] = new TH1D("h_trigger_track_smear","h_trigger_track_smear",N_leading_pt_bins+1,0,N_leading_pt_bins+1);
    h_trigger_track[2] = new TH1D("h_trigger_track_all","h_trigger_track_all",N_leading_pt_bins+1,0,N_leading_pt_bins+1);
    h_trigger_track[3] = new TH1D("h_trigger_track_smear_all","h_trigger_track_smear_all",N_leading_pt_bins+1,0,N_leading_pt_bins+1);
    h_trigger_track_vs_global_bin[0] = new TH2D("h_trigger_track_vs_global_bin","h_trigger_track_vs_global_bin",N_global_bin,0,N_global_bin,N_leading_pt_bins+1,0,N_leading_pt_bins+1);
    h_trigger_track_vs_global_bin[1] = new TH2D("h_trigger_track_vs_global_bin_smear","h_trigger_track_vs_global_bin_smear",N_global_bin,0,N_global_bin,N_leading_pt_bins+1,0,N_leading_pt_bins+1);
    h_trigger_track_vs_global_bin[2] = new TH2D("h_trigger_track_vs_global_bin_all","h_trigger_track_vs_global_bin_all",N_global_bin,0,N_global_bin,N_leading_pt_bins+1,0,N_leading_pt_bins+1);
    h_trigger_track_vs_global_bin[3] = new TH2D("h_trigger_track_vs_global_bin_smear_all","h_trigger_track_vs_global_bin_smear_all",N_global_bin,0,N_global_bin,N_leading_pt_bins+1,0,N_leading_pt_bins+1);

    h_jet_rho_vs_Et                  ->Sumw2();
    h_trigger_track[0]               ->Sumw2();
    h_trigger_track[1]               ->Sumw2();
    h_trigger_track[2]               ->Sumw2();
    h_trigger_track[3]               ->Sumw2();
    h_trigger_track_vs_global_bin[0] ->Sumw2();
    h_trigger_track_vs_global_bin[1] ->Sumw2();
    h_trigger_track_vs_global_bin[2] ->Sumw2();
    h_trigger_track_vs_global_bin[3] ->Sumw2();


    for(Int_t i = 0; i < 2; i++)
    {
        HistName = "p_Array_leading_pt_bins_";
        HistName += i;
        p_Array_leading_pt_bins[i] = new TProfile(HistName.Data(),HistName.Data(),N_leading_pt_bins+1,0,N_leading_pt_bins+1);
        for(Int_t ipt_lead = 0; ipt_lead < (N_leading_pt_bins+1); ipt_lead++)
        {
            p_Array_leading_pt_bins[i] ->Fill(ipt_lead,Array_leading_pt_bins[i][ipt_lead]);
        }
    }

    for(Int_t ipt_lead = 0; ipt_lead < (N_leading_pt_bins+1); ipt_lead++)
    {
        for(Int_t i_orig_smear = 0; i_orig_smear < N_orig_smear; i_orig_smear++)
        {
            HistName = "h_N_tracks_dijet_pt_";
            HistName += ipt_lead;
            if(i_orig_smear == 1) HistName += "_smear";
            h_N_tracks_dijet[ipt_lead][i_orig_smear]  = new TH1D(HistName.Data(),HistName.Data(),100,mult_start_stop_delta[eBeamTimeNum][0],mult_start_stop_delta[eBeamTimeNum][1]);
            h_N_tracks_dijet[ipt_lead][i_orig_smear]  ->Sumw2();

            HistName = "h2D_mult_vs_global_bin_pt_";
            HistName += ipt_lead;
            if(i_orig_smear == 1) HistName += "_smear";
            h2D_mult_vs_global_bin[ipt_lead][i_orig_smear] = new TH2D(HistName.Data(),HistName.Data(),N_global_bin,0,N_global_bin,100,mult_start_stop_delta[eBeamTimeNum][0],mult_start_stop_delta[eBeamTimeNum][1]);
            h2D_mult_vs_global_bin[ipt_lead][i_orig_smear] ->Sumw2();
        }
    }

    if(eMode == 311 || eMode == 312)
    {
        HistName = "h_1D_det_diff_part_lvl_area";
        h_1D_det_diff_part_lvl_area = new TH1D(HistName.Data(),HistName.Data(),100,-0.6,0.6);

        HistName = "h_1D_det_diff_part_lvl_area_lt5GeV";
        h_1D_det_diff_part_lvl_area_lt5GeV = new TH1D(HistName.Data(),HistName.Data(),100,-0.6,0.6);

        HistName = "h_2D_det_vs_part_lvl_area";
        h_2D_det_vs_part_lvl_area = new TH2D(HistName.Data(),HistName.Data(),100,0,0.6,100,0,0.6);

        HistName = "h_2D_det_vs_part_lvl_area_lt5GeV";
        h_2D_det_vs_part_lvl_area_lt5GeV = new TH2D(HistName.Data(),HistName.Data(),100,0,0.6,100,0,0.6);
      

        for(Int_t i = 0; i < N_jet_areas; i++)
        {
            HistName = "h_Sim_original_pT_cut_area";
            HistName += i;
            h_Sim_original_pT_cut[i]          = new TH1D(HistName.Data(),HistName.Data(),280,0,70);
            HistName = "h_Sim_original_pT_all_area";
            HistName += i;
            h_Sim_original_pT_all[i]          = new TH1D(HistName.Data(),HistName.Data(),280,0,70);
            HistName = "h_Sim_smeared_pT_not_matched";
            HistName += i;
            h_Sim_smeared_pT_not_matched[i]          = new TH1D(HistName.Data(),HistName.Data(),280,0,70);
            HistName = "h2D_Sim_matched_pT_vs_original_pT_area";
            HistName += i;
            h2D_Sim_matched_pT_vs_original_pT[i]          = new TH2D(HistName.Data(),HistName.Data(),280,0,70,280,0,70);
            HistName = "h2D_Sim_original_pT_vs_matched_pT_area";
            HistName += i;
            h2D_Sim_original_pT_vs_matched_pT[i]          = new TH2D(HistName.Data(),HistName.Data(),280,0,70,280,0,70);
            HistName = "h_matched_tracks_fraction_area";
            HistName += i;
            h_matched_tracks_fraction[i]                  = new TH1D(HistName.Data(),HistName.Data(),150,0,1.5);
            HistName = "h2D_matched_tracks_fraction_vs_original_pT_area";
            HistName += i;
            h2D_matched_tracks_fraction_vs_original_pT[i] = new TH2D(HistName.Data(),HistName.Data(),280,0,70,150,0,1.5);

            h_1D_det_diff_part_lvl_area                   ->Sumw2();
            h_1D_det_diff_part_lvl_area_lt5GeV            ->Sumw2();
            h_2D_det_vs_part_lvl_area                     ->Sumw2();
            h_2D_det_vs_part_lvl_area_lt5GeV              ->Sumw2();
            h_Sim_original_pT_cut[i]                      ->Sumw2();
            h_Sim_original_pT_all[i]                      ->Sumw2();
            h_Sim_smeared_pT_not_matched[i]               ->Sumw2();
            h2D_Sim_matched_pT_vs_original_pT[i]          ->Sumw2();
            h2D_Sim_original_pT_vs_matched_pT[i]          ->Sumw2();
            h_matched_tracks_fraction[i]                  ->Sumw2();
            h2D_matched_tracks_fraction_vs_original_pT[i] ->Sumw2();
        }
 
    }

    if(eMode == 312)
    {
        h_PYTHIA_hard_bin_weigh_factors    = new TH1D("h_PYTHIA_hard_bin_weigh_factors","h_PYTHIA_hard_bin_weigh_factors",11,0,11);
        h_PYTHIA_hard_bin_high_pT_N_events = new TH1D("h_PYTHIA_hard_bin_high_pT_N_events","h_PYTHIA_hard_bin_high_pT_N_events",11,0,11);
        for(Int_t PYTHIA_file = 0; PYTHIA_file < 11; PYTHIA_file++)
        {
            h_PYTHIA_hard_bin_weigh_factors   ->SetBinContent(PYTHIA_file+1,PYTHIA_hard_bin_weigh_factors[PYTHIA_file]);
            h_PYTHIA_hard_bin_high_pT_N_events->SetBinContent(PYTHIA_file+1,PYTHIA_hard_bin_high_pT_N_events[PYTHIA_file]);

            //cout << "PYTHIA_file: " << PYTHIA_file << ", weight: " << PYTHIA_hard_bin_weigh_factors[PYTHIA_file]
            //    << ", weight2: " << h_PYTHIA_hard_bin_weigh_factors->GetBinContent(PYTHIA_file + 1) << endl;
        }
    }

    if(eMode == 42 || eMode == 312)
    {
        p_v2_vs_pt     = new TProfile("p_v2_vs_pt","p_v2_ps_pt",200,0,20);
        p_v2_vs_pt_jet = new TProfile("p_v2_vs_pt_jet","p_v2_ps_pt_jet",200,0,20);
        h_Psi_Full     = new TH1D("h_Psi_Full","h_Psi_Full",100,-TMath::Pi(),TMath::Pi());
        h_phi          = new TH1D("h_phi","h_phi",200,-2.0*TMath::Pi(),2.0*TMath::Pi());
        h_Psi_etapos   = new TH1D("h_Psi_etapos","h_Psi_etapos",100,-TMath::Pi(),TMath::Pi());
        h_Psi_etaneg   = new TH1D("h_Psi_etaneg","h_Psi_etaneg",100,-TMath::Pi(),TMath::Pi());

        for(Int_t i_orig_smear = 0; i_orig_smear < N_orig_smear; i_orig_smear++)
        {
            for(Int_t i = 0; i < N_jet_areas; i++)
            {
                for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
                {
                    HistName = "h_Delta_pt_vs_embed_pt_";
                    HistName += i;
                    HistName += "_P_";
                    HistName += i_Psi;
                    if(i_orig_smear == 1) HistName += "_smear";
                    const Int_t    N_bins_pt_embed = (Int_t)(max_pt_val_embed*4.0);
                    const Double_t Pt_range_embed  = (Double_t)(N_bins_pt_embed/4.0);
                    h_Delta_pt_vs_embed_pt[i][i_orig_smear][i_Psi] = new TH2D(HistName.Data(),HistName.Data(),N_bins_pt_embed,0,Pt_range_embed,160,-30,50.0);
                    h_Delta_pt_vs_embed_pt[i][i_orig_smear][i_Psi] ->Sumw2();

                    if(eMode == 42)
                    {
                        HistName = "h_Delta_pt_vs_embed_pt_weight_";
                        HistName += i;
                        HistName += "_P_";
                        HistName += i_Psi;
                        if(i_orig_smear == 1) HistName += "_smear";
                        const Int_t    N_bins_pt_embed = (Int_t)(max_pt_val_embed*4.0);
                        const Double_t Pt_range_embed  = (Double_t)(N_bins_pt_embed/4.0);
                        h_Delta_pt_vs_embed_pt_weight[i][i_orig_smear][i_Psi] = new TH2D(HistName.Data(),HistName.Data(),N_bins_pt_embed,0,Pt_range_embed,160,-30,50.0);
                        h_Delta_pt_vs_embed_pt_weight[i][i_orig_smear][i_Psi] ->Sumw2();
                    }
                }
            }
        }
    }

    HistName = "h2D_dijet_pt_sub_test";
    h2D_dijet_pt_sub_test = new TH2D(HistName.Data(),HistName.Data(),100,-0.5*Pi,1.5*Pi,350,-40,100.0);
    HistName = "h2D_dijet_pt_sub_weight_test";
    h2D_dijet_pt_sub_weight_test = new TH2D(HistName.Data(),HistName.Data(),100,-0.5*Pi,1.5*Pi,350,-40,100.0);
    for(Int_t i_orig_smear = 0; i_orig_smear < 4; i_orig_smear++)
    {
        for(Int_t i = 0; i < N_jet_areas; i++)
        {
            for(Int_t ipt_lead = 0; ipt_lead < (N_leading_pt_bins+1); ipt_lead++)
            {

                HistName = "h_jet_pt_";
                HistName += i;
                HistName += "_lead_pt_";
                HistName += ipt_lead;
                if(i_orig_smear == 1) HistName += "_smear";
                if(i_orig_smear > 1)
                {
                    HistName += "_";
                    HistName += i_orig_smear;
                }
                h_jet_pt[i][ipt_lead][i_orig_smear]        = new TH1D(HistName.Data(),HistName.Data(),1400,-40,100.0);
                h_jet_pt[i][ipt_lead][i_orig_smear]   ->Sumw2();

                HistName = "h_jet_pt_sub_";
                HistName += i;
                HistName += "_lead_pt_";
                HistName += ipt_lead;
                if(i_orig_smear == 1) HistName += "_smear";
                if(i_orig_smear > 1)
                {
                    HistName += "_";
                    HistName += i_orig_smear;
                }
                h_jet_pt_sub[i][ipt_lead][i_orig_smear]    = new TH1D(HistName.Data(),HistName.Data(),1400,-40,100.0);
                h_jet_pt_sub[i][ipt_lead][i_orig_smear]    ->Sumw2();

                HistName = "h2D_dijet_pt_";
                HistName += i;
                HistName += "_lead_pt_";
                HistName += ipt_lead;
                if(i_orig_smear == 1) HistName += "_smear";
                if(i_orig_smear > 1)
                {
                    HistName += "_";
                    HistName += i_orig_smear;
                }
                h2D_dijet_pt[i][ipt_lead][i_orig_smear]        = new TH2D(HistName.Data(),HistName.Data(),100,-0.5*Pi,1.5*Pi,350,-40,100.0);
                h2D_dijet_pt[i][ipt_lead][i_orig_smear]        ->Sumw2();

                HistName = "h2D_dijet_pt_sub_";
                HistName += i;
                HistName += "_lead_pt_";
                HistName += ipt_lead;
                if(i_orig_smear == 1) HistName += "_smear";
                if(i_orig_smear > 1)
                {
                    HistName += "_";
                    HistName += i_orig_smear;
                }
                h2D_dijet_pt_sub[i][ipt_lead][i_orig_smear]    = new TH2D(HistName.Data(),HistName.Data(),100,-0.5*Pi,1.5*Pi,350,-40,100.0);
                h2D_dijet_pt_sub[i][ipt_lead][i_orig_smear]    ->Sumw2();
            }

            if(i_orig_smear < N_orig_smear)
            {
                HistName = "h_jet_area_";
                HistName += i;
                if(i_orig_smear == 1) HistName += "_smear";
                h_jet_area[i][i_orig_smear]      = new TH1D(HistName.Data(),HistName.Data(),1000,0.0,2.0);
                h_jet_area[i][i_orig_smear]      ->Sumw2();

                HistName = "h_jet_rho_";
                HistName += i;
                if(i_orig_smear == 1) HistName += "_smear";
                h_jet_rho[i][i_orig_smear]      = new TH1D(HistName.Data(),HistName.Data(),1000,0.0,100);
                h_jet_rho[i][i_orig_smear]      ->Sumw2();

                HistName = "h_jet_per_event_";
                HistName += i;
                if(i_orig_smear == 1) HistName += "_smear";
                h_jet_per_event[i][i_orig_smear] = new TH1D(HistName.Data(),HistName.Data(),100,0.0,100.0);
                h_jet_per_event[i][i_orig_smear] ->Sumw2();

                HistName = "h_dijet_per_event_";
                HistName += i;
                if(i_orig_smear == 1) HistName += "_smear";
                h_dijet_per_event[i][i_orig_smear] = new TH1D(HistName.Data(),HistName.Data(),100,0.0,100.0);
                h_dijet_per_event[i][i_orig_smear] ->Sumw2();
            }
        }
    }

    if(eMode == 0 || eMode == 1 || eMode == 2 || eMode == 4)
    {
        for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
        {
            HistName = "h_Phi_vs_eta_random_phi_";
            HistName += i_z;
            h_Phi_vs_eta_random_phi[i_z] = new TH2D(HistName.Data(),HistName.Data(),200,-1.5,1.5,200,-1.0*TMath::Pi(),1.0*TMath::Pi());
   
        }
    }

    Float_t delta_z = (vertex_z_start_stop_delta[eBeamTimeNum][1] - vertex_z_start_stop_delta[eBeamTimeNum][0])/((Float_t)N_z_vertex_bins);
    vertex_z_start_stop_delta[eBeamTimeNum][2] = delta_z;

    Float_t delta_mult = (mult_start_stop_delta[eBeamTimeNum][1] - mult_start_stop_delta[eBeamTimeNum][0])/((Float_t)N_mult_bins);
    mult_start_stop_delta[eBeamTimeNum][2] = delta_mult;

    Float_t delta_Psi = (Psi_start_stop_delta[eBeamTimeNum][1] - Psi_start_stop_delta[eBeamTimeNum][0])/((Float_t)N_Psi_bins);
    Psi_start_stop_delta[eBeamTimeNum][2] = delta_Psi;
    //----------------------------------------------------------------------------------------------------
    */

    return 1;
}
//------------------------------------------------------------------------------------------------------------------
#endif



//------------------------------------------------------------------------------------------------------------------
Int_t StJetAnalysis::Make()
{
    cout << "Make started" << endl;
    r3.SetSeed(0);
    gRandom->SetSeed(0);
    cout << "Seed = " << r3.GetSeed() << endl;

    ran_gen.SetSeed(0);
    ran.SetSeed(0);



    //----------------------------------------------------------------------------------------------------
    cout << "Define bad run list" << endl;
    Int_t bad_run_numbers[n_bad_run_numbers[eBeamTimeNum]];
    for(Int_t j = 0; j < n_bad_run_numbers[eBeamTimeNum]; j++)
    {
        if(eBeamTimeNum == 0) bad_run_numbers[j] = bad_run_list_7GeV[j];
        if(eBeamTimeNum == 1) bad_run_numbers[j] = bad_run_list_11GeV[j];
        if(eBeamTimeNum == 2) bad_run_numbers[j] = bad_run_list_39GeV[j];
        if(eBeamTimeNum == 3) bad_run_numbers[j] = bad_run_list_62GeV[j];
        if(eBeamTimeNum == 4) bad_run_numbers[j] = bad_run_list_19GeV[j];
        if(eBeamTimeNum == 5) bad_run_numbers[j] = bad_run_list_27GeV[j];
        if(eBeamTimeNum == 6) bad_run_numbers[j] = bad_run_list_200GeV[j];
        if(eBeamTimeNum == 7) bad_run_numbers[j] = bad_run_list_15GeV[j];
    }
    //----------------------------------------------------------------------------------------------------



    if(eMode == 0 || eMode == 1 || eMode == 11 || eMode == 3 || eMode == 31 || eMode == 310 || eMode == 32 || eMode == 312 || eMode == 42 || eMode == 4 || eMode == 311)  // eMode 2 is mixed event
    {
        if(eMode == 1 || eMode == 2 || eMode == 4)
        {
            F_mixed_event[ez_bin][emult_bin][ePsi_bin] ->cd();
        }

        for(Int_t i_SE_ME = 0; i_SE_ME < 1; i_SE_ME++)
        {
            if(i_SE_ME == 0) cout << "Same event loop started" << endl;
            if(i_SE_ME == 1) cout << "Mixed event loop started" << endl;

            Int_t vertex_pos_counter = 0;

            Long64_t stop_event_use_loop = stop_event_use;
            if(stop_event_use_loop > file_entries_SE_ME[i_SE_ME]) stop_event_use_loop = file_entries_SE_ME[i_SE_ME];
            for(Long64_t counter = start_event_use; counter < stop_event_use_loop; counter++)
            {
                if (counter != 0  &&  counter % 100 == 0)
                    cout << "." << flush;
                if (counter != 0  &&  counter % 1000 == 0)
                {
                    if((stop_event_use_loop-start_event_use) > 0)
                    {
                        Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use_loop-start_event_use));
                        cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data " << flush;
                    }
                }

                if (!input_SE_ME[i_SE_ME]->GetEntry( counter )) // take the event -> information is stored in event
                    break;  // end of data chunk


                //-----------------------------------------------------------------------------
                // Event information
                Float_t  prim_vertex_x   = JetTrackEvent->getx();
                Float_t  prim_vertex_y   = JetTrackEvent->gety();
                Float_t  prim_vertex_z   = JetTrackEvent->getz();
                Int_t    RunId           = JetTrackEvent->getid();
                Float_t  refMult         = JetTrackEvent->getmult();
                Float_t  n_prim          = JetTrackEvent->getn_prim();
                Float_t  n_non_prim      = JetTrackEvent->getn_non_prim();
                Int_t    n_tofmatch_prim = JetTrackEvent->getn_tof_prim();
                Int_t    SE_ME_flag      = JetTrackEvent->getSE_ME_flag();
                Float_t  ZDCx            = JetTrackEvent->getZDCx();
                Float_t  BBCx            = JetTrackEvent->getBBCx();
                Float_t  vzVPD           = JetTrackEvent->getvzVpd();
                Int_t    N_Particles     = JetTrackEvent->getNumParticle();
                TVector2 QvecEtaPos      = JetTrackEvent->getQvecEtaPos();
                TVector2 QvecEtaNeg      = JetTrackEvent->getQvecEtaNeg();
                Int_t    cent9           = JetTrackEvent->getcent9();

                if(fabs(prim_vertex_z) > z_acceptance[eBeamTimeNum]) continue;

                Double_t reweight = 1.0;
                if(eMode != 311) // not PYTHIA
                {
                    //refmultCorrUtil->init(RunId); // ALEX
                    //refmultCorrUtil->initEvent(refMult, prim_vertex_z, ZDCx); // ALEX

                    // Get centrality bins
                    //   see StRefMultCorr.h for the definition of centrality bins
                    //erefMult_bin16   = refmultCorrUtil->getCentralityBin16();
                    //erefMult_bin     = refmultCorrUtil->getCentralityBin9();

                    //reweight = refmultCorrUtil->getWeight(); // ALEX
                    reweight = 1.0; // ALEX
                    //cout << "Centrality reweight: " << reweight << endl;
                }

                Int_t cent9_eff_array[9] = {8,7,6,5,4,3,2,1,0}; // index for centrality was inverted for efficiency functions
                Int_t cent9_eff = 0;
                if(cent9 >= 0 && cent9 < 9) cent9_eff = cent9_eff_array[cent9];

                Double_t Psi2 = 0.0;

                Double_t EP_eta_pos_corr    = -999.0; // angle
                Double_t EP_eta_neg_corr    = -999.0; // angle
                Double_t EP_Qx_eta_pos_corr = -999.0;
                Double_t EP_Qy_eta_pos_corr = -999.0;
                Double_t EP_Qx_eta_neg_corr = -999.0;
                Double_t EP_Qy_eta_neg_corr = -999.0;

                if(eMode == 1) // Do a re-centering correction for the Q-vectors
                {
                    Calc_Corr_EventPlane_Angle(JetTrackEvent, EP_eta_pos_corr, EP_eta_neg_corr, EP_Qx_eta_pos_corr,
                                               EP_Qy_eta_pos_corr, EP_Qx_eta_neg_corr, EP_Qy_eta_neg_corr);

                    QvecEtaPos.Set(EP_Qx_eta_pos_corr,EP_Qy_eta_pos_corr);
                    QvecEtaNeg.Set(EP_Qx_eta_neg_corr,EP_Qy_eta_neg_corr);


                    Double_t EP_eta_full = TMath::ATan2(EP_Qy_eta_pos_corr+EP_Qy_eta_neg_corr,EP_Qx_eta_pos_corr+EP_Qx_eta_neg_corr);
                    EP_eta_full /= 2.0;

                    Psi2 = EP_eta_full;
                }
                else
                {
                    EP_eta_pos_corr = TMath::ATan2(QvecEtaPos.Y(),QvecEtaPos.X());
                    EP_eta_pos_corr /= 2.0;
                    EP_eta_neg_corr = TMath::ATan2(QvecEtaNeg.Y(),QvecEtaNeg.X());
                    EP_eta_neg_corr /= 2.0;

                    Double_t EP_eta_full = TMath::ATan2(QvecEtaPos.Y()+QvecEtaNeg.Y(),QvecEtaPos.X()+QvecEtaNeg.X());
                    EP_eta_full /= 2.0;

                    Psi2 = EP_eta_full;
                }
                //-----------------------------------------------------------------------------



                //---------------------------------------------------------------------------
                // Good/Bad run selection
                Int_t flag_good_run = 1;
                for(Int_t bad_run = 0; bad_run < n_bad_run_numbers[eBeamTimeNum]; bad_run++)
                {
                    if(bad_run_numbers[bad_run] == (Int_t)RunId)
                    {
                        flag_good_run = 0;
                        break;
                    }
                }
                //if(!flag_good_run) cout << "bad run: " << (Int_t)RunId << endl;
                //---------------------------------------------------------------------------



                // Pileup protection
                if(
                   fabs(vzVPD - prim_vertex_z) < 3.0 &&
                   n_tofmatch_prim > 1 &&
                   flag_good_run
                  )
                {
                    //cout << "" << endl;
                    //cout << "New Event" << endl;


                    Int_t z_bin = -1;
                    if(prim_vertex_z > vertex_z_start_stop_delta[eBeamTimeNum][0] && prim_vertex_z < vertex_z_start_stop_delta[eBeamTimeNum][1])
                    {
                        z_bin = (Int_t)((prim_vertex_z-vertex_z_start_stop_delta[eBeamTimeNum][0])/vertex_z_start_stop_delta[eBeamTimeNum][2]);
                    }
                    Int_t mult_bin = -1;
                    if(N_Particles > mult_start_stop_delta[eBeamTimeNum][0] && N_Particles < mult_start_stop_delta[eBeamTimeNum][1])
                    {
                        mult_bin = (Int_t)((N_Particles-mult_start_stop_delta[eBeamTimeNum][0])/mult_start_stop_delta[eBeamTimeNum][2]);
                    }
                    Int_t Psi_bin = -1;
                    if(Psi2 > Psi_start_stop_delta[eBeamTimeNum][0] && Psi2 < Psi_start_stop_delta[eBeamTimeNum][1])
                    {
                        Psi_bin = (Int_t)((Psi2-Psi_start_stop_delta[eBeamTimeNum][0])/Psi_start_stop_delta[eBeamTimeNum][2]);
                    }

                    //if(mult_bin == 7 && Psi_bin == 0 && z_bin == 10)
                    //{
                    //cout << "N_Particles = " << N_Particles << ", refMult = " << refMult << ", mult_bin = " << mult_bin << ", Psi_bin = " << Psi_bin << ", z_bin = " << z_bin << endl;
                    //cout << "start mult: " << mult_start_stop_delta[eBeamTimeNum][0] << ", Delta mult: " << mult_start_stop_delta[eBeamTimeNum][2] << endl;
                    //QvecEtaPos.Print();
                    //QvecEtaNeg.Print();
                    //}

                    Int_t global_bin = z_bin + mult_bin*N_z_vertex_bins + Psi_bin*N_z_vertex_bins*N_mult_bins;
                    if(eMode != 1) global_bin = ez_bin + emult_bin*N_z_vertex_bins + ePsi_bin*N_z_vertex_bins*N_mult_bins;

                    if(eMode == 11
                       //&& z_bin != -1 && mult_bin != -1 && Psi_bin != -1
                      )
                    {
                        Int_t N_Particles_below_threshold = 0;
                        Long64_t N_tracks_above_threshold = 0;
                        Long64_t N_tracks_above_threshold_array[N2D_tracks_above_threshold];
                        for(Int_t i_hist = 0; i_hist < N2D_tracks_above_threshold; i_hist++)
                        {
                            N_tracks_above_threshold_array[i_hist] = 0;
                        }
                        Long64_t N_tracks_above_trigger   = 0;
                        for(Int_t i_Particle = 0; i_Particle < N_Particles; i_Particle++)
                        {
                            // Particle information
                            JetTrackParticle            = JetTrackEvent   ->getParticle(i_Particle);
                            Float_t dca                 = JetTrackParticle->get_dca_to_prim();
                            Float_t m2                  = JetTrackParticle->get_Particle_m2 ();
                            Float_t nSPi                = JetTrackParticle->get_Particle_nSigmaPi();
                            Float_t nSK                 = JetTrackParticle->get_Particle_nSigmaK();
                            Float_t nSP                 = JetTrackParticle->get_Particle_nSigmaP();
                            Float_t qp                  = JetTrackParticle->get_Particle_qp();
                            Float_t nhitsfit            = JetTrackParticle->get_Particle_hits_fit();
                            TLorentzVector TLV_Particle_prim = JetTrackParticle->get_TLV_Particle_prim();
                            TLorentzVector TLV_Particle_glob = JetTrackParticle->get_TLV_Particle_glob();

                            TLorentzVector TLV_Particle_use = TLV_Particle_prim;
                            if(eflab_prim_glob == 1) TLV_Particle_use = TLV_Particle_glob;

                            Double_t track_pT  = TLV_Particle_use.Pt();
                            if(track_pT != track_pT) continue; // that is a NaN test. It always fails if track_pT = nan.
                            Double_t track_eta = TLV_Particle_use.PseudoRapidity();
                            Double_t track_phi = TLV_Particle_use.Phi();
                            
                            if(dca < 1.0 && nhitsfit > 14)
                            {
                                h_track_pT_cut ->Fill(TLV_Particle_use.Pt());
                                if(TLV_Particle_use.Pt() > max_pt_threshold) N_tracks_above_threshold++;
                                for(Int_t i_hist = 0; i_hist < N2D_tracks_above_threshold; i_hist++)
                                {
                                    if(TLV_Particle_use.Pt() > (max_pt_threshold + 2.0*((Double_t)i_hist))) N_tracks_above_threshold_array[i_hist]++;
                                }
                                if(TLV_Particle_use.Pt() > 9.0 && TLV_Particle_use.Pt() <= max_pt_threshold) N_tracks_above_trigger++;
                            }
                     

                            Int_t epT_bin = 0;
                            for(Int_t i_pT = 0; i_pT < N_track_pt_bins_eta_phi; i_pT++)
                            {
                                if(track_pT <= array_pt_bins_eta_phi[i_pT])
                                {
                                    epT_bin = i_pT;
                                    break;
                                }
                            }


                            //if(nSPi == -999.0) cout << "i_Particle = " << i_Particle << ", mother_pT = " << nSK << ", pT = " << TLV_Particle_use.Pt() << endl;


                            if(
                               (nSPi == -999.0 && nSK < max_pt_threshold) ||  // nSK stores in eMode = 4 the orignal pT before splitting
                               (nSPi != -999.0 && TLV_Particle_use.Pt() < max_pt_threshold)
                              ) 
                            {
                                h2D_track_eta_vs_phi[ez_bin][emult_bin][ePsi_bin][epT_bin]->Fill(track_phi,track_eta);
                                if(epT_bin == 0)
                                {
                                    //cout << "track_pT: " << track_pT << ", eta = " << track_eta << ", phi = " << track_phi << endl;
                                    h_pT_eta_phi ->Fill(track_pT);
                                }
                                N_Particles_below_threshold++;
                            }


                            //else
                            //{
                            //    cout << "i_Particle = " << i_Particle << ", nSPi = " << nSPi << ", mother_pT = " << nSK << ", pT = " << TLV_Particle_use.Pt() << endl;
                            //}

                        }

                        h_tracks_above_threshold_per_event            ->Fill(N_tracks_above_threshold);
                        for(Int_t i_hist = 0; i_hist < N2D_tracks_above_threshold; i_hist++)
                        {
                            h2D_tracks_above_threshold[i_hist] ->Fill(N_tracks_above_trigger,N_tracks_above_threshold_array[i_hist]);
                        }
                        h_tracks_vs_z_vertex_array[ez_bin][emult_bin] ->Fill(prim_vertex_z,N_Particles_below_threshold);
                        h_Psi_vs_z_vertex_array[ez_bin][ePsi_bin]     ->Fill(prim_vertex_z,Psi2);

                    }

                    if(eMode == 0 && z_bin != -1)
                    {
                        //cout << "prim_vertex_z = " << prim_vertex_z << ", N_Particles = " << N_Particles << endl;

                        //h_tracks_vs_z_vertex->Fill(prim_vertex_z,N_Particles);
                        //h_Psi2 ->Fill(Psi2);

                        //-----------------------------------------------------------------------------
                        // Particle loop
                        for(Int_t i_Particle = 0; i_Particle < N_Particles; i_Particle++)
                        {
                            // Particle information
                            JetTrackParticle            = JetTrackEvent   ->getParticle(i_Particle);
                            Float_t dca                 = JetTrackParticle->get_dca_to_prim();
                            Float_t m2                  = JetTrackParticle->get_Particle_m2 ();
                            Float_t nSPi                = JetTrackParticle->get_Particle_nSigmaPi();
                            Float_t nSK                 = JetTrackParticle->get_Particle_nSigmaK();
                            Float_t nSP                 = JetTrackParticle->get_Particle_nSigmaP();
                            Float_t qp                  = JetTrackParticle->get_Particle_qp();
                            Float_t nhitsfit            = JetTrackParticle->get_Particle_hits_fit();
                            TLorentzVector TLV_Particle_prim = JetTrackParticle->get_TLV_Particle_prim();
                            TLorentzVector TLV_Particle_glob = JetTrackParticle->get_TLV_Particle_glob();

                            TLorentzVector TLV_Particle_use = TLV_Particle_prim;
                            if(eflab_prim_glob == 1) TLV_Particle_use = TLV_Particle_glob;

                            Double_t track_pT  = TLV_Particle_use.Pt();
                            if(track_pT != track_pT) continue; // that is a NaN test. It always fails if track_pT = nan.
                            Double_t track_eta = TLV_Particle_use.PseudoRapidity();
                            Double_t track_phi = TLV_Particle_use.Phi();

                            //h_Phi_vs_eta[z_bin] ->Fill(Eta,Phi);

                            Double_t ran_angle = ran_gen.Rndm()*TMath::Pi()*2.0;
                            TLV_Particle_use.SetPhi(ran_angle);
                            track_phi = TLV_Particle_use.Phi();
                            h_Phi_vs_eta_random_phi[z_bin] ->Fill(track_eta,track_phi);
                        }
                        //-----------------------------------------------------------------------------

                    }


                    //-----------------------------------------------------------------------------
                    if(eMode == 1)
                    {
                        // Fill event plane QA histograms
                        h_PsiA_vs_PsiB ->Fill(EP_eta_neg_corr,EP_eta_pos_corr);
                        h_PsiA         ->Fill(EP_eta_pos_corr);
                        h_PsiB         ->Fill(EP_eta_neg_corr);
                        h_Psi_Full     ->Fill(Psi2);
                    }
                    //-----------------------------------------------------------------------------


                    if(eMode == 1 &&
                       z_bin    == ez_bin &&
                       mult_bin == emult_bin &&
                       Psi_bin  == ePsi_bin

                       //z_bin    >= 0 && z_bin    < N_z_vertex_bins &&
                       //mult_bin >= 0 && mult_bin < N_mult_bins &&
                       //Psi_bin  >= 0 && Psi_bin  < N_Psi_bins
                       //Psi_bin  == ePsi_bin
                      )
                    {
                        //F_mixed_event[z_bin][mult_bin][Psi_bin] ->cd();

                        // Calculate biased Et
                        Double_t Et_total = 0.0;
                        for(Int_t i_Particle = 0; i_Particle < N_Particles; i_Particle++)
                        {
                            // Particle information
                            JetTrackParticle                 = JetTrackEvent->getParticle(i_Particle);
                            TLorentzVector TLV_Particle_prim = JetTrackParticle->get_TLV_Particle_prim();
                            TLorentzVector TLV_Particle_glob = JetTrackParticle->get_TLV_Particle_glob();

                            TLorentzVector TLV_Particle_use = TLV_Particle_prim;
                            if(eflab_prim_glob == 1) TLV_Particle_use = TLV_Particle_glob;

                            //if(TLV_Particle_use.Pt() < 2.5)
                            {
                                Et_total += TLV_Particle_use.Et();
                            }
                        }

                        h_Et[z_bin][mult_bin][Psi_bin]                 ->Fill(Et_total);
                        h_tracks_vs_z_vertex[z_bin][mult_bin][Psi_bin] ->Fill(prim_vertex_z,N_Particles);
                        h_Psi2[z_bin][mult_bin][Psi_bin]               ->Fill(Psi2);

                        // Fill event information for d4s
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].clearParticleList();
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setx(prim_vertex_x);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].sety(prim_vertex_y);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setz(prim_vertex_z);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setid(RunId);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setmult(refMult);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setn_prim(n_prim);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setn_non_prim(n_non_prim);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setn_tof_prim(n_tofmatch_prim);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setSE_ME_flag(SE_ME_flag);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setZDCx(ZDCx);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setBBCx(BBCx);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setvzVpd(vzVPD);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setQvecEtaPos(QvecEtaPos);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setQvecEtaNeg(QvecEtaNeg);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setcent9(cent9);

                        for(Int_t i_Particle = 0; i_Particle < N_Particles; i_Particle++)
                        {
                            // Particle information
                            JetTrackParticle            = JetTrackEvent   ->getParticle(i_Particle);
                            Float_t dca                 = JetTrackParticle->get_dca_to_prim();
                            Float_t m2                  = JetTrackParticle->get_Particle_m2 ();
                            Float_t nSPi                = JetTrackParticle->get_Particle_nSigmaPi();
                            Float_t nSK                 = JetTrackParticle->get_Particle_nSigmaK();
                            Float_t nSP                 = JetTrackParticle->get_Particle_nSigmaP();
                            Float_t qp                  = JetTrackParticle->get_Particle_qp();
                            Float_t nhitsfit            = JetTrackParticle->get_Particle_hits_fit();
                            TLorentzVector TLV_Particle_prim = JetTrackParticle->get_TLV_Particle_prim();
                            TLorentzVector TLV_Particle_glob = JetTrackParticle->get_TLV_Particle_glob();

                            TLorentzVector TLV_Particle_use = TLV_Particle_prim;
                            if(eflab_prim_glob == 1) TLV_Particle_use = TLV_Particle_glob;

                            Double_t track_pT  = TLV_Particle_use.Pt();
                            if(track_pT != track_pT) continue; // that is a NaN test. It always fails if track_pT = nan.
                        
                            Double_t Eta = TLV_Particle_use.PseudoRapidity();
                            Double_t Phi = TLV_Particle_use.Phi();

                            h_Phi_vs_eta[z_bin][mult_bin][Psi_bin]        ->Fill(Eta,Phi);
                            h_Momentum[z_bin][mult_bin][Psi_bin]          ->Fill(qp);

                            JetTrackParticle_Fill = JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].createParticle();
                            JetTrackParticle_Fill ->set_dca_to_prim(dca);
                            JetTrackParticle_Fill ->set_Particle_m2(m2);
                            JetTrackParticle_Fill ->set_Particle_nSigmaPi(nSPi);
                            JetTrackParticle_Fill ->set_Particle_nSigmaK(nSK);
                            JetTrackParticle_Fill ->set_Particle_nSigmaP(nSP);
                            JetTrackParticle_Fill ->set_Particle_qp(qp);
                            JetTrackParticle_Fill ->set_TLV_Particle_prim(TLV_Particle_prim);
                            JetTrackParticle_Fill ->set_TLV_Particle_glob(TLV_Particle_glob);
                            JetTrackParticle_Fill ->set_Particle_hits_fit(nhitsfit);
                        }

                        Tree_JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin]  ->Fill();
                    }


                    if(eMode == 4 &&
                       z_bin    == ez_bin &&
                       mult_bin == emult_bin &&
                       Psi_bin  == ePsi_bin
                      )
                    {
                        // Calculate biased Et
                        Double_t Et_total = 0.0;
                        for(Int_t i_Particle = 0; i_Particle < N_Particles; i_Particle++)
                        {
                            // Particle information
                            JetTrackParticle                 = JetTrackEvent->getParticle(i_Particle);
                            TLorentzVector TLV_Particle_prim = JetTrackParticle->get_TLV_Particle_prim();
                            TLorentzVector TLV_Particle_glob = JetTrackParticle->get_TLV_Particle_glob();

                            TLorentzVector TLV_Particle_use = TLV_Particle_prim;
                            if(eflab_prim_glob == 1) TLV_Particle_use = TLV_Particle_glob;

                            //if(TLV_Particle_use.Pt() < 2.5)
                            {
                                Et_total += TLV_Particle_use.Et();
                            }
                        }

                        h_Et[z_bin][mult_bin][Psi_bin]                 ->Fill(Et_total);
                        h_tracks_vs_z_vertex[z_bin][mult_bin][Psi_bin] ->Fill(prim_vertex_z,N_Particles);
                        h_Psi2[z_bin][mult_bin][Psi_bin]               ->Fill(Psi2);

                        // Fill event information for d4s
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].clearParticleList();
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setx(prim_vertex_x);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].sety(prim_vertex_y);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setz(prim_vertex_z);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setid(RunId);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setmult(refMult);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setn_prim(n_prim);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setn_non_prim(n_non_prim);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setn_tof_prim(n_tofmatch_prim);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setSE_ME_flag(SE_ME_flag);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setZDCx(ZDCx);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setBBCx(BBCx);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setvzVpd(vzVPD);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setQvecEtaPos(QvecEtaPos);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setQvecEtaNeg(QvecEtaNeg);
                        JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].setcent9(cent9);

                        for(Int_t i_Particle = 0; i_Particle < N_Particles; i_Particle++)
                        {
                            // Particle information
                            JetTrackParticle            = JetTrackEvent->getParticle(i_Particle);
                            Float_t dca                 = JetTrackParticle->get_dca_to_prim();
                            Float_t m2                  = JetTrackParticle->get_Particle_m2 ();
                            Float_t nSPi                = JetTrackParticle->get_Particle_nSigmaPi();
                            Float_t nSK                 = JetTrackParticle->get_Particle_nSigmaK();
                            Float_t nSP                 = JetTrackParticle->get_Particle_nSigmaP();
                            Float_t qp                  = JetTrackParticle->get_Particle_qp();
                            Float_t nhitsfit            = JetTrackParticle->get_Particle_hits_fit();
                            TLorentzVector TLV_Particle_prim = JetTrackParticle->get_TLV_Particle_prim();
                            TLorentzVector TLV_Particle_glob = JetTrackParticle->get_TLV_Particle_glob();

                            TLorentzVector TLV_Particle_use = TLV_Particle_prim;
                            if(eflab_prim_glob == 1) TLV_Particle_use = TLV_Particle_glob;

                            Double_t track_pT  = TLV_Particle_use.Pt();
                            if(track_pT != track_pT) continue; // that is a NaN test. It always fails if track_pT = nan.

                            Float_t Phi    = TLV_Particle_use.Phi();
                            Float_t Eta    = TLV_Particle_use.PseudoRapidity();
                            Float_t Energy = TLV_Particle_use.E();


                            h_Phi_vs_eta[z_bin][mult_bin][Psi_bin]        ->Fill(Eta,Phi);
                            h_Momentum[z_bin][mult_bin][Psi_bin]          ->Fill(qp);

                            Float_t Pt = track_pT;

                            //if(Pt > 10.0)
                            //{
                            //    cout << "Pt = " << Pt << ", Eta = " << Eta << ", Phi = " << Phi << ", m2 = " << m2
                            //        << ", nSPi = " << nSPi << ", qp = "  << qp << ", dca = " << dca << endl;
                            //}


                            if(Pt >= leading_pt_split_cut) // split the track into smaller pieces
                            {
                                Double_t N_d_split_tracks = Pt/leading_pt_split_val;
                                Int_t    N_i_split_tracks = (Int_t)N_d_split_tracks; // number of split tracks
                                Double_t pt_split_val     = Pt/((Double_t)N_i_split_tracks); // real pt value of split tracks

                                //cout << "Pt = " << Pt << ", Eta = " << Eta << ", Phi = " << Phi << ", N_i_split_tracks = " << N_i_split_tracks << ", pt_split_val = " << pt_split_val << endl;

                                for(Int_t i_split = 0; i_split < N_i_split_tracks; i_split++)
                                {
                                    TLorentzVector TLV_Particle_split_prim = TLV_Particle_use;
                                    Double_t Delta_eta = (ran_gen.Rndm()-0.5)*0.2;
                                    Double_t Delta_phi = (ran_gen.Rndm()-0.5)*0.2;
                                    Double_t Eta_split = TLV_Particle_use.PseudoRapidity() + Delta_eta;
                                    Double_t Phi_split = TLV_Particle_use.Phi() + Delta_phi;
                                    if(Eta_split > 1.2) Eta_split = 1.2 - ran_gen.Rndm()*0.1;
                                    if(Eta_split < -1.2) Eta_split = -1.2 + ran_gen.Rndm()*0.1;
                                    if(Phi_split > TMath::Pi()) Phi_split -= 2.0*TMath::Pi();
                                    if(Phi_split < (-TMath::Pi())) Phi_split += 2.0*TMath::Pi();
                                    Double_t E_split   = TLV_Particle_use.E()/((Double_t)N_d_split_tracks);
                                    TLV_Particle_split_prim.SetPtEtaPhiE(pt_split_val,Eta_split,Phi_split,E_split);

                                    TLorentzVector TLV_Particle_split_glob = TLV_Particle_glob;
                                    Delta_eta = (ran_gen.Rndm()-0.5)*0.2;
                                    Delta_phi = (ran_gen.Rndm()-0.5)*0.2;
                                    Eta_split = TLV_Particle_glob.PseudoRapidity() + Delta_eta;
                                    Phi_split = TLV_Particle_glob.Phi() + Delta_phi;
                                    if(Eta_split > 1.2) Eta_split = 1.2 - ran_gen.Rndm()*0.1;
                                    if(Eta_split < -1.2) Eta_split = -1.2 + ran_gen.Rndm()*0.1;
                                    if(Phi_split > TMath::Pi()) Phi_split -= 2.0*TMath::Pi();
                                    if(Phi_split < (-TMath::Pi())) Phi_split += 2.0*TMath::Pi();
                                    E_split   = TLV_Particle_glob.E()/((Double_t)N_d_split_tracks);
                                    TLV_Particle_split_glob.SetPtEtaPhiE(pt_split_val,Eta_split,Phi_split,E_split);

                                    //cout << "i_split = " << i_split << ", Eta_split = " << Eta_split << ", Phi_split = " << Phi_split << endl;

                                    JetTrackParticle_Fill = JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].createParticle();
                                    JetTrackParticle_Fill ->set_dca_to_prim(dca);
                                    JetTrackParticle_Fill ->set_Particle_m2(Phi);
                                    JetTrackParticle_Fill ->set_Particle_nSigmaPi(-999.0);
                                    JetTrackParticle_Fill ->set_Particle_nSigmaK(Pt);
                                    JetTrackParticle_Fill ->set_Particle_nSigmaP(N_d_split_tracks);
                                    JetTrackParticle_Fill ->set_Particle_qp(Eta);
                                    JetTrackParticle_Fill ->set_TLV_Particle_prim(TLV_Particle_split_prim);
                                    JetTrackParticle_Fill ->set_TLV_Particle_glob(TLV_Particle_split_glob);
                                    JetTrackParticle_Fill ->set_Particle_hits_fit(nhitsfit);
                                }
                            }
                            else
                            {
                                JetTrackParticle_Fill = JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin].createParticle();
                                JetTrackParticle_Fill ->set_dca_to_prim(dca);
                                JetTrackParticle_Fill ->set_Particle_m2(m2);
                                JetTrackParticle_Fill ->set_Particle_nSigmaPi(nSPi);
                                JetTrackParticle_Fill ->set_Particle_nSigmaK(nSK);
                                JetTrackParticle_Fill ->set_Particle_nSigmaP(nSP);
                                JetTrackParticle_Fill ->set_Particle_qp(qp);
                                JetTrackParticle_Fill ->set_TLV_Particle_prim(TLV_Particle_prim);
                                JetTrackParticle_Fill ->set_TLV_Particle_glob(TLV_Particle_glob);
                                JetTrackParticle_Fill ->set_Particle_hits_fit(nhitsfit);
                            }
                        }

                        Tree_JetTrackEvent_Fill[z_bin][mult_bin][Psi_bin]  ->Fill();
                    }


                    if(eMode == 3 || eMode == 31 || eMode == 310 || eMode == 32 || eMode == 312 || eMode == 311 || eMode == 42)
                    {
                        //-----------------------------------------------------------------------------
                        // Particle loop
                        std::vector< std::vector<PseudoJet> > particles; // original, smeared for PYTHIA
                        std::vector< std::vector<Double_t> > particle_info;
                        particle_info.resize(9);
                        particles.resize(2);
                        vector<Float_t>   particles_info;
                        Double_t Et_total = 0.0;
                        Int_t i_Particle_use = 0;

                        Double_t pt_Delta_pt = 0.0;
                        Double_t qp_Embed    = 1.0;
                        if(eMode == 42) // Get embedding pt value for Delta pt calculation
                        {
                            //Int_t pt_bin_Delta_pt = ran_gen.Integer(N_Delta_pt_bins);
                            //pt_Delta_pt = array_pt_bins_Delta_pt[pt_bin_Delta_pt];
                            pt_Delta_pt = ran_gen.Rndm()*max_pt_val_embed;
                        }


                        std::vector< std::vector<Double_t> > vec_trigger_tracks; // eta, phi, pT, index

                        Long64_t N_tracks_above_threshold = 0;
                        Long64_t N_tracks_above_threshold_array[N2D_tracks_above_threshold];
                        for(Int_t i_hist = 0; i_hist < N2D_tracks_above_threshold; i_hist++)
                        {
                            N_tracks_above_threshold_array[i_hist] = 0;
                        }
                        Long64_t N_tracks_above_trigger   = 0;

                        Int_t N_pos_neg_charges_reconstructed[2]     = {0,0};
                        Int_t N_pos_neg_charges_reconstructed_Ach[2] = {0,0};

                        for(Int_t i_Particle = 0; i_Particle < N_Particles; i_Particle++)
                        {
                            // Particle information
                            JetTrackParticle            = JetTrackEvent->getParticle(i_Particle);
                            Float_t dca                 = JetTrackParticle->get_dca_to_prim();
                            Float_t m2                  = JetTrackParticle->get_Particle_m2 ();
                            Float_t nSPi                = JetTrackParticle->get_Particle_nSigmaPi();
                            Float_t nSK                 = JetTrackParticle->get_Particle_nSigmaK();
                            Float_t nSP                 = JetTrackParticle->get_Particle_nSigmaP();
                            Float_t qp                  = JetTrackParticle->get_Particle_qp();
                            Float_t nhitsfit            = JetTrackParticle->get_Particle_hits_fit();
                            TLorentzVector TLV_Particle_prim = JetTrackParticle->get_TLV_Particle_prim();
                            TLorentzVector TLV_Particle_glob = JetTrackParticle->get_TLV_Particle_glob();

                            TLorentzVector TLV_Particle_use = TLV_Particle_prim;
                            if(eflab_prim_glob == 1) TLV_Particle_use = TLV_Particle_glob;

                            Double_t track_phi    = TLV_Particle_use.Phi();
                            Double_t track_pT     = TLV_Particle_use.Pt();

                            if(track_pT != track_pT) continue; // that is a NaN test. It always fails if track_pT = nan.
                            Double_t track_eta    = TLV_Particle_use.PseudoRapidity();

                            if(nhitsfit < 15) cout << "track_eta: " << track_eta << ", track_pT: " << track_pT << ", track_phi: " << track_phi << ", nhitsfit: " << nhitsfit << endl;
                            h_track_pT            ->Fill(TLV_Particle_use.Pt());
                            if(m2 > 0.01 && m2 < 0.03 && qp > 0.0) // fill only pi+
                            {
                                h_track_pT_piplus     ->Fill(TLV_Particle_use.Pt());
                            }


#if 0
                            //---------------------------------------------------------------
                            // A_ch analysis -> CMW effect (nothing to do with jets)
                            // Cuts to determine Ach
                            Int_t charge = 0;
                            if(qp < 0.0) charge = 1;

                            if(
                               fabs(track_eta) < 1.0
                               && track_pT > 0.15
                               && track_pT < 12.0
                               && !(track_pT < 0.4 && fabs(nSP) < 3.0)
                               && dca < 1.0
                               && nhitsfit >= 15
                               && N_pos_neg_charges_reconstructed_Ach[charge] < max_tracks
                              )
                            {
                                N_pos_neg_charges_reconstructed_Ach[charge]++;
                            }

                            // Cuts to determine pion pT
                            if(
                               fabs(track_eta) < 1.0
                               && track_pT > 0.15
                               && track_pT < 0.5
                               && dca < 1.0
                               && nhitsfit >= 15
                               && N_pos_neg_charges_reconstructed[charge] < max_tracks
                               && fabs(nSPi) < 2.0
                              )
                            {
                                pt_tracks[charge][N_pos_neg_charges_reconstructed[charge]] = track_pT;
                                N_pos_neg_charges_reconstructed[charge]++;
                            }
                            //---------------------------------------------------------------

#endif

                            if(dca < 1.0 && nhitsfit > 14)
                            {
                                h_track_pT_cut ->Fill(TLV_Particle_use.Pt());
                                if(TLV_Particle_use.Pt() > max_pt_threshold) N_tracks_above_threshold++;
                                for(Int_t i_hist = 0; i_hist < N2D_tracks_above_threshold; i_hist++)
                                {
                                    if(TLV_Particle_use.Pt() > (max_pt_threshold + 2.0*((Double_t)i_hist))) N_tracks_above_threshold_array[i_hist]++;
                                }
                                if(TLV_Particle_use.Pt() > 9.0 && TLV_Particle_use.Pt() <= max_pt_threshold) N_tracks_above_trigger++;


                                if(EP_eta_pos_corr > -999.0 && EP_eta_neg_corr > -999.0)
                                {
                                    Double_t Psi2_use = EP_eta_pos_corr;
                                    if(track_eta > 0.0) Psi2_use = EP_eta_neg_corr;
                                    Double_t phi_use  = track_phi;

                                    // -pi/2..pi/2 -> 0..pi
                                    if(Psi2_use < 0.0) Psi2_use += TMath::Pi();

                                    // -pi..pi -> 0..pi
                                    if(phi_use < 0.0) phi_use = phi_use + TMath::Pi();

                                    // -pi..pi, delta_phi_angle: 0..pi
                                    Double_t delta_phi_angle = phi_use - Psi2_use;
                                    if(phi_use >= 0.0 && delta_phi_angle >= 0.0) delta_phi_angle = delta_phi_angle;
                                    if(phi_use >= 0.0 && delta_phi_angle < 0.0)  delta_phi_angle += TMath::Pi();
                                    if(phi_use < 0.0  && delta_phi_angle >= -TMath::Pi()) delta_phi_angle += TMath::Pi();
                                    if(phi_use < 0.0  && delta_phi_angle < -TMath::Pi())  delta_phi_angle += 2.0*TMath::Pi();

                                    Double_t v2_val = TMath::Cos(2.0*delta_phi_angle);

                                    //--------------------------------------------------------
                                    if(eMode == 42)
                                    {
                                        p_v2_vs_pt ->Fill(track_pT,v2_val);
                                    }
                                }

                            }
                            // if(fabs(track_eta) > 1.0) cout << "i_Particle: " << i_Particle << ", Track pT: " << TLV_Particle_use.Pt() << ", Track eta: " << track_eta << endl;

                            qp_Embed = qp;

                            //if(eRandom == 1)
                            //{
                            //Double_t ran_angle = ran_gen.Rndm()*TMath::Pi()*2.0;
                            //TLV_Particle_use.SetPhi(ran_angle);
                            //}

                            if(
                               (nSPi == -999.0 && nSK < max_pt_threshold) ||  // nSK stores in eMode = 4 the orignal pT before splitting
                               (nSPi == -999.0 && eIn_Mode == 24) ||
                               (nSPi != -999.0 && TLV_Particle_use.Pt() < max_pt_threshold) // remove tracks above max pT threshold
                               //(nSPi != -999.0)
                              )
                            {
                                if(dca < 1.0)
                                {
                                    //----------------------------------------------
                                    // Trigger particles
                                    if(track_pT > 0.2 && track_pT < max_pt_threshold && fabs(track_eta) < 1.0)
                                    {
                                        std::vector<Double_t> vec_in;
                                        vec_in.resize(4);
                                        vec_in[0] = track_eta;
                                        vec_in[1] = track_phi;
                                        vec_in[2] = track_pT;
                                        vec_in[3] = i_Particle_use;

                                        vec_trigger_tracks.push_back(vec_in);
                                    }
                                    //----------------------------------------------


                                    if((eMode == 32 || eMode == 312) && TLV_Particle_use.Pt() > max_pt_downscale_threshold)
                                    {
                                        //cout << "Before: p = {" << TLV_Particle_use.Px() << ", " <<  TLV_Particle_use.Py() << ", " << TLV_Particle_use.Pz() << ", " << TLV_Particle_use.E() << "}" << endl;
                                        TLV_Particle_use *= downscale_factor;
                                        //cout << "After: p = {" << TLV_Particle_use.Px() << ", " <<  TLV_Particle_use.Py() << ", " << TLV_Particle_use.Pz() << ", " << TLV_Particle_use.E() << "}" << endl;
                                    }
                                    PseudoJet Fill_PseudoJet(TLV_Particle_use.Px(),TLV_Particle_use.Py(),TLV_Particle_use.Pz(),TLV_Particle_use.E());
                                    Fill_PseudoJet.set_user_index(i_Particle_use);
                                    particle_info[0].push_back(m2);
                                    particle_info[1].push_back(dca);
                                    particle_info[2].push_back(nSPi);
                                    particle_info[3].push_back(nSK);
                                    particle_info[4].push_back(nSP);
                                    particle_info[5].push_back(qp);
                                    particle_info[6].push_back(nhitsfit);
                                    particle_info[7].push_back(TLV_Particle_use.Eta());
                                    particle_info[8].push_back(TLV_Particle_use.Phi());

                                    //particles.push_back( PseudoJet(TLV_Particle_use.Px(),TLV_Particle_use.Py(),TLV_Particle_use.Pz(),TLV_Particle_use.E()) );
                                    if(eMode != 312) particles[0].push_back(Fill_PseudoJet); // mode 312: save only original PYTHIA tracks (later) for hparticles[0] -> Calculate Delta pT with PYTHIA jets

                                    if(eMode == 312)
                                    {
                                        particles[1].push_back(Fill_PseudoJet); // ME
                                    }

                                    particles_info.push_back(qp);
                                    Et_total += TLV_Particle_use.Et();
                                    //cout << "index = " << i_Particle_use << ", qp = " << qp << ", Et_total = " << Et_total << ", m2: " << m2 << endl;

                                    // Do momentum smearing and apply track efficiencies for PYTHIA
                                    if(eMode == 311) // PYTHIA
                                    {
                                        // Apply momentum smearing and track reconstruction efficiency
                                        PseudoJet Fill_PseudoJet_smear;
                                        Int_t track_acc = Apply_mom_smearing_and_efficiency(eflab_prim_glob,qp,eCentrality,i_Particle_use,m2,ePYTHIA_eff_factor,TLV_Particle_use,Fill_PseudoJet_smear,f_EfficiencyVsPt);

                                        if(track_acc)
                                        {
                                            particles[1].push_back(Fill_PseudoJet_smear);
                                        }
                                    }


                                    i_Particle_use++;
                                }
                            }
                        } // end of track loop


#if 0
                        //---------------------------------------------------------------
                        // A_ch analysis -> CMW effect (nothing to do with jets)
                        if(N_pos_neg_charges_reconstructed_Ach[0] + N_pos_neg_charges_reconstructed_Ach[1] > 0)
                        {
                            Double_t Ach  = ((Double_t)(N_pos_neg_charges_reconstructed_Ach[0] - N_pos_neg_charges_reconstructed_Ach[1]))/((Double_t)(N_pos_neg_charges_reconstructed_Ach[0] + N_pos_neg_charges_reconstructed_Ach[1]));
                            Double_t Mult = (Double_t)(N_pos_neg_charges_reconstructed_Ach[0] + N_pos_neg_charges_reconstructed_Ach[1]);

                            h_Ach           ->Fill(Ach);

                            for(Int_t i_pos_neg = 0; i_pos_neg < 2; i_pos_neg++) // loop first over all positive and then over all negative charges
                            {
                                for(Int_t i_tracks = 0; i_tracks < N_pos_neg_charges_reconstructed[i_pos_neg]; i_tracks++)
                                {
                                    p_pt_Ach[i_pos_neg] ->Fill(Ach,pt_tracks[i_pos_neg][i_tracks]);
                                }
                            }
                        }
                        //---------------------------------------------------------------
#endif


                        h_tracks_above_threshold_per_event->Fill(N_tracks_above_threshold);
                        for(Int_t i_hist = 0; i_hist < N2D_tracks_above_threshold; i_hist++)
                        {
                            h2D_tracks_above_threshold[i_hist] ->Fill(N_tracks_above_trigger,N_tracks_above_threshold_array[i_hist]);
                        }

                        Int_t Embed_user_index = i_Particle_use+1;
                        //---------------------------------------------------
                        // Single particle embedding for delta pt calculation
                        if(eMode == 42) // Delta pt calculation -> add one additional track
                        {
                            // Calculate pT bin for sampling distribution
                            Int_t epT_bin = 0;
                            for(Int_t i_pT = 0; i_pT < N_track_pt_bins_eta_phi; i_pT++)
                            {
                                if(pt_Delta_pt <= array_pt_bins_eta_phi[i_pT])
                                {
                                    epT_bin = i_pT;
                                    break;
                                }
                            }

                            // Sample eta and phi of the track based on same event distribution
                            Double_t track_phi;
                            Double_t track_eta;
                            h2D_track_eta_vs_phi[ez_bin][emult_bin][ePsi_bin][epT_bin]->GetRandom2(track_phi,track_eta);
                            //cout << "Embed pt = " << pt_Delta_pt << ", epT_bin = " << epT_bin << ", phi = " << track_phi << ", eta = " << track_eta << endl;
                            TLorentzVector TLV_Particle_Embed;
                            TLV_Particle_Embed.SetPtEtaPhiM(pt_Delta_pt,track_eta,track_phi,1.0);

                            PseudoJet Fill_PseudoJet(TLV_Particle_Embed.Px(),TLV_Particle_Embed.Py(),TLV_Particle_Embed.Pz(),TLV_Particle_Embed.E());
                            Fill_PseudoJet.set_user_index(Embed_user_index);
                            particles[0].push_back(Fill_PseudoJet);
                            particles_info.push_back(qp_Embed);  // taken from previous real track
                        }
                        //---------------------------------------------------


                        //---------------------------------------------------
                        Double_t PYTHIA_hard_bin_index_weight = 1.0;
                        if(eMode == 312)
                        {
                            // PYTHIA embedding for closure test

                            // YYY
                            // Get random hard bin
                            Int_t PYTHIA_hard_bin_index  = (Int_t)(h_PYTHIA_hard_bin_high_pT_N_events->GetRandom()); // Number of PYTHIA events with at least one 9 GeV/c track
                            if(eIn_Mode == 2)
                            {
                                while(PYTHIA_hard_bin_index < 5)
                                {
                                    PYTHIA_hard_bin_index  = (Int_t)(h_PYTHIA_hard_bin_high_pT_N_events->GetRandom()); // Number of PYTHIA events with at least one 9 GeV/c track
                                }
                            }
                            //Int_t PYTHIA_hard_bin_index  = ran_gen.Integer(11);

                            // Get PYTHIA hard bin weighting factor
                            PYTHIA_hard_bin_index_weight = h_PYTHIA_hard_bin_weigh_factors->GetBinContent(PYTHIA_hard_bin_index + 1);

                            //cout << "PYTHIA embedding, index: " << PYTHIA_hard_bin_index << ", PYTHIA_hard_bin_index_weight: " << PYTHIA_hard_bin_index_weight << endl;

                            if(PYTHIA_hard_bin_index >= 0 && PYTHIA_hard_bin_index < 11)
                            {
                                if(counter_PYTHIA[PYTHIA_hard_bin_index] >= N_PYTHIA_events[PYTHIA_hard_bin_index]) // start from the first PYTHIA event if all events were already used
                                {
                                    counter_PYTHIA[PYTHIA_hard_bin_index] = 0;
                                }
                                input_PYTHIA[PYTHIA_hard_bin_index]->GetEntry(counter_PYTHIA[PYTHIA_hard_bin_index]);

#if 0
                                cout << "" << endl;
                                cout << "PYTHIA embedding, index: " << PYTHIA_hard_bin_index << ", counter: " << counter_PYTHIA[PYTHIA_hard_bin_index]
                                    << ", PYTHIA_hard_bin_index_weight: " << PYTHIA_hard_bin_index_weight << endl;
#endif

                                Float_t  N_Particles_PYTHIA   = JetTrackEvent_PYTHIA[PYTHIA_hard_bin_index]->getNumParticle();

                                // Loop over the PYTHIA tracks
                                Int_t PYTHIA_user_index = 10000;
                                for(Int_t i_Particle_PYTHIA = 0; i_Particle_PYTHIA < N_Particles_PYTHIA; i_Particle_PYTHIA++)
                                {
                                    // Particle information
                                    JetTrackParticle_PYTHIA[PYTHIA_hard_bin_index] = JetTrackEvent_PYTHIA[PYTHIA_hard_bin_index]->getParticle(i_Particle_PYTHIA);
                                    Float_t dca_PYTHIA                 = JetTrackParticle_PYTHIA[PYTHIA_hard_bin_index]->get_dca_to_prim();
                                    Float_t m2_PYTHIA                  = JetTrackParticle_PYTHIA[PYTHIA_hard_bin_index]->get_Particle_m2 ();
                                    Float_t nSPi_PYTHIA                = JetTrackParticle_PYTHIA[PYTHIA_hard_bin_index]->get_Particle_nSigmaPi();
                                    Float_t nSK_PYTHIA                 = JetTrackParticle_PYTHIA[PYTHIA_hard_bin_index]->get_Particle_nSigmaK();
                                    Float_t nSP_PYTHIA                 = JetTrackParticle_PYTHIA[PYTHIA_hard_bin_index]->get_Particle_nSigmaP();
                                    Float_t qp_PYTHIA                  = JetTrackParticle_PYTHIA[PYTHIA_hard_bin_index]->get_Particle_qp();
                                    Float_t nhitsfit_PYTHIA            = JetTrackParticle_PYTHIA[PYTHIA_hard_bin_index]->get_Particle_hits_fit();
                                    TLorentzVector TLV_Particle_prim_PYTHIA = JetTrackParticle_PYTHIA[PYTHIA_hard_bin_index]->get_TLV_Particle_prim();
                                    TLorentzVector TLV_Particle_glob_PYTHIA = JetTrackParticle_PYTHIA[PYTHIA_hard_bin_index]->get_TLV_Particle_glob();

                                    TLorentzVector TLV_Particle_use_PYTHIA = TLV_Particle_prim_PYTHIA;
                                    if(eflab_prim_glob == 1) TLV_Particle_use_PYTHIA = TLV_Particle_glob_PYTHIA;

                                    Double_t track_pT_PYTHIA  = TLV_Particle_use_PYTHIA.Pt();
                                    if(track_pT_PYTHIA != track_pT_PYTHIA) continue; // that is a NaN test. It always fails if track_pT = nan.
                                    Double_t track_eta_PYTHIA = TLV_Particle_use_PYTHIA.PseudoRapidity();
                                    Double_t track_phi_PYTHIA = TLV_Particle_use_PYTHIA.Phi();

                                    TLorentzVector TLV_Particle_Embed;
                                    TLV_Particle_Embed.SetPtEtaPhiM(track_pT_PYTHIA,track_eta_PYTHIA,track_phi_PYTHIA,1.0);

                                    PseudoJet Fill_PseudoJet(TLV_Particle_Embed.Px(),TLV_Particle_Embed.Py(),TLV_Particle_Embed.Pz(),TLV_Particle_Embed.E());
                                    Fill_PseudoJet.set_user_index(PYTHIA_user_index);
                                    //particles[0].push_back(Fill_PseudoJet);
                                    particles_info.push_back(qp_PYTHIA);  // taken from previous real track


                                    // Apply momentum smearing and track reconstruction efficiency
                                    PseudoJet Fill_PseudoJet_smear_PYTHIA;

                                    Int_t track_acc = Apply_mom_smearing_and_efficiency(eflab_prim_glob,qp_PYTHIA,eCentrality,PYTHIA_user_index,m2_PYTHIA,ePYTHIA_eff_factor,TLV_Particle_use_PYTHIA,Fill_PseudoJet_smear_PYTHIA,f_EfficiencyVsPt);

                                    if(track_acc)
                                    {
                                        // Only difference for mode 312 is that [0] has no underlying heavy-ion event
                                        if(eIn_Mode == 2 && eJet_R > 0.45)
                                        {
                                            particles[0].push_back(Fill_PseudoJet); // only PYTHIA USE only for R = 0.5 closure test
                                        }
                                        else
                                        {
                                            particles[0].push_back(Fill_PseudoJet_smear_PYTHIA); // only PYTHIA
                                        }
                                        particles[1].push_back(Fill_PseudoJet_smear_PYTHIA); // PYTHIA + ME

                                        //----------------------------------------------
                                        // Trigger particles
                                        if(track_pT_PYTHIA > 0.2 && track_pT_PYTHIA < max_pt_threshold && fabs(track_eta_PYTHIA) < 1.0)
                                        {
                                            std::vector<Double_t> vec_in;
                                            vec_in.resize(4);
                                            vec_in[0] = track_eta_PYTHIA;
                                            vec_in[1] = track_phi_PYTHIA;
                                            vec_in[2] = track_pT_PYTHIA;
                                            vec_in[3] = PYTHIA_user_index;

                                            vec_trigger_tracks.push_back(vec_in);
                                        }
                                        //----------------------------------------------
                                    }

                                    //cout << "i_Particle_PYTHIA: " << i_Particle_PYTHIA << ", pt: " << track_pT_PYTHIA << endl;

                                    PYTHIA_user_index++;
                                    Embed_user_index++;
                                }

                                counter_PYTHIA[PYTHIA_hard_bin_index]++;
                            }
                            else
                            {
                                cout << "ERROR: PYTHIA_hard_bin_index out of range!" << endl;
                                continue;
                            }
                        }
                        //---------------------------------------------------



                        //-----------------------------------------------------------------------------


                        //-----------------------------------------------------------------------------
                        // choose a jet definition
                        JetDefinition jet_def(antikt_algorithm, jet_R);

                        // jet area definition
                        Double_t ghost_maxrap = 1.0; // Fiducial cut for background estimation
                        GhostedAreaSpec area_spec(ghost_maxrap);
                        //AreaDefinition area_def(active_area, area_spec);
                        AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));

                        // Loop only activates for PYTHIA
                        vector<PseudoJet> jets[2]; // [not smeared, smeared] only for PYTHIA
                        Double_t jet_rho_array[2];

                        std::vector< std::vector< std::vector<Float_t> > > vec_jet_array;
                        std::vector< std::vector< std::vector< std::vector<Float_t> > > > vec_jet_const_array; // jet constituent information

                        vec_jet_array.resize(2);       // [not smeared, smeared] only for PYTHIA
                        vec_jet_const_array.resize(2); // [not smeared, smeared] only for PYTHIA

#if 0
                        cout << "tracks original PYTHIA: " << (Int_t)particles[0].size() << ", tracks PYTHIA (smeared+eff): " << (Int_t)particles[1].size() << endl;
#endif

                        std::vector< std::vector<Int_t> > random_trigger_index_array;
                        random_trigger_index_array.resize(2); // i_orig_smear
                        std::vector<Int_t>    vec_flag_low_pt_PYTHIA;
                        std::vector<Double_t> vec_flag_low_pt_PYTHIA_jet_pt;

                        for(Int_t i_orig_smear = 0; i_orig_smear < N_orig_smear; i_orig_smear++)
                        {
                            random_trigger_index_array[i_orig_smear].resize(N_leading_pt_bins);
                            for(Int_t lead_pt_bin = 0; lead_pt_bin < N_leading_pt_bins; lead_pt_bin++)
                            {
                                random_trigger_index_array[i_orig_smear][lead_pt_bin] = -1;
                            }

                            //AreaDefinition area_def(active_area, GhostedAreaSpec(ghost_maxrap));
                            //ClusterSequenceArea clust_seq(particles, jet_def, area_def);
                            ClusterSequenceArea clust_seq_hard(particles[i_orig_smear], jet_def, area_def);

                            // run the clustering, extract the jets
                            //ClusterSequence clust_seq(particles, jet_def);
                            double ptmin = 0.2;
                            vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));
                            Selector Fiducial_cut_selector = SelectorAbsEtaMax(1.0 - jet_R); // Fiducial cut for jets
                            jets[i_orig_smear] = Fiducial_cut_selector(jets_all);
                            //vector<PseudoJet> jets[i_orig_smear] = sorted_by_pt(clust_seq.inclusive_jets());

                            // print out some info
                            //cout << "Clustered with " << jet_def.description() << endl;



                            // background estimation
                            //cout << "Define JetDefinition" << endl;
                            //JetDefinition jet_def_bkgd(kt_algorithm, 0.4);
                            JetDefinition jet_def_bkgd(kt_algorithm, jet_R_background); // <--
                            //JetDefinition jet_def_bkgd(antikt_algorithm, jet_R); // test
                            //cout << "Define AreaDefinition" << endl;
                            AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
                            //AreaDefinition area_def_bkgd(active_area,GhostedAreaSpec(ghost_maxrap,1,0.005));
                            //cout << "Define selector" << endl;
                            //Selector selector = SelectorAbsRapMax(1.0) * (!SelectorNHardest(Remove_N_hardest)); // 2
                            Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(Remove_N_hardest)); // <--
                            //Selector selector = SelectorAbsEtaMax(1.0 - jet_R); // test


                            //cout << "Define JetMedianBackgroundEstimator" << endl;
                            JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd); // <--
                            //JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def, area_def); // test
                            //cout << "Define Subtractor" << endl;
                            Subtractor subtractor(&bkgd_estimator);
                            //cout << "Define bkgd_estimator" << endl;
                            bkgd_estimator.set_particles(particles[i_orig_smear]);

                            //cout << "Calculate jet_rho and jet_sigma" << endl;
                            //Double_t jet_rho   = bkgd_estimator.rho() - 0.212;   // TO BE CHANGED, equivallent to a 60 MeV shift for R = 0.3 (rho*pi*R*2 = 60 MeV)
                            Double_t jet_rho   = bkgd_estimator.rho();
                            jet_rho_array[i_orig_smear] = jet_rho;
                            //cout << "jet_sigma" << endl;
                            Double_t jet_sigma = bkgd_estimator.sigma();

                            //cout << "jet_rho = " << jet_rho << ", jet_sigma = " << jet_sigma << endl;


                            //h_jet_per_event[0] ->Fill((Int_t)jets.size());

                            //----------------------------------------------------------------------------------------
                            // Calculate Et weight (SE and ME are different due to statistical average effect)
                            Float_t SE_ME_Et_val[2];
                            Float_t Et_weight_factor = 1.0;
                            if(eRandom != -1)
                            {
                                for(Int_t iSE_ME_Et = 0; iSE_ME_Et < 2; iSE_ME_Et++)
                                {
                                    //SE_ME_Et_val[iSE_ME_Et] = h_jet_Et_array_In[ez_bin][emult_bin][ePsi_bin][iSE_ME_Et] ->GetBinContent(h_jet_Et_array_In[ez_bin][emult_bin][ePsi_bin][iSE_ME_Et]->FindBin(Et_total));
                                    //SE_ME_Et_val[iSE_ME_Et] = h_tracks_per_event_array_In[ez_bin][emult_bin][ePsi_bin][iSE_ME_Et]  ->GetBinContent(h_tracks_per_event_array_In[ez_bin][emult_bin][ePsi_bin][iSE_ME_Et] ->FindBin(N_Particles));
                                    SE_ME_Et_val[iSE_ME_Et] = h2D_jet_rho_vs_mult_array_In[ez_bin][emult_bin][ePsi_bin][iSE_ME_Et] ->GetBinContent(h2D_jet_rho_vs_mult_array_In[ez_bin][emult_bin][ePsi_bin][iSE_ME_Et]->FindBin(N_Particles,jet_rho));
                                }
                                //cout << "SE_ME_Et_val[0] = " << SE_ME_Et_val[0] << ", SE_ME_Et_val[1] = " << SE_ME_Et_val[1] << endl;

                                if(SE_ME_Et_val[0] > 0.0 && SE_ME_Et_val[1] > 0.0)
                                {
                                    Et_weight_factor = SE_ME_Et_val[0]/SE_ME_Et_val[1];
                                }
                            }
                            if(!(eIn_Mode == 2 || eIn_Mode == 24)) Et_weight_factor = 1.0; // only relevant for mixed event

                            //Et_weight_factor = 1.0;
                            Et_weight_factor *= reweight;
                            //cout << "Et_weight_factor = " << Et_weight_factor << endl;
                            //----------------------------------------------------------------------------------------

                            h_jet_rho[0][i_orig_smear]       ->Fill(jet_rho);
                            if(i_orig_smear == 0)
                            {
                                h2D_jet_rho_vs_mult_array[ez_bin][emult_bin][ePsi_bin] ->Fill(N_Particles,jet_rho,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                h_jet_rho_array[ez_bin][emult_bin][ePsi_bin]           ->Fill(jet_rho,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                h_jet_Et_array[ez_bin][emult_bin][ePsi_bin]            ->Fill(Et_total,PYTHIA_hard_bin_index_weight);
                                h_jet_Et_array_weight[ez_bin][emult_bin][ePsi_bin]     ->Fill(Et_total,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                h_jet_rho_vs_Et                                        ->Fill(Et_total,jet_rho,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                h_jet_per_event_array[ez_bin][emult_bin][ePsi_bin]     ->Fill((Int_t)jets[i_orig_smear].size(),Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                h_tracks_per_event_array[ez_bin][emult_bin][ePsi_bin]  ->Fill(N_Particles,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                            }

                            //cout << "jet_rho = " << jet_rho << ", N_Particles = " << N_Particles << ", i_orig_smear: " << i_orig_smear << endl;

                            // get the subtracted jets
                            Int_t jets_per_event[N_jet_areas];
                            Int_t dijets_per_event[N_jet_areas];
                            for(Int_t i_area = 0; i_area < N_jet_areas; i_area++)
                            {
                                jets_per_event[i_area]   = 0;
                                dijets_per_event[i_area] = 0;
                            }

                            if(jet_rho >= 0.0)
                            {
                                //cout << "Accepted" << endl;
                                Int_t flag_save_to_ntuple  = 0;
                                Int_t jets_above_threshold = 0;
                                for(Int_t i = 0; i < (Int_t)jets[i_orig_smear].size(); i++)
                                {
                                    Float_t jet_pt     = jets[i_orig_smear][i].perp();
                                    Float_t jet_area   = jets[i_orig_smear][i].area();
                                    Float_t jet_pt_sub = jets[i_orig_smear][i].perp() - jet_rho*jet_area;
                                    Float_t jet_eta    = jets[i_orig_smear][i].eta();
                                    Float_t jet_phi    = jets[i_orig_smear][i].phi();

                                    //cout << "i: " << i << ", jet_pt_sub: " << jet_pt_sub << ", jet_pt: " << jet_pt << endl;
                                    if(jet_pt_sub > 9.0) jets_above_threshold++;

                                    //if(jet_pt_sub > 2500.0) flag_save_to_ntuple = 1; // at least one jet candidate with high pt in event
                                    //flag_save_to_ntuple = 1;
                                }

                                if(jets_above_threshold > 1)
                                {
                                    Int_t N_good_triggers_ntuple = 0;
                                    for(Int_t i_val = 0; i_val < (Int_t)vec_trigger_tracks.size(); i_val++)
                                    {
                                        if(vec_trigger_tracks[i_val][2] >= 9.0) N_good_triggers_ntuple++;
                                    }
                                    if(N_good_triggers_ntuple > 0)
                                    {
                                        flag_save_to_ntuple = 1;
                                        //cout << "Jet candidate with " << jets_above_threshold << " jet(s) above threshold found in event and " <<
                                        //    N_good_triggers_ntuple <<  " triggers, " << counter << " -> will be stored in Ntuple" << endl;
                                    }
                                }
                                flag_save_to_ntuple = 0;

                                //h_jet_per_event[1] ->Fill((Int_t)jets[i_orig_smear].size());

#if 0
                                if(i_orig_smear == 0)
                                {
                                    cout << "RunId: " << RunId << ", rho: " << jet_rho << ", N total particles: " << N_Particles << endl;
                                    for(Int_t i_track = 0; i_track < (Int_t)particles[0].size(); i_track++)
                                    {
                                        cout << "i_track: " << i_track << ", pT: " << particles[0][i_track].perp() << ", eta: " << particles[0][i_track].pseudorapidity() << endl;
                                    }
                                }
#endif

                                //----------------------------------------------------------------------------------------
                                // Jet A loop
                                vec_jet_array[i_orig_smear].resize((Int_t)jets[i_orig_smear].size());
                                vec_jet_const_array[i_orig_smear].resize((Int_t)jets[i_orig_smear].size());


                                //-------------------
                                // for mode 312, save user index of tracks for jets with very low pT
                                // those jets will be downscaled for the closure test since they
                                // are not similar to what is seen in data.
                                if(eMode == 312 && i_orig_smear == 0)
                                {
                                    for(Int_t i = 0; i < (Int_t)jets[i_orig_smear].size(); i++)
                                    {
                                        Float_t jet_pt     = jets[i_orig_smear][i].perp();
                                        vector<PseudoJet> jet_constituents = jets[i_orig_smear][i].constituents();

                                        if(jet_pt < jet_pt_cut_closure)
                                        {
                                            for(Int_t j = 0; j < (Int_t)jet_constituents.size(); j++)
                                            {
                                                Float_t jet_const_pt  = jet_constituents[j].perp();
                                                Float_t jet_const_phi = jet_constituents[j].phi();
                                                Float_t jet_const_eta = jet_constituents[j].eta();
                                                Int_t   user_index    = jet_constituents[j].user_index();
                                                //cout << "jet: " << i << ", track: " << j << ", jet_pt: " << jet_pt
                                                //    << ", track pt: " << jet_const_pt << ", user_index: " << user_index << endl;
                                                if(user_index >= 10000) // PYTHIA embedded track
                                                {
                                                    vec_flag_low_pt_PYTHIA.push_back(user_index);
                                                    vec_flag_low_pt_PYTHIA_jet_pt.push_back(jet_pt);
                                                }
                                            }
                                            //cout << "START, jet_pt: " << jet_pt << endl;
                                            //for(Int_t i_flagged = 0; i_flagged < (Int_t)vec_flag_low_pt_PYTHIA.size(); i_flagged++)
                                            //{
                                            //    cout << "i_flagged: " << i_flagged << ", user_index: " << vec_flag_low_pt_PYTHIA[i_flagged] << endl;
                                            //}

                                        }
                                    }
                                }
                                //-------------------


                                for(Int_t i = 0; i < (Int_t)jets[i_orig_smear].size(); i++)
                                {
                                    Float_t jet_pt     = jets[i_orig_smear][i].perp();
                                    Float_t jet_area   = jets[i_orig_smear][i].area();
                                    Float_t jet_pt_sub = jets[i_orig_smear][i].perp() - jet_rho*jet_area;
                                    Float_t jet_eta    = jets[i_orig_smear][i].eta();
                                    Float_t jet_phi    = jets[i_orig_smear][i].phi();

                                    vec_jet_array[i_orig_smear][i].push_back(jet_pt); // save the jet information in the array for later PYTHIA matching studies
                                    vec_jet_array[i_orig_smear][i].push_back(jet_area); // save the jet information in the array for later PYTHIA matching studies
                                    vec_jet_array[i_orig_smear][i].push_back(jet_pt_sub); // save the jet information in the array for later PYTHIA matching studies
                                    vec_jet_array[i_orig_smear][i].push_back(jet_eta); // save the jet information in the array for later PYTHIA matching studies
                                    vec_jet_array[i_orig_smear][i].push_back(jet_phi); // save the jet information in the array for later PYTHIA matching studies

                                    if((i_orig_smear == 0 && eMode != 312) || (eMode == 312 && i_orig_smear == 1))
                                    {
                                        h_area_vs_jet_pt ->Fill(jet_pt_sub,jet_area,PYTHIA_hard_bin_index_weight);
                                        h_area           ->Fill(jet_area,PYTHIA_hard_bin_index_weight);
                                    }

                                    //cout << "Event: " << counter
                                    //    << ", i_orig_smear: " << i_orig_smear << ", jet: " << i << ", jet_area: " << jet_area << ", jet_rho: " << jet_rho << endl;

                                    vector<PseudoJet> jet_constituents = jets[i_orig_smear][i].constituents();
                                    //cout << "jet " << i << ", jet_pt = " << jet_pt << ", number of constituents = " << (Int_t)jet_constituents.size() << endl;
                                    Float_t leading_pt      = 0.0;
                                    Float_t sub_leading_pt  = 0.0;
                                    Float_t leading_phi     = 0.0;
                                    Float_t leading_eta     = 0.0;
                                    Int_t   N_tracks_above_pT_assoc_threshold = 0;


                                    // Determine leading and sub_leading pt values
                                    vec_jet_const_array[i_orig_smear][i].resize((Int_t)jet_constituents.size());
                                    for(Int_t j = 0; j < (Int_t)jet_constituents.size(); j++)
                                    {
                                        Float_t jet_const_pt  = jet_constituents[j].perp();
                                        Float_t jet_const_phi = jet_constituents[j].phi();
                                        Float_t jet_const_eta = jet_constituents[j].eta();
                                        Int_t   user_index    = jet_constituents[j].user_index();

                                        if(jet_const_pt > 0.05)
                                        {
                                            //cout << "jet: " << j << ", user_index: " << user_index << ", jet_const_pt: " << jet_const_pt
                                            //    << ", m2: " << particle_info[0][user_index] << endl;
                                        }

                                        vec_jet_const_array[i_orig_smear][i][j].push_back(jet_const_pt);
                                        vec_jet_const_array[i_orig_smear][i][j].push_back(jet_const_phi);
                                        vec_jet_const_array[i_orig_smear][i][j].push_back(jet_const_eta);
                                        vec_jet_const_array[i_orig_smear][i][j].push_back((Float_t)user_index);

                                        Double_t Psi2_use = EP_eta_pos_corr;
                                        if(jet_const_eta > 0.0) Psi2_use = EP_eta_neg_corr;

                                        Double_t phi_use  = jet_const_phi;

                                        // -pi/2..pi/2 -> 0..pi
                                        if(Psi2_use < 0.0) Psi2_use += TMath::Pi();

                                        // -pi..pi -> 0..pi
                                        if(phi_use < 0.0) phi_use = phi_use + TMath::Pi();

                                        // -pi..pi, delta_phi_angle: 0..pi
                                        Double_t delta_phi_angle = phi_use - Psi2_use;
                                        if(phi_use >= 0.0 && delta_phi_angle >= 0.0) delta_phi_angle = delta_phi_angle;
                                        if(phi_use >= 0.0 && delta_phi_angle < 0.0)  delta_phi_angle += TMath::Pi();
                                        if(phi_use < 0.0  && delta_phi_angle >= -TMath::Pi()) delta_phi_angle += TMath::Pi();
                                        if(phi_use < 0.0  && delta_phi_angle < -TMath::Pi())  delta_phi_angle += 2.0*TMath::Pi();

                                        Double_t v2_val = TMath::Cos(2.0*delta_phi_angle);

                                        //--------------------------------------------------------
                                        if(eMode == 42)
                                        {
                                            if(i_orig_smear == 0 && jet_const_pt > 0.15)
                                            {
                                                if(i == 0 && j == 0)
                                                {
                                                    h_Psi_Full   ->Fill(Psi2_use);
                                                    h_Psi_etapos ->Fill(EP_eta_pos_corr);
                                                    h_Psi_etaneg ->Fill(EP_eta_neg_corr);
                                                }
                                                p_v2_vs_pt_jet ->Fill(jet_const_pt,v2_val);
                                                h_phi          ->Fill(jet_const_phi);
                                            }
                                            // check if embedded particle is in the jet
                                            if(user_index == Embed_user_index)
                                            {
                                                //cout << "Embedded particle: pt = " << jet_const_pt << ", phi = " << jet_const_phi << ", eta = " << jet_const_eta << endl;
                                                Double_t Delta_pt = jet_pt_sub - jet_const_pt;

                                                for(Int_t i_area = 0; i_area < N_jet_areas; i_area++)
                                                {
                                                    if(jet_area > array_areas[i_area])
                                                    {
                                                        // Define weighting based on elliptic flow
                                                        Double_t v2_max     = 0.1; // maximum v2 for 200 GeV 0-10%
                                                        if(eCentrality == 1) v2_max = 0.25;
                                                        Double_t EP_res     = 0.465; // ~eta sub EP resolution for 0-10% and 60-80%
                                                        Double_t v2_high_pt = 0.04; // estimate of v2 at high pT
                                                        if(eCentrality == 1) v2_high_pt = 0.08;
                                                        Double_t v2_use = v2_max*(1.0/EP_res)*jet_const_pt/3.0; // linear increase of v2 with pT up to 3.0 GeV
                                                        if(jet_const_pt > 3.0 && jet_const_pt <= 7.0) // dropping v2 from max value to v2 at high pT value
                                                        {
                                                            Double_t slope_val = (v2_high_pt*(1.0/EP_res) - v2_max*(1.0/EP_res))/(7.0 - 3.0);
                                                            Double_t t_val     = v2_max*(1.0/EP_res) - slope_val * 3.0;
                                                            v2_use             = slope_val*jet_const_pt + t_val;
                                                        }
                                                        if(jet_const_pt > 7.0) // constant v2 at high pT
                                                        {
                                                            v2_use = v2_high_pt*(1.0/EP_res);
                                                        }
                                                        Double_t dNdDeltaphi = 1.0 + 2.0*0.1*TMath::Cos(2.0*jet_const_phi-2.0*Psi2_use); // weighting factor for embedded tracks with assumed v2 of 0.05 (EP resolution ~0.5 -> v2 factor here is 0.1)

                                                        //if(i_area == 0) cout << "pt: " << jet_const_pt << ", v2_use: " << v2_use << ", weight: " << dNdDeltaphi << endl;

                                                        h_Delta_pt_vs_embed_pt[i_area][i_orig_smear][ePsi_bin]        ->Fill(jet_const_pt,Delta_pt,PYTHIA_hard_bin_index_weight);
                                                        h_Delta_pt_vs_embed_pt_weight[i_area][i_orig_smear][ePsi_bin] ->Fill(jet_const_pt,Delta_pt,PYTHIA_hard_bin_index_weight*dNdDeltaphi);
                                                        //cout << "i_area = " << i_area << ", epT_bin = " << epT_bin << ", i_orig_smear = " << i_orig_smear << ", ePsi_bin = " << ePsi_bin << ", Delta_pt = " << Delta_pt << endl;
                                                    }
                                                }
                                            }
                                        }
                                        //--------------------------------------------------------


                                        //cout << "jet counter = " << i << ", track counter = " << j << ", track pt = " << jet_const_pt << endl;

                                        Float_t qp_val = 1.0;
                                        Int_t part_info_size = (Int_t)particles_info.size();
                                        if(user_index >= 0 && user_index < part_info_size)
                                        {
                                            qp_val = particles_info[user_index];
                                            if(qp_val > 0.0)
                                            {
                                                qp_val = 1.0;
                                            }
                                            else
                                            {
                                                qp_val = -1.0;
                                            }

                                        }

                                        //cout << "user_index = " << user_index << ", qp_val = " << qp_val << ", part_info_size = " << part_info_size << endl;

                                        if(flag_save_to_ntuple && i_orig_smear == 0)
                                        {
                                            // EventId:JetId:rho:area:Jetphi:Jeteta:Jetpt:TrackId:eta:phi:pt
                                            ReCoil_Jet_NTDataArray[0]  = (Float_t)counter;
                                            ReCoil_Jet_NTDataArray[1]  = (Float_t)i;
                                            ReCoil_Jet_NTDataArray[2]  = (Float_t)jet_rho;
                                            ReCoil_Jet_NTDataArray[3]  = (Float_t)jet_area;
                                            ReCoil_Jet_NTDataArray[4]  = (Float_t)jet_phi;
                                            ReCoil_Jet_NTDataArray[5]  = (Float_t)jet_eta;
                                            ReCoil_Jet_NTDataArray[6]  = (Float_t)jet_pt_sub;
                                            ReCoil_Jet_NTDataArray[7]  = (Float_t)j;
                                            ReCoil_Jet_NTDataArray[8]  = (Float_t)jet_const_eta;
                                            ReCoil_Jet_NTDataArray[9]  = (Float_t)jet_const_phi;
                                            ReCoil_Jet_NTDataArray[10] = (Float_t)qp_val*jet_const_pt;
                                            ReCoil_Jet_NTDataArray[11] = (Float_t)prim_vertex_x;
                                            ReCoil_Jet_NTDataArray[12] = (Float_t)prim_vertex_y;
                                            ReCoil_Jet_NTDataArray[13] = (Float_t)prim_vertex_z;

                                            NT_ReCoil_Jet->Fill(ReCoil_Jet_NTDataArray);
                                        }

                                        if(jet_const_pt > leading_pt)
                                        {
                                            sub_leading_pt = leading_pt;
                                            leading_pt     = jet_const_pt;
                                            leading_phi    = jet_const_phi;
                                            leading_eta    = jet_const_eta;
                                        }
                                        //cout << "track " << j << ", track pt = " << jet_const_pt << ", phi = " << jet_const_phi << ", eta = " << jet_const_eta << endl;
                                        if(jet_const_pt > track_pT_assoc_threshold)
                                        {
                                            N_tracks_above_pT_assoc_threshold++;
                                        }
                                    }

                                    if(eMode == 32 || eMode == 312) // for eMode == 32 the high pT tracks are downscaled. Make sure all trigger pT bins are filled.
                                    {
                                        leading_pt *= ran_gen.Rndm()*max_pt_threshold; // Not used anymore
                                    }

                                    //cout << "leading_pt = " << leading_pt << ", sub_leading_pt = " << sub_leading_pt << endl;
                                    if(sub_leading_pt > 0.0 && leading_pt > 0.0 && i_orig_smear == 0)
                                    {
                                        h_ratio_sub_lead_to_lead_pt_vs_lead_pt ->Fill(leading_pt,sub_leading_pt/leading_pt,PYTHIA_hard_bin_index_weight);
                                        h_sub_lead_vs_lead_pt ->Fill(leading_pt,sub_leading_pt,PYTHIA_hard_bin_index_weight);
                                    }


                                    //----------------------------------------------------------------
                                    if(eMode == 31 || eMode == 310 || eMode == 32 || eMode == 312 || eMode == 311 || eMode == 3)
                                    {
                                        if(i == 0)
                                        {
#if 0
                                            for(Int_t i_val = 0; i_val < (Int_t)vec_trigger_tracks.size(); i_val++) cout << "i before: " << i_val << ", eta: " << vec_trigger_tracks[i_val][0] << ", phi: " << vec_trigger_tracks[i_val][1] << ", pT: " << vec_trigger_tracks[i_val][2] << endl;
#endif
                                            // sort the pT values
                                            std::sort (vec_trigger_tracks.begin(), vec_trigger_tracks.end(), sortFunc);

#if 0
                                            for(Int_t i_val = 0; i_val < (Int_t)vec_trigger_tracks.size(); i_val++) cout << "i after: " << i_val << ", eta: " << vec_trigger_tracks[i_val][0] << ", phi: " << vec_trigger_tracks[i_val][1] << ", pT: " << vec_trigger_tracks[i_val][2] << endl;
#endif

                                            // loop over all defined trigger pt intervals
                                            for(Int_t lead_pt_bin = 0; lead_pt_bin < N_leading_pt_bins; lead_pt_bin++)
                                            {
                                                if(lead_pt_bin >= 0)
                                                {
                                                    h_N_tracks_dijet[lead_pt_bin][i_orig_smear]           ->Fill(N_Particles,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                    h_N_tracks_dijet[N_leading_pt_bins][i_orig_smear]     ->Fill(N_Particles,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                }
                                                if(global_bin >= 0 && global_bin < N_global_bin && lead_pt_bin >= 0)
                                                {
                                                    h2D_mult_vs_global_bin[lead_pt_bin][i_orig_smear]       ->Fill(global_bin,N_Particles,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                    h2D_mult_vs_global_bin[N_leading_pt_bins][i_orig_smear] ->Fill(global_bin,N_Particles,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                }

                                                // Find first index of sorted pT vector which fits into trigger pT range
                                                Int_t trigger_index_min = -1;
                                                for(Int_t i_trigger = 0; i_trigger < (Int_t)vec_trigger_tracks.size(); i_trigger++)
                                                {
                                                    if(vec_trigger_tracks[i_trigger][2] >= Array_leading_pt_bins[0][lead_pt_bin])
                                                    {
                                                        trigger_index_min = i_trigger;
                                                        break;
                                                    }
                                                }

                                                //cout << "lead_pt_bin: " << lead_pt_bin << ", random_trigger_index: " << random_trigger_index << endl;
                                                if(trigger_index_min == -1 && !(eMode == 32))
                                                {
                                                    random_trigger_index = -1;
                                                    //cout << "Set index to -1" << endl;
                                                    continue; // Skip event if it doesn't fit into trigger range
                                                }

                                                // Choose a random track within trigger pT range as trigger
                                                //Int_t random_trigger_index;

                                                Int_t trigger_index_range = (Int_t)vec_trigger_tracks.size()-trigger_index_min;

                                                if(eMode == 32) // for eMode == 32 the high pT tracks are downscaled. Make sure all trigger pT bins are filled.
                                                {
                                                    random_trigger_index = ran_gen.Integer((Int_t)vec_trigger_tracks.size());
                                                    trigger_index_range  = 1;
                                                    trigger_index_min    = 0;
                                                }
                                                else
                                                {
                                                    random_trigger_index = ran_gen.Integer(trigger_index_range); // [0,trigger_index_range-1]
                                                    random_trigger_index += trigger_index_min;
                                                }

#if 0
                                                cout << "jetA: " << i << ", lead_pt_bin: " << lead_pt_bin << ", low val: " << Array_leading_pt_bins[0][lead_pt_bin]
                                                    << ", trigger_index_min: " <<  trigger_index_min
                                                    << ", pt val: " << vec_trigger_tracks[trigger_index_min][2]
                                                    << ", trigger index range: " << (Int_t)vec_trigger_tracks.size()-trigger_index_min
                                                    << ", random_trigger_index: " << random_trigger_index
                                                    << endl;

#endif

                                                if(random_trigger_index >= (Int_t)vec_trigger_tracks.size() || random_trigger_index < 0)
                                                {
                                                    cout << "ERROR: random_trigger_index out of range!, random_trigger_index: " << random_trigger_index << endl;
                                                    continue;
                                                }

                                                random_trigger_index_array[i_orig_smear][lead_pt_bin] = random_trigger_index;
                                                if(i_orig_smear == 1) random_trigger_index = random_trigger_index_array[0][lead_pt_bin];

                                                if(random_trigger_index < 0) // for i_orig_smear == 1, can happen since its copied from random_trigger_index_array
                                                {
                                                    continue;
                                                }


#if 1
                                                //------------------------------------------------------------------------------------------------------------
                                                // Determine angle between trigger and its parent jet
                                                for(Int_t i_jet_trig = 0; i_jet_trig < (Int_t)jets[i_orig_smear].size(); i_jet_trig++)
                                                {
                                                    Float_t jet_pt     = jets[i_orig_smear][i_jet_trig].perp();
                                                    Float_t jet_area   = jets[i_orig_smear][i_jet_trig].area();
                                                    Float_t jet_eta    = jets[i_orig_smear][i_jet_trig].eta();
                                                    Float_t jet_phi    = jets[i_orig_smear][i_jet_trig].phi();

                                                    vector<PseudoJet> jet_constituents_trig = jets[i_orig_smear][i_jet_trig].constituents();
                                                    for(Int_t j_trig = 0; j_trig < (Int_t)jet_constituents_trig.size(); j_trig++)
                                                    {
                                                        Float_t jet_const_pt  = jet_constituents_trig[j_trig].perp();
                                                        Float_t jet_const_phi = jet_constituents_trig[j_trig].phi();
                                                        Float_t jet_const_eta = jet_constituents_trig[j_trig].eta();
                                                        Int_t   user_index    = jet_constituents_trig[j_trig].user_index();

                                                        //vec_trigger_tracks[random_trigger_index][0]; // eta
                                                        Double_t trigger_phi     = vec_trigger_tracks[random_trigger_index][1]; // phi
                                                        Int_t trigger_user_index = (Int_t)vec_trigger_tracks[random_trigger_index][3]; // user index

                                                        if(user_index == trigger_user_index && i_orig_smear == 0 && lead_pt_bin == 9)
                                                        {
                                                            Double_t delta_phi_trigger_jet = fabs(trigger_phi - jet_phi);
                                                            if(delta_phi_trigger_jet > 2.0*Pi)  delta_phi_trigger_jet -= 2.0*Pi;
                                                            h_Delta_phi_trig_jet ->Fill(delta_phi_trigger_jet);
                                                            //cout << "i_orig_smear: " << i_orig_smear << ", i_jet_trig: " << i_jet_trig << ", j_trig: " << j_trig
                                                            //    << ", lead_pt_bin: " << lead_pt_bin <<
                                                            //    ", user_index: " << user_index << ", trigger_user_index: " << trigger_user_index << ", phi_track: " << jet_const_phi <<
                                                            //    ", phi_trigger: " << trigger_phi << ", delta_phi_trigger_jet: " << delta_phi_trigger_jet << endl;
                                                            break;
                                                        }
                                                    }
                                                }
                                                //------------------------------------------------------------------------------------------------------------
#endif


                                                // Loop over ALL possible triggers in event, for i_trigger_use == 0 use the random trigger defined above
                                                for(Int_t i_trigger_use = 0; i_trigger_use < trigger_index_range+1; i_trigger_use++)
                                                {
                                                    //cout << "i_trigger_use: " << i_trigger_use << endl;
                                                    if(i_trigger_use > 0 && lead_pt_bin < 4) break; // Do not fill the histograms for the low pT bins -> not needed, takes too much time
                                                    Int_t trigger_orig_smear = 0;
                                                    if(i_trigger_use > 0) // not the random trigger, loop overa ALL triggers in event
                                                    {
                                                        trigger_orig_smear = 2; // That is for all triggers, index = 0 and index 1 (i_orig_smear) is for random single trigger
                                                        random_trigger_index = trigger_index_min + i_trigger_use - 1; // loop over all good triggers
                                                    }


                                                    h_trigger_track[i_orig_smear+trigger_orig_smear] ->Fill(N_leading_pt_bins,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                    h_trigger_track[i_orig_smear+trigger_orig_smear] ->Fill(lead_pt_bin,Et_weight_factor*PYTHIA_hard_bin_index_weight);

                                                    h_trigger_track_vs_global_bin[i_orig_smear+trigger_orig_smear]->Fill(global_bin,N_leading_pt_bins,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                    h_trigger_track_vs_global_bin[i_orig_smear+trigger_orig_smear]->Fill(global_bin,lead_pt_bin,Et_weight_factor*PYTHIA_hard_bin_index_weight);



                                                    // Jet B loop
                                                    Int_t N_accepted_jets[N_jet_areas];
                                                    for(Int_t i_area_acc = 0; i_area_acc < N_jet_areas; i_area_acc++)
                                                    {
                                                        N_accepted_jets[i_area_acc] = 0;
                                                    }

                                                    for(Int_t iB = 0; iB < (Int_t)jets[i_orig_smear].size(); iB++)
                                                    {
                                                        //if(iB == i) continue; // Don't use the same jet candidate twice
                                                        Float_t jet_ptB     = jets[i_orig_smear][iB].perp();
                                                        Float_t jet_areaB   = jets[i_orig_smear][iB].area();
                                                        Float_t jet_pt_subB = jets[i_orig_smear][iB].perp() - jet_rho*jet_areaB;
                                                        Float_t jet_etaB    = jets[i_orig_smear][iB].eta();
                                                        Float_t jet_phiB    = jets[i_orig_smear][iB].phi();


                                                        vector<PseudoJet> jet_constituentsB = jets[i_orig_smear][iB].constituents();


                                                        //------------------
                                                        // Remove very low pT jets (put them close to 0). Unphysical.
                                                        if(eMode == 312)
                                                        {
                                                            Double_t jet_pt_original = 0.0;
                                                            Double_t weight_closure  = 1.0;
                                                            Int_t flag_low_pt = 0;
                                                            for(Int_t j = 0; j < (Int_t)jet_constituentsB.size(); j++)
                                                            {
                                                                Float_t jet_const_pt  = jet_constituentsB[j].perp();
                                                                Float_t jet_const_phi = jet_constituentsB[j].phi();
                                                                Float_t jet_const_eta = jet_constituentsB[j].eta();
                                                                Int_t   user_index    = jet_constituentsB[j].user_index();
                                                                if(user_index >= 10000) // PYTHIA embedded track
                                                                {
                                                                    for(Int_t i_flagged = 0; i_flagged < (Int_t)vec_flag_low_pt_PYTHIA.size(); i_flagged++)
                                                                    {
                                                                        if(user_index == vec_flag_low_pt_PYTHIA[i_flagged])
                                                                        {
                                                                            jet_pt_original = vec_flag_low_pt_PYTHIA_jet_pt[i_flagged];
                                                                            weight_closure = 0.65*TMath::Sqrt(jet_pt_original)/TMath::Sqrt(jet_pt_cut_closure) + 0.05;
                                                                            flag_low_pt = 1;

                                                                            if(ran_gen.Rndm() > weight_closure)
                                                                            {
                                                                                jet_ptB     -= jet_const_pt;
                                                                                jet_pt_subB -= jet_const_pt;
                                                                                //if(i_orig_smear == 0) cout << "jet_ptB: " << jet_ptB << endl;
                                                                            }

                                                                            //cout << "jetA: " << i << ", jetB: " << iB << ", lead_pt_bin: " << lead_pt_bin << ", track: " << j << ", jet_pt: " << jet_ptB
                                                                            //    << ", original jet pt: " << jet_pt_original << ", track pt: " << jet_const_pt << ", user_index: " << user_index << endl;

                                                                            break;
                                                                        }
                                                                    }
                                                                    //if(flag_low_pt)
                                                                    //{
                                                                    //    break;
                                                                    //}
                                                                }
                                                            }
                                                            //if(flag_low_pt) cout << "jet_ptB: " << jet_ptB << ", jet_pt_original: " << jet_pt_original << ", weight_closure: " << weight_closure << endl;
                                                        }
#if 0
                                                        if(jet_ptB < 3.0 && i_orig_smear == 0 && trigger_orig_smear == 0)
                                                        {
                                                            cout << "jet_ptB: " << jet_ptB << ", size: " << (Int_t)vec_flag_low_pt_PYTHIA.size() << endl;
                                                            for(Int_t i_flagged = 0; i_flagged < (Int_t)vec_flag_low_pt_PYTHIA.size(); i_flagged++)
                                                            {
                                                                cout << "i_flagged: " << i_flagged << ", user_index: " << vec_flag_low_pt_PYTHIA[i_flagged] << endl;
                                                            }
                                                            for(Int_t j = 0; j < (Int_t)jet_constituentsB.size(); j++)
                                                            {
                                                                Float_t jet_const_pt  = jet_constituentsB[j].perp();
                                                                Float_t jet_const_phi = jet_constituentsB[j].phi();
                                                                Float_t jet_const_eta = jet_constituentsB[j].eta();
                                                                Int_t   user_index    = jet_constituentsB[j].user_index();
                                                                cout << "track: " << j << ", user_index: " << user_index << ", pt: " << jet_const_pt << endl;
                                                            }
                                                        }
#endif
                                                        //------------------


                                                        //cout << "i_trigger_use: " << i_trigger_use << ", trigger_index_min: " << trigger_index_min << ", iB: " << iB << ", random_trigger_index: " << random_trigger_index << endl;
                                                        Float_t jet_delta_eta = fabs(jet_etaB + vec_trigger_tracks[random_trigger_index][0]);
                                                        Float_t jet_delta_phi = fabs(jet_phiB - vec_trigger_tracks[random_trigger_index][1]);
                                                        if(jet_delta_phi > 2.0*Pi)  jet_delta_phi -= 2.0*Pi;

#if 0
                                                        if(i_trigger_use == 0 && lead_pt_bin == 9)
                                                        {
                                                            cout << "random_trigger_index: " << random_trigger_index
                                                                << ", i: " << i << ", iB: " << iB << ", event: " << counter << ", i_orig_smear: " << i_orig_smear << endl;
                                                        }
#endif
                                                        Float_t dijet_delta_phi = jet_phiB-vec_trigger_tracks[random_trigger_index][1]; // -2*Pi..2*Pi
                                                        if(dijet_delta_phi < 0.0) dijet_delta_phi += 2.0*Pi; // 0..2*Pi
                                                        if(dijet_delta_phi > 1.5*Pi)
                                                        {
                                                            dijet_delta_phi = -0.5*Pi + (dijet_delta_phi-1.5*Pi); // -0.5*Pi..1.5*Pi
                                                        }

                                                        // WWW
                                                        if(lead_pt_bin == 9 && i_trigger_use == 0 && jet_delta_eta < jet_delta_eta_cut && jet_areaB > array_areas[4])
                                                        {
                                                            if(jet_pt_subB < -5.6)
                                                            {
                                                                h_delta_phi_trigger_recoil[i_orig_smear]   ->Fill(dijet_delta_phi);
                                                                h_delta_phi_trigger_recoil[i_orig_smear+2] ->Fill(dijet_delta_phi,Et_weight_factor*PYTHIA_hard_bin_index_weight);

                                                                //if(fabs(Pi-jet_delta_phi) < jet_delta_phi_cut)
                                                                {
                                                                    h_eta_recoil[i_orig_smear]                 ->Fill(jet_etaB);
                                                                    h_phi_recoil[i_orig_smear]                 ->Fill(jet_phiB);
                                                                    h_eta_recoil[i_orig_smear+2]               ->Fill(jet_etaB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                    h_phi_recoil[i_orig_smear+2]               ->Fill(jet_phiB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                    if(iB == 0)
                                                                    {
                                                                        h_eta_trigger[i_orig_smear]   ->Fill(vec_trigger_tracks[random_trigger_index][0]);
                                                                        h_phi_trigger[i_orig_smear]   ->Fill(vec_trigger_tracks[random_trigger_index][1]);
                                                                        h_eta_trigger[i_orig_smear+2] ->Fill(vec_trigger_tracks[random_trigger_index][0],Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                        h_phi_trigger[i_orig_smear+2] ->Fill(vec_trigger_tracks[random_trigger_index][1],Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                    }
                                                                }
                                                            }
                                                        }

                                                        //cout << "leading_phi = " << leading_phi << ", jet_phiB = " << jet_phiB << ", dijet_delta_phi = " << dijet_delta_phi << endl;
                                                        //cout << "phiA = " << jet_phi << ", phiB =" << jet_phiB << ", jet_delta_phi = " << jet_delta_phi
                                                        //    << ", etaA = " << jet_eta << ", etaB = " << jet_etaB << ", jet_delta_eta = " << jet_delta_eta << endl;

                                                        //--------------------------------------------------
                                                        // Di-jet histograms
                                                        if(
                                                           jet_delta_eta < jet_delta_eta_cut
                                                          )
                                                        {
                                                            for(Int_t i_area = 0; i_area < N_jet_areas; i_area++)
                                                            {
                                                                // Fill spectra
                                                                if(jet_areaB > array_areas[i_area])
                                                                {
                                                                    h2D_dijet_pt_sub[i_area][N_leading_pt_bins][i_orig_smear+trigger_orig_smear] ->Fill(dijet_delta_phi,jet_pt_subB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                    h2D_dijet_pt[i_area][N_leading_pt_bins][i_orig_smear+trigger_orig_smear]     ->Fill(dijet_delta_phi,jet_ptB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                    if(lead_pt_bin != -1)
                                                                    {
                                                                        //-----------------------------------
                                                                        // For test purposes
                                                                        if(i_area == 4 && lead_pt_bin == 9 && i_orig_smear == 1 && trigger_orig_smear == 0)
                                                                        {
                                                                            h2D_dijet_pt_sub_test        ->Fill(dijet_delta_phi,jet_pt_subB);
                                                                            h2D_dijet_pt_sub_weight_test ->Fill(dijet_delta_phi,jet_pt_subB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                        }
                                                                        //-----------------------------------
                                                                        h2D_dijet_pt_sub[i_area][lead_pt_bin][i_orig_smear+trigger_orig_smear]   ->Fill(dijet_delta_phi,jet_pt_subB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                        h2D_dijet_pt[i_area][lead_pt_bin][i_orig_smear+trigger_orig_smear]       ->Fill(dijet_delta_phi,jet_ptB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                    }
                                                                    dijets_per_event[i_area]++;
                                                                }
                                                            }
                                                        }
                                                        //--------------------------------------------------



                                                        //--------------------------------------------------
                                                        // Recoil jet histograms
                                                        // M2S
                                                        if(
                                                           jet_delta_eta < jet_delta_eta_cut
                                                           && fabs(Pi-jet_delta_phi) < (25.0*(2.0*Pi/360.0))
                                                           // && iB != i // Don't use the same jet candidate twice
                                                          )
                                                        {
                                                            if(eMode == 310 && trigger_orig_smear == 0 && lead_pt_bin == 9)
                                                            {
                                                                vector<PseudoJet> jet_constituentsB = jets[i_orig_smear][iB].constituents();
                                                                for(Int_t const_iB = 0; const_iB < (Int_t)jet_constituentsB.size(); const_iB++)
                                                                {
                                                                    Int_t   user_indexB    = jet_constituentsB[const_iB].user_index();
                                                                    Float_t jet_const_ptB  = jet_constituentsB[const_iB].perp();

                                                                    cout << "jet_const_ptB: " << jet_const_ptB << endl;
                                                                    if(jet_const_ptB > 0.05)
                                                                    {
                                                                        Double_t m2B       = particle_info[0][user_indexB];
                                                                        Double_t dcaB      = particle_info[1][user_indexB];
                                                                        Double_t nSPiB     = particle_info[2][user_indexB];
                                                                        Double_t nSKB      = particle_info[3][user_indexB];
                                                                        Double_t nSPB      = particle_info[4][user_indexB];
                                                                        Double_t qpB       = particle_info[5][user_indexB];
                                                                        Double_t nhitsfitB = particle_info[6][user_indexB];
                                                                        Double_t etaB      = particle_info[7][user_indexB];
                                                                        Double_t phiB      = particle_info[8][user_indexB];

                                                                        Int_t chargeB      = 0;
                                                                        if(qpB < 0.0) chargeB = 1;

                                                                        Float_t jet_delta_eta_trackB = fabs(etaB - vec_trigger_tracks[random_trigger_index][0]);
                                                                        Float_t jet_delta_phi_trackB = fabs(phiB - vec_trigger_tracks[random_trigger_index][1]);

                                                                        Float_t delta_phiB = fabs(Pi-jet_delta_phi_trackB);

                                                                        Double_t radial_distance = TMath::Sqrt(TMath::Power(delta_phiB,2.0) + TMath::Power(jet_delta_eta_trackB,2.0));

#if 0
                                                                        cout << "event: " << counter << ", lead_pt_bin: " << lead_pt_bin << ", lead_pt: " << vec_trigger_tracks[random_trigger_index][2] << ", i: " << i << ", iB: " << iB << ", i_orig_smear: " << i_orig_smear
                                                                            << ", random_trigger_index: " << random_trigger_index
                                                                            << ", i_trigger_use: " << i_trigger_use << ", trigger_index_min: " << trigger_index_min << endl;
#endif
#if 0
                                                                        cout << "event: " << counter << ", i_orig_smear: " << i_orig_smear
                                                                            << ", random_trigger_index:" << random_trigger_index
                                                                            <<  ", phi_track: " << phiB << ", phi_trigger: " << vec_trigger_tracks[random_trigger_index][1]
                                                                            << ", delta_phiB: " << delta_phiB << ", eta_track: " << etaB
                                                                            << ", eta_trigger: " << vec_trigger_tracks[random_trigger_index][0]
                                                                            << ", delta_eta: " << jet_delta_eta_trackB
                                                                            << ", radial_distance: " << radial_distance << ", m2B: " << m2B << endl;
#endif
                                                                        if(radial_distance > Array_m2_radial_dist[0] || m2B < -99.0) continue;

                                                                        for(Int_t i_radial = 0; i_radial < N_m2_radial_dist; i_radial++)
                                                                        {
                                                                            if(radial_distance < Array_m2_radial_dist[i_radial])
                                                                            {
                                                                                cout << "radial_distance: " << radial_distance << ", jet_const_ptB: " << jet_const_ptB << ", m2B: " << m2B << endl;
                                                                                h2D_m2_vs_pT[i_radial][chargeB] ->Fill(jet_const_ptB,m2B);
                                                                            }
                                                                        }

                                                                    }
                                                                }
                                                            }


                                                            //cout << "Et_weight_factor: " << Et_weight_factor << endl;
                                                            h_jet_area_array[ez_bin][emult_bin][ePsi_bin]    ->Fill(jet_areaB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                            h_jet_rhoarea_array[ez_bin][emult_bin][ePsi_bin] ->Fill(jet_areaB*jet_rho,Et_weight_factor*PYTHIA_hard_bin_index_weight);

                                                            if(i_orig_smear == 0 && trigger_orig_smear == 0)
                                                            {
                                                                //if(jet_rho > 0.001) cout << "iB: " << iB << ", rho: " << jet_rho << ", jet_areaB: " << jet_areaB << ", jet_pt_subB: " << jet_pt_subB << endl;
                                                                h_area_vs_recoil_jet_pt[N_leading_pt_bins] ->Fill(jet_pt_subB,jet_areaB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                h_rho[N_leading_pt_bins]                   ->Fill(jet_rho,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                if(lead_pt_bin != -1)
                                                                {
                                                                    h_area_vs_recoil_jet_pt[lead_pt_bin] ->Fill(jet_pt_subB,jet_areaB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                    h_rho[lead_pt_bin]                   ->Fill(jet_rho,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                }
                                                            }

                                                            for(Int_t i_area = 0; i_area < N_jet_areas; i_area++)
                                                            {
                                                                // Fill spectra
                                                                if(jet_areaB > array_areas[i_area])
                                                                {
                                                                    //cout << "jet_ptB = " << jet_ptB << endl;

                                                                    N_accepted_jets[i_area]++;

                                                                    //cout << "Weight factor: " << Et_weight_factor*PYTHIA_hard_bin_index_weight << ", Et_weight_factor: " << Et_weight_factor << ", PYTHIA_hard_bin_index_weight: " << PYTHIA_hard_bin_index_weight << endl;
                                                                    h_jet_pt_sub[i_area][N_leading_pt_bins][i_orig_smear+trigger_orig_smear] ->Fill(jet_pt_subB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                    h_jet_pt[i_area][N_leading_pt_bins][i_orig_smear+trigger_orig_smear]     ->Fill(jet_ptB,Et_weight_factor*PYTHIA_hard_bin_index_weight);

                                                                    if(lead_pt_bin != -1)
                                                                    {
                                                                        //cout << "Fill spectra, PYTHIA_hard_bin_index_weight: " << PYTHIA_hard_bin_index_weight << endl;
                                                                        h_jet_pt_sub[i_area][lead_pt_bin][i_orig_smear+trigger_orig_smear]   ->Fill(jet_pt_subB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                        h_jet_pt[i_area][lead_pt_bin][i_orig_smear+trigger_orig_smear]       ->Fill(jet_ptB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
#if 0
                                                                        if(i_area == 4 && i_trigger_use == 0 && lead_pt_bin == 9) cout << "Event: " << counter
                                                                            << ", i_orig_smear: " << i_orig_smear << ", jet_etaB: " << jet_etaB << ", jet_phiB: " << jet_phiB << ", jet_ptB: " << jet_ptB
                                                                                << ", jet_pt_subB: " << jet_pt_subB << ", jet_areaB: " << jet_areaB << ", jet_rho: " << jet_rho << endl;
#endif

                                                                    }
                                                                    h_jet_area[i_area][i_orig_smear+trigger_orig_smear]                      ->Fill(jet_areaB,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                    //if(i_area == 4) cout << "i_orig_smear: " << i_orig_smear << ", i_area: " << i_area << ", lead_pt_bin: " << lead_pt_bin << ", jet_pt_subB: " << jet_pt_subB << ", jets_per_event: " << jets_per_event[i_area] << endl;
                                                                    jets_per_event[i_area]++;
                                                                }
                                                            }
                                                        }
                                                        //--------------------------------------------------


                                                    } // end of jet B loop

                                                    if(i_trigger_use == 0)
                                                    {
                                                        for(Int_t i_area_acc = 0; i_area_acc < N_jet_areas; i_area_acc++)
                                                        {
                                                            h_N_accepted_recoil_jets[i_area_acc][lead_pt_bin]       ->Fill(N_accepted_jets[i_area_acc]);
                                                            h_N_accepted_recoil_jets[i_area_acc][N_leading_pt_bins] ->Fill(N_accepted_jets[i_area_acc]);
#if 0
                                                            if(i_area_acc == 4) cout << "event: " << counter << ", i_area_acc: " << i_area_acc << ", N_recoil_jets: " << N_accepted_jets[i_area_acc]
                                                                << ", jetA: " << i << ", lead_pt_bin: " << lead_pt_bin << ", i_orig_smear: " << i_orig_smear << endl;
#endif
                                                        }
                                                    }

                                                } // end of trigger loop


                                                //----------------------------------------------------------------
                                                if(eMode == 3)
                                                {
                                                    //cout << "leading_pt = " << leading_pt << ", lead_pt_bin = " << lead_pt_bin << endl;

                                                    h_jet_area_array[ez_bin][emult_bin][ePsi_bin]    ->Fill(jet_area,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                    h_jet_rhoarea_array[ez_bin][emult_bin][ePsi_bin] ->Fill(jet_area*jet_rho,Et_weight_factor*PYTHIA_hard_bin_index_weight);

                                                    for(Int_t i_area = 0; i_area < N_jet_areas; i_area++)
                                                    {
                                                        if(jet_area > array_areas[i_area])
                                                        {
                                                            h_jet_pt_sub[i_area][N_leading_pt_bins][i_orig_smear] ->Fill(jet_pt_sub,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                            h_jet_pt[i_area][N_leading_pt_bins][i_orig_smear]     ->Fill(jet_pt,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                            if(lead_pt_bin != -1)
                                                            {
                                                                h_jet_pt_sub[i_area][lead_pt_bin][i_orig_smear]   ->Fill(jet_pt_sub,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                                h_jet_pt[i_area][lead_pt_bin][i_orig_smear]       ->Fill(jet_pt,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                            }
                                                            h_jet_area[i_area][i_orig_smear]                      ->Fill(jet_area,Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                                            jets_per_event[i_area]++;
                                                        }
                                                    }
                                                }
                                                // end eMode = 3
                                                //----------------------------------------------------------------


                                            } // end of trigger track loop
                                        }
                                    }
                                    //----------------------------------------------------------------

                                }
                                // End jet A loop
                                //----------------------------------------------------------------------------------------


                            }


                            for(Int_t k = 0; k < N_jet_areas; k++)
                            {
                                h_jet_per_event[k][i_orig_smear]   ->Fill(jets_per_event[k],Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                h_dijet_per_event[k][i_orig_smear] ->Fill(dijets_per_event[k],Et_weight_factor*PYTHIA_hard_bin_index_weight);
                                //cout << "area: " << k << ", i_orig_smear: " << i_orig_smear << ", jet_per_event[area]: " << jets_per_event[k] << endl;
                            }
                        } // end i_orig_smear
                        vec_flag_low_pt_PYTHIA.clear();
                        vec_flag_low_pt_PYTHIA_jet_pt.clear();


                        //----------------------------------
                        // Calculate the jet matching/reconstrution efficiencies - only for mode 311 -> PYTHIA
                        if(eMode == 311 || eMode == 312)
                        {
                            //-------------------------------------------------
                            // Jet matching cuts
                            const Double_t min_fraction_of_matched_pt      = 0.15; // require at least 15% of matched track pt, otherwise count it as a not matched jet --> inefficiency
                            const Double_t max_ratio_jet_pt                = 1.5;  // avoid merged jets (real) into one eff jet -> if ratio eff jet / real jet > 1.5 skip it --> inefficiency
                            //-------------------------------------------------


                            //cout << "rho[0]: " << jet_rho_array[0] << ", rho[1]: " << jet_rho_array[1] << endl;
                            if(jet_rho_array[0] > 0.0 && jet_rho_array[1] > 0.0)
                            {
                                Int_t N_jets_PYTHIA[2] = {(Int_t)vec_jet_array[0].size(),(Int_t)vec_jet_array[1].size()};
#if 0
                                cout << "" << endl;
                                cout << "----------------------------------" << endl;
                                cout << "PYTHIA jets: " << N_jets_PYTHIA[0] << ", PYTHIA jets (smeared+eff): " << N_jets_PYTHIA[1] <<
                                    ", rho: " << jet_rho_array[0] << ", rho(smeared+eff): " << jet_rho_array[1] << endl;
#endif

                                //----------------------------------------------------------------------------------------
                                // Vectors for smeared jets
                                std::vector< std::vector<Int_t> > Index_smeared_jets_matched_not_matched;
                                std::vector< std::vector< std::vector<Double_t> > > Vec_info_smeared_jets;
                                const Int_t N_variables_jet = 11;
                                Index_smeared_jets_matched_not_matched.resize(N_variables_jet);
                                Vec_info_smeared_jets.resize(N_variables_jet);
                                for(Int_t i_index_matched = 0; i_index_matched < N_variables_jet; i_index_matched++)
                                {
                                    Index_smeared_jets_matched_not_matched[i_index_matched].resize(N_jets_PYTHIA[1]);
                                    Vec_info_smeared_jets[i_index_matched].resize(N_jets_PYTHIA[1]);
                                    for(Int_t i_jet_smeared = 0; i_jet_smeared < N_jets_PYTHIA[1]; i_jet_smeared++)
                                    {
                                        Index_smeared_jets_matched_not_matched[i_index_matched][i_jet_smeared] = 0;
                                    }
                                }
                                //----------------------------------------------------------------------------------------



                                //----------------------------------------------------------------------------------------
                                // Vectors for original jets
                                std::vector< std::vector<Double_t> > Vec_info_original_jets;
                                Vec_info_original_jets.resize(N_variables_jet);
                                for(Int_t i_index_matched = 0; i_index_matched < N_variables_jet; i_index_matched++)
                                {
                                    Vec_info_original_jets[i_index_matched].resize(N_jets_PYTHIA[0]);
                                }
                                //----------------------------------------------------------------------------------------


                                // Loop over the original PYTHIA jets
                                for(Int_t i_jet = 0; i_jet < N_jets_PYTHIA[0]; i_jet++)
                                {
                                    Float_t jet_pt     = vec_jet_array[0][i_jet][0];
                                    Float_t jet_area   = vec_jet_array[0][i_jet][1];
                                    Float_t jet_pt_sub = vec_jet_array[0][i_jet][2];
                                    Float_t jet_eta    = vec_jet_array[0][i_jet][3];
                                    Float_t jet_phi    = vec_jet_array[0][i_jet][4];

                                    //cout << "Matching, random_trigger_index[original]: " << random_trigger_index_array[0][9] << ", random_trigger_index[eff]: " << random_trigger_index_array[1][9] << endl;
                                    Float_t jet_delta_eta_orig = -1.0;
                                    Float_t jet_delta_phi_orig = -1.0;
                                    random_trigger_index = random_trigger_index_array[0][9];
                                    if(random_trigger_index >= 0)
                                    {
                                        if((Int_t)vec_trigger_tracks.size() <= 0) cout << "vec_trigger_tracks has size 0" << endl;
                                        if(random_trigger_index >= (Int_t)vec_trigger_tracks.size()) cout << "random_trigger_index larger than vec_trigger_tracks size" << endl;
                                        jet_delta_eta_orig = fabs(jet_eta + vec_trigger_tracks[random_trigger_index][0]);
                                        jet_delta_phi_orig = fabs(jet_phi - vec_trigger_tracks[random_trigger_index][1]);
                                    }

                                    // Loop over the original PYTHIA jet constituent tracks
                                    Float_t jet_sum_track_pt = 0.0;
                                    vector<Int_t> vec_user_index;
                                    for(Int_t j_track = 0; j_track < (Int_t)vec_jet_const_array[0][i_jet].size(); j_track++)
                                    {
                                        Float_t jet_const_pt  = vec_jet_const_array[0][i_jet][j_track][0];
                                        if(jet_const_pt < 0.15) continue; // remove ghost tracks
                                        jet_sum_track_pt += jet_const_pt; // jet pt based only on the tracks which are later on matched
                                        Float_t jet_const_phi = vec_jet_const_array[0][i_jet][j_track][1];
                                        Float_t jet_const_eta = vec_jet_const_array[0][i_jet][j_track][2];
                                        Int_t   user_index    = (Int_t)vec_jet_const_array[0][i_jet][j_track][3];
                                        vec_user_index.push_back(user_index);
                                    }

#if 0
                                    cout << "PYTHIA jet: " << i_jet << ", jet_pT: " << jet_pt << ", tracks: " << (Int_t)vec_user_index.size() << endl;
#endif
                                    // Loop over the smeared+eff PYTHIA jets
                                    vector<Int_t> vec_index_counter; // stores the number of matched tracks for every smeared+eff PYTHIA jet
                                    //std::vector< std::vector<Double_t> > vec_match_jet_info;
                                    //vec_match_jet_info.resize(N_jets_PYTHIA[1]);
                                    Double_t best_fraction_of_matched_tracks = -999.0;  // keeps always the best fraction of matched tracks
                                    Double_t best_jet_pt_smeared             = -999.0;  // keeps always the highest matched jet pT
                                    Double_t best_jet_area_smeared           = -999.0;
                                    Double_t best_jet_pt_sub_smeared         = -999.0;
                                    Double_t best_jet_delta_phi              = -999.0;
                                    Int_t    Index_smeared_matched           = -999;
                                    Double_t best_jet_pt_match_fraction      = -999.0;
                                    Double_t best_jet_pt_ratio               = -999.0; // jet pt (eff)/ jet pt (original) -> low pt (original) can be matched to high pt (eff) (jet merging)


                                    // Loop over the momentum smeared and efficiency applied PYTHIA jets
                                    for(Int_t i_jet_smeared = 0; i_jet_smeared < N_jets_PYTHIA[1]; i_jet_smeared++)
                                    {
                                        Float_t jet_pt_smeared     = vec_jet_array[1][i_jet_smeared][0];
                                        Float_t jet_area_smeared   = vec_jet_array[1][i_jet_smeared][1];
                                        Float_t jet_pt_sub_smeared = vec_jet_array[1][i_jet_smeared][2];
                                        Float_t jet_eta_smeared    = vec_jet_array[1][i_jet_smeared][3];
                                        Float_t jet_phi_smeared    = vec_jet_array[1][i_jet_smeared][4];

                                        Float_t jet_delta_eta_smear = -1.0;
                                        Float_t jet_delta_phi_smear = -1.0;
                                        if(random_trigger_index >= 0)
                                        {
                                            jet_delta_eta_smear = fabs(jet_eta_smeared + vec_trigger_tracks[random_trigger_index][0]);
                                            jet_delta_phi_smear = fabs(jet_phi_smeared - vec_trigger_tracks[random_trigger_index][1]);
                                        }



#if 0
                                        if(i_jet == 0)
                                        {
                                            cout << "i_jet_smeared: " << i_jet_smeared << ", rho: " << jet_rho_array[1] << ", area: " << jet_area_smeared
                                                << ", rho*area: " << jet_rho_array[1]*jet_area_smeared << ", jet_pt_smeared: " << jet_pt_smeared << ", jet_pt_sub_smeared: " << jet_pt_sub_smeared << endl;
                                        }
#endif

                                        // Loop over the smeared+eff PYTHIA jet constituent tracks
                                        Int_t Index_counter = 0; // counts for this particular original PYTHIA jet - smeared+eff PYTHIA jet combination how many tracks are matched
                                        Float_t jet_sum_track_pt_smeared = 0.0;
                                        for(Int_t j_track_smeared = 0; j_track_smeared < (Int_t)vec_jet_const_array[1][i_jet_smeared].size(); j_track_smeared++)
                                        {
                                            Float_t jet_const_pt_smeared  = vec_jet_const_array[1][i_jet_smeared][j_track_smeared][0];
                                            if(jet_const_pt_smeared < 0.15) continue; // remove ghost tracks
                                            Float_t jet_const_phi_smeared = vec_jet_const_array[1][i_jet_smeared][j_track_smeared][1];
                                            Float_t jet_const_eta_smeared = vec_jet_const_array[1][i_jet_smeared][j_track_smeared][2];
                                            Int_t   user_index_smeared    = (Int_t)vec_jet_const_array[1][i_jet_smeared][j_track_smeared][3];

                                            // check if that user_index is in original PYTHIA jet
                                            for(Int_t i_index = 0; i_index < (Int_t)vec_user_index.size(); i_index++)
                                            {
                                                if(user_index_smeared == vec_user_index[i_index]) // same track was found
                                                {
                                                    Index_counter++;
                                                    jet_sum_track_pt_smeared += jet_const_pt_smeared; // use only matched tracks for the matched pt comparison
                                                    //cout << "user_index_smeared: " << user_index_smeared << endl;
                                                }
                                            }
                                        }

                                        vec_index_counter.push_back(Index_counter);
                                        Double_t fraction_of_matched_tracks = -9999.0;
                                        Double_t fraction_of_matched_pt     = -9999.0;
                                        Double_t jet_pt_ratio               = -9999.0;
                                        if((Int_t)vec_user_index.size() > 0 && Index_counter > 0 && jet_sum_track_pt > 0.0 && jet_pt > 0.0)
                                        {
                                            fraction_of_matched_tracks = ((Double_t)Index_counter)/((Double_t)vec_user_index.size());
                                            fraction_of_matched_pt     = jet_sum_track_pt_smeared/jet_sum_track_pt;
                                            jet_pt_ratio               = jet_pt_smeared/jet_pt;

#if 0

                                            cout << "PYTHIA jet (smeared+eff): " << i_jet_smeared << ", fraction of matched tracks: " << fraction_of_matched_tracks
                                                << ", matched tracks: " << Index_counter << ", tracks in original jet: " << (Int_t)vec_user_index.size() << ", pt-rho*A (smeared+eff): " << jet_pt_sub_smeared <<
                                                ", original jet pT: " << jet_sum_track_pt << ", matched tracks jet pT: " << jet_sum_track_pt_smeared << endl;
#endif
                                        }

                                        // Find the reconstructed jet with most matched tracks and highest matched energy
                                        if(
                                           //fraction_of_matched_tracks  >= best_fraction_of_matched_tracks
                                           fraction_of_matched_tracks  > 0
                                           && jet_sum_track_pt_smeared >= best_jet_pt_smeared
                                          )
                                        {
                                            best_fraction_of_matched_tracks = fraction_of_matched_tracks;
                                            best_jet_pt_smeared             = jet_sum_track_pt_smeared;
                                            best_jet_area_smeared           = jet_area_smeared;
                                            best_jet_pt_sub_smeared         = jet_pt_sub_smeared;
                                            best_jet_delta_phi              = jet_delta_phi_smear;
                                            Index_smeared_matched           = i_jet_smeared;
                                            best_jet_pt_match_fraction      = fraction_of_matched_pt;
                                            best_jet_pt_ratio               = jet_pt_ratio;
                                            //best_jet_pt_match_fraction

#if 0
                                            cout << "Index (eff): " << i_jet_smeared << ", match: " << i_jet
                                                << ", frac. match. tracks: " << fraction_of_matched_tracks
                                                << ", match. tracks: " << Index_counter << ", tracks in original jet: " << (Int_t)vec_user_index.size()
                                                << ", pt-rho*A (smeared+eff): " << jet_pt_sub_smeared <<
                                                ", original jet pT: " << jet_sum_track_pt << ", frac. match. jet pT: " << fraction_of_matched_pt << endl;
#endif
                                        }

                                        //vec_match_jet_info[i_jet_smeared].push_back(jet_pt_sub,jet_pt_sub_smeared,fraction_of_matched_tracks);
                                    } // End of detector level loop


                                    h_1D_det_diff_part_lvl_area      ->Fill(jet_area-best_jet_area_smeared);
                                    h_2D_det_vs_part_lvl_area        ->Fill(jet_area,best_jet_area_smeared);
                                    if(best_jet_pt_smeared > 5.0)
                                    {
                                        h_1D_det_diff_part_lvl_area_lt5GeV ->Fill(jet_area-best_jet_area_smeared);
                                        h_2D_det_vs_part_lvl_area_lt5GeV   ->Fill(jet_area,best_jet_area_smeared);
                                    }

                                    if(Index_smeared_matched >= 0)
                                    {
                                        // Count how often a particular PYTHIA smeared + efficiency applied jet was matched to a original PYTHIA jet
                                        Index_smeared_jets_matched_not_matched[0][Index_smeared_matched]++;

                                        // Sub information of smeared matched jets
                                        Vec_info_smeared_jets[0][Index_smeared_matched].push_back(i_jet);
                                        Vec_info_smeared_jets[1][Index_smeared_matched].push_back(best_fraction_of_matched_tracks);
                                        Vec_info_smeared_jets[2][Index_smeared_matched].push_back(best_jet_pt_match_fraction);
                                        Vec_info_smeared_jets[3][Index_smeared_matched].push_back(best_jet_pt_ratio);
                                        Vec_info_smeared_jets[4][Index_smeared_matched].push_back(best_jet_delta_phi);
                                        Vec_info_smeared_jets[5][Index_smeared_matched].push_back(best_jet_pt_smeared);
                                        Vec_info_smeared_jets[6][Index_smeared_matched].push_back(best_jet_pt_sub_smeared);
                                        Vec_info_smeared_jets[7][Index_smeared_matched].push_back(best_jet_area_smeared);
                                        Vec_info_smeared_jets[8][Index_smeared_matched].push_back(jet_sum_track_pt);
                                        Vec_info_smeared_jets[9][Index_smeared_matched].push_back(jet_delta_phi_orig);
                                        Vec_info_smeared_jets[10][Index_smeared_matched].push_back(jet_area);

                                        Float_t jet_pt     = vec_jet_array[0][i_jet][0];
                                        Float_t jet_area   = vec_jet_array[0][i_jet][1];
                                        Float_t jet_pt_sub = vec_jet_array[0][i_jet][2];
                                        Float_t jet_eta    = vec_jet_array[0][i_jet][3];
                                        Float_t jet_phi    = vec_jet_array[0][i_jet][4];
                                    }

                                    // Sub information of original matched jets
                                    Vec_info_original_jets[0][i_jet]  = Index_smeared_matched;
                                    Vec_info_original_jets[1][i_jet]  = best_fraction_of_matched_tracks;
                                    Vec_info_original_jets[2][i_jet]  = best_jet_pt_match_fraction;
                                    Vec_info_original_jets[3][i_jet]  = best_jet_pt_ratio;
                                    Vec_info_original_jets[4][i_jet]  = jet_delta_phi_orig;
                                    Vec_info_original_jets[5][i_jet]  = jet_pt;
                                    Vec_info_original_jets[6][i_jet]  = jet_pt_sub;
                                    Vec_info_original_jets[7][i_jet]  = jet_area;
                                    Vec_info_original_jets[8][i_jet]  = jet_sum_track_pt;
                                    Vec_info_original_jets[9][i_jet]  = best_jet_delta_phi;
                                    Vec_info_original_jets[10][i_jet] = best_jet_area_smeared;

                                    Double_t Delta_pt = best_jet_pt_sub_smeared - jet_pt; // reconstructed jet pT - embedded jet pT

                                    for(Int_t i_area = 0; i_area < N_jet_areas; i_area++)
                                    {
                                        //cout << "Original PYTHIA jet: " << i_jet << endl;
                                        if(best_jet_area_smeared > array_areas[i_area])
                                        {
                                            // Fill histogram for matched pT
                                            if(eMode == 312)
                                            {
                                                h_Delta_pt_vs_embed_pt[i_area][0][ePsi_bin] ->Fill(jet_pt,Delta_pt,PYTHIA_hard_bin_index_weight);
#if 0
                                                if(i_area == 0) cout << "jet_pt (embed): " << jet_pt << ", rec jet_pt (matched): " << best_jet_pt_smeared << ", Delta_pt: " << Delta_pt << ", best_jet_pt_sub_smeared: "
                                                    << best_jet_pt_sub_smeared <<  endl;
#endif
                                            }
                                        }
                                    }
#if 0
                                    cout << "i_jet: " << i_jet << ", N_tracks: " << (Int_t)vec_user_index.size() << ", best fraction of matched track: "
                                        << best_fraction_of_matched_tracks << ", matched pT: " << best_jet_pt_smeared << ", out of: " << jet_sum_track_pt << endl;
#endif
                                } // End of particle level loop
                                //----------------------------------
#if 0
                                cout << "" << endl;
                                cout << "----------------------" << endl;
#endif
                                // Loop over original jet information --> for final efficiency
                                for(Int_t i_jet = 0; i_jet < N_jets_PYTHIA[0]; i_jet++)
                                {
#if 0
                                    cout << "Original jet: " << i_jet << ", matched with smeared jet: " << Vec_info_original_jets[0][i_jet]
                                        << ", % tracks: " << Vec_info_original_jets[1][i_jet]
                                        << ", % pt: " << Vec_info_original_jets[2][i_jet]
                                        << ", jet pt(eff)/jet pt(orig.): " << Vec_info_original_jets[3][i_jet] << endl;

#endif
                                    if(fabs(Pi-Vec_info_original_jets[4][i_jet]) < jet_delta_phi_cut) // check if original Delta phi is in acceptance
                                    {
                                        for(Int_t i_area = 0; i_area < N_jet_areas; i_area++)
                                        {
                                            if(Vec_info_original_jets[7][i_jet] > array_areas[i_area]) // check if original track area is accepted
                                            {
                                                h_Sim_original_pT_all[i_area] ->Fill(Vec_info_original_jets[8][i_jet]);
                                                if( // check if smeared Delta phi and area are accepted
                                                   fabs(Pi-Vec_info_original_jets[9][i_jet]) < jet_delta_phi_cut // only smeared part
                                                   && Vec_info_original_jets[10][i_jet]      > array_areas[i_area] // only smeared part
                                                   && Vec_info_original_jets[3][i_jet]       < max_ratio_jet_pt // avoid merged jets (real) into one eff jet -> if ratio eff jet / real jet > 1.5 skip it --> inefficiency
                                                   && Vec_info_original_jets[2][i_jet]       > min_fraction_of_matched_pt // require at least 15% of matched track pt, otherwise count it as a not matched jet --> inefficiency
                                                  )
                                                {
                                                    h_Sim_original_pT_cut[i_area] ->Fill(Vec_info_original_jets[8][i_jet]); // Some jets are completely lost -> final correction all/cut
                                                }
                                            }
                                        }
                                    }
                                }
#if 0
                                cout << "----------------------" << endl;
                                cout << "" << endl;

                                cout << "----------------------" << endl;
#endif
                                // Loop over smeared jet information
                                for(Int_t i_jet_smeared = 0; i_jet_smeared < N_jets_PYTHIA[1]; i_jet_smeared++)
                                {
#if 0
                                    cout << "Index (eff): " << i_jet_smeared << ", match counter PYTHIA smeared: " << Index_smeared_jets_matched_not_matched[0][i_jet_smeared] << endl;
#endif

                                    // Check all smeared jets which were not matched to a real jet
                                    if(Index_smeared_jets_matched_not_matched[0][i_jet_smeared] == 0) // no real jet matched to this smeared jet
                                    {
                                        Float_t jet_pt_smeared     = vec_jet_array[1][i_jet_smeared][0];
                                        Float_t jet_area_smeared   = vec_jet_array[1][i_jet_smeared][1];
                                        Float_t jet_pt_sub_smeared = vec_jet_array[1][i_jet_smeared][2];
                                        Float_t jet_eta_smeared    = vec_jet_array[1][i_jet_smeared][3];
                                        Float_t jet_phi_smeared    = vec_jet_array[1][i_jet_smeared][4];

                                        Float_t jet_delta_eta_smear = -1.0;
                                        Float_t jet_delta_phi_smear = -1.0;
                                        if(random_trigger_index >= 0)
                                        {
                                            jet_delta_eta_smear = fabs(jet_eta_smeared + vec_trigger_tracks[random_trigger_index][0]);
                                            jet_delta_phi_smear = fabs(jet_phi_smeared - vec_trigger_tracks[random_trigger_index][1]);
                                        }

                                        if(fabs(Pi-jet_delta_phi_smear) < jet_delta_phi_cut) // check if smeared Delta phi is in acceptance
                                        {
                                            for(Int_t i_area = 0; i_area < N_jet_areas; i_area++)
                                            {
                                                if(jet_area_smeared > array_areas[i_area])
                                                {
                                                    h_Sim_smeared_pT_not_matched[i_area] ->Fill(jet_pt_smeared);
                                                }
                                            }
                                        }
                                    }

                                    // Sub information for multiple matched real jets
                                    for(Int_t i_real_index = 0; i_real_index < (Int_t)Vec_info_smeared_jets[0][i_jet_smeared].size(); i_real_index++)
                                    {
#if 0
                                        cout << "    PYTHIA original index: " << Vec_info_smeared_jets[0][i_jet_smeared][i_real_index]
                                            << ", % tracks: " << Vec_info_smeared_jets[1][i_jet_smeared][i_real_index]
                                            << ", % pt: " << Vec_info_smeared_jets[2][i_jet_smeared][i_real_index]
                                            << ", jet pt(eff)/jet pt(orig.): " << Vec_info_smeared_jets[3][i_jet_smeared][i_real_index] << endl;
#endif

                                        if(fabs(Pi-Vec_info_smeared_jets[4][i_jet_smeared][i_real_index]) < jet_delta_phi_cut) // check if smeared Delta phi is in acceptance
                                        {
                                            for(Int_t i_area = 0; i_area < N_jet_areas; i_area++)
                                            {
                                                if(Vec_info_smeared_jets[7][i_jet_smeared][i_real_index]    > array_areas[i_area]
                                                   && Vec_info_smeared_jets[3][i_jet_smeared][i_real_index] < max_ratio_jet_pt // avoid merged jets (real) into one eff jet -> if ratio eff jet / real jet > 1.5 skip it --> inefficiency
                                                   && Vec_info_smeared_jets[2][i_jet_smeared][i_real_index] > min_fraction_of_matched_pt // require at least 15% of matched track pt, otherwise count it as a not matched jet --> inefficiency
                                                  )
                                                {
                                                    // Fill histogram for matched pT
                                                    h2D_Sim_matched_pT_vs_original_pT[i_area]          ->Fill(Vec_info_smeared_jets[8][i_jet_smeared][i_real_index],Vec_info_smeared_jets[5][i_jet_smeared][i_real_index]);
                                                    h2D_Sim_original_pT_vs_matched_pT[i_area]          ->Fill(Vec_info_smeared_jets[5][i_jet_smeared][i_real_index],Vec_info_smeared_jets[8][i_jet_smeared][i_real_index]);
                                                    h_matched_tracks_fraction[i_area]                  ->Fill(Vec_info_smeared_jets[1][i_jet_smeared][i_real_index]);
                                                    h2D_matched_tracks_fraction_vs_original_pT[i_area] ->Fill(Vec_info_smeared_jets[8][i_jet_smeared][i_real_index],Vec_info_smeared_jets[1][i_jet_smeared][i_real_index]);
                                                }
                                            }
                                        }
                                    }
                                }

#if 0
                                cout << "----------------------" << endl;
                                //----------------------------------

                                cout << "----------------------------------" << endl;
                                cout << "" << endl;
#endif
                            }
                            else
                            {
                                cout << "Jet matching couldn't be calculated, rho <= 0" << endl;
                            }
                        }
                        // End of PYTHIA jet reconstruction efficiency calculation
                        //----------------------------------



                    }
                    //-----------------------------------------------------------------------------

                }
            }
        }
    }

    if(eMode == 2) // eMode 2 is mixed event
    {
        cout << "Mixed event mode started" << endl;
        Int_t vertex_pos_counter = 0;
        Int_t i_SE_ME = 0;

        Long64_t stop_event_use_loop = stop_event_use;
        if(stop_event_use_loop > file_entries_SE_ME[i_SE_ME]) stop_event_use_loop = file_entries_SE_ME[i_SE_ME];

        Int_t ME_event_counter = 0;
        Two_dim_Int_vector ME_track_number_vector; // stores the track ids for every event which is mixed
        ME_track_number_vector.resize(N_max_events);

        for(Long64_t counter = start_event_use; counter < stop_event_use_loop; counter++)
        {
            if (counter != 0  &&  counter % 100 == 0)
                cout << "." << flush;
            if (counter != 0  &&  counter % 1000 == 0)
            {
                if((stop_event_use_loop-start_event_use) > 0)
                {
                    Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use_loop-start_event_use));
                    cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data " << flush;
                }
            }

            if (!input_SE_ME[i_SE_ME]->GetEntry( counter )) // take the event -> information is stored in event
                break;  // end of data chunk


            //-----------------------------------------------------------------------------
            // Event information
            Float_t prim_vertex_x   = JetTrackEvent->getx();
            Float_t prim_vertex_y   = JetTrackEvent->gety();
            Float_t prim_vertex_z   = JetTrackEvent->getz();
            Int_t   RunId           = JetTrackEvent->getid();
            Float_t refMult         = JetTrackEvent->getmult();
            Float_t n_prim          = JetTrackEvent->getn_prim();
            Float_t n_non_prim      = JetTrackEvent->getn_non_prim();
            Int_t   n_tofmatch_prim = JetTrackEvent->getn_tof_prim();
            Int_t   SE_ME_flag      = JetTrackEvent->getSE_ME_flag();
            Float_t ZDCx            = JetTrackEvent->getZDCx();
            Float_t BBCx            = JetTrackEvent->getBBCx();
            Float_t vzVPD           = JetTrackEvent->getvzVpd();
            Int_t   N_Particles     = JetTrackEvent->getNumParticle();
            TVector2 QvecEtaPos     = JetTrackEvent->getQvecEtaPos();
            TVector2 QvecEtaNeg     = JetTrackEvent->getQvecEtaNeg();
            Int_t    cent9          = JetTrackEvent->getcent9();
            ME_track_number_vector[ME_event_counter].resize(N_Particles);
            for(Int_t i_track = 0; i_track < N_Particles; i_track++)
            {
                ME_track_number_vector[ME_event_counter][i_track] = i_track;
            }

            if(fabs(prim_vertex_z) > z_acceptance[eBeamTimeNum]) continue;

            Double_t reweight = 1.0;
            if(eMode != 311) // not PYTHIA
            {
                //refmultCorrUtil->init(RunId); // ALEX
                // refmultCorrUtil->initEvent(refMult, prim_vertex_z, ZDCx); // ALEX

                // Get centrality bins
                //   see StRefMultCorr.h for the definition of centrality bins
                //erefMult_bin16   = refmultCorrUtil->getCentralityBin16();
                //erefMult_bin     = refmultCorrUtil->getCentralityBin9();

                // reweight = refmultCorrUtil->getWeight();  // ALEX
                reweight = 1.0; // ALEX
            }

            //cout << "ME_event_counter = " << ME_event_counter << ", array size = " << (Int_t)ME_track_number_vector[ME_event_counter].size() << endl;
            //-----------------------------------------------------------------------------



            //---------------------------------------------------------------------------
            // Good/Bad run selection
            Int_t flag_good_run = 1;
            for(Int_t bad_run = 0; bad_run < n_bad_run_numbers[eBeamTimeNum]; bad_run++)
            {
                if(bad_run_numbers[bad_run] == (Int_t)RunId)
                {
                    flag_good_run = 0;
                    break;
                }
            }
            if(flag_good_run == 0) continue;
            //---------------------------------------------------------------------------



            // Fill event information for jet
            JetTrackEvent_ME[ME_event_counter]->clearParticleList();
            JetTrackEvent_ME[ME_event_counter]->setx(prim_vertex_x);
            JetTrackEvent_ME[ME_event_counter]->sety(prim_vertex_y);
            JetTrackEvent_ME[ME_event_counter]->setz(prim_vertex_z);
            JetTrackEvent_ME[ME_event_counter]->setid(RunId);
            JetTrackEvent_ME[ME_event_counter]->setmult(refMult);
            JetTrackEvent_ME[ME_event_counter]->setn_prim(n_prim);
            JetTrackEvent_ME[ME_event_counter]->setn_non_prim(n_non_prim);
            JetTrackEvent_ME[ME_event_counter]->setn_tof_prim(n_tofmatch_prim);
            JetTrackEvent_ME[ME_event_counter]->setSE_ME_flag(SE_ME_flag);
            JetTrackEvent_ME[ME_event_counter]->setZDCx(ZDCx);
            JetTrackEvent_ME[ME_event_counter]->setBBCx(BBCx);
            JetTrackEvent_ME[ME_event_counter]->setvzVpd(vzVPD);
            JetTrackEvent_ME[ME_event_counter]->setQvecEtaPos(QvecEtaPos);
            JetTrackEvent_ME[ME_event_counter]->setQvecEtaNeg(QvecEtaNeg);
            JetTrackEvent_ME[ME_event_counter]->setcent9(cent9);

            Double_t Psi2 = 0.0;

            for(Int_t i_track = 0; i_track < N_Particles; i_track++)
            {
                // Particle information
                JetTrackParticle            = JetTrackEvent->getParticle(i_track);
                Float_t dca                 = JetTrackParticle->get_dca_to_prim();
                Float_t m2                  = JetTrackParticle->get_Particle_m2 ();
                Float_t nSPi                = JetTrackParticle->get_Particle_nSigmaPi();
                Float_t nSK                 = JetTrackParticle->get_Particle_nSigmaK();
                Float_t nSP                 = JetTrackParticle->get_Particle_nSigmaP();
                Float_t qp                  = JetTrackParticle->get_Particle_qp();
                Float_t nhitsfit            = JetTrackParticle->get_Particle_hits_fit();
                TLorentzVector TLV_Particle_prim = JetTrackParticle->get_TLV_Particle_prim();
                TLorentzVector TLV_Particle_glob = JetTrackParticle->get_TLV_Particle_glob();

                TLorentzVector TLV_Particle_use = TLV_Particle_prim;
                if(eflab_prim_glob == 1) TLV_Particle_use = TLV_Particle_glob;

                JetTrackParticle_ME = JetTrackEvent_ME[ME_event_counter]->createParticle();
                JetTrackParticle_ME ->set_dca_to_prim(dca);
                JetTrackParticle_ME ->set_Particle_m2(m2);
                JetTrackParticle_ME ->set_Particle_nSigmaPi(nSPi);
                JetTrackParticle_ME ->set_Particle_nSigmaK(nSK);
                JetTrackParticle_ME ->set_Particle_nSigmaP(nSP);
                JetTrackParticle_ME ->set_Particle_qp(qp);
                JetTrackParticle_ME ->set_TLV_Particle_prim(TLV_Particle_prim);
                JetTrackParticle_ME ->set_TLV_Particle_glob(TLV_Particle_glob);
                JetTrackParticle_ME ->set_Particle_hits_fit(nhitsfit);
            }

            //cout << "ME_event_counter (A) = " << ME_event_counter << ", N_Particles = " << N_Particles << endl;

            ME_event_counter++;
            if(ME_event_counter == N_max_events) // event buffer is full, start to create mixed events
            {
                for(Int_t mix_loop = 0; mix_loop < N_max_events; mix_loop++)
                {
                    if (mix_loop != 0  &&  mix_loop % 10 == 0)
                        cout << "-" << flush;

                    // Sample track number and z-vertex value
                    Double_t N_tracks_sample_d, N_z_vertex_sample, N_Psi_sample;
                    h_tracks_vs_z_vertex_array[ez_bin][emult_bin] ->GetRandom2(N_z_vertex_sample,N_tracks_sample_d);
                    h_Psi_vs_z_vertex_array[ez_bin][ePsi_bin]     ->GetRandom2(N_z_vertex_sample,N_Psi_sample);
                    //h_tracks_vs_z_vertex_sample->GetRandom2(N_z_vertex_sample,N_tracks_sample_d);
                    Int_t N_tracks_sample = (Int_t)N_tracks_sample_d;
                    //cout << "global event = " << counter << ", mix_loop = " << mix_loop <<  ", N_z_vertex_sample = " << N_z_vertex_sample << ", N_tracks_sample = " << N_tracks_sample << endl;

                    if(N_tracks_sample > N_max_events)
                    {
                        cout << "Error: N_tracks_sample > N_max_events" << endl;
                        break;
                    }

                    Int_t z_bin = -1;
                    if(N_z_vertex_sample > vertex_z_start_stop_delta[eBeamTimeNum][0] && N_z_vertex_sample < vertex_z_start_stop_delta[eBeamTimeNum][1])
                    {
                        z_bin = (Int_t)((N_z_vertex_sample-vertex_z_start_stop_delta[eBeamTimeNum][0])/vertex_z_start_stop_delta[eBeamTimeNum][2]);
                    }
                    Int_t mult_bin = -1;
                    if(N_tracks_sample_d > mult_start_stop_delta[eBeamTimeNum][0] && N_tracks_sample_d < mult_start_stop_delta[eBeamTimeNum][1])
                    {
                        mult_bin = (Int_t)((N_tracks_sample_d-mult_start_stop_delta[eBeamTimeNum][0])/mult_start_stop_delta[eBeamTimeNum][2]);
                    }
                    Int_t Psi_bin = -1;
                    if(Psi2 > Psi_start_stop_delta[eBeamTimeNum][0] && N_Psi_sample < Psi_start_stop_delta[eBeamTimeNum][1])
                    {
                        Psi_bin = (Int_t)((N_Psi_sample-Psi_start_stop_delta[eBeamTimeNum][0])/Psi_start_stop_delta[eBeamTimeNum][2]);
                    }

                    //cout << "z_bin = " << z_bin << ", mult_bin = " << mult_bin << ", Psi_bin = " << Psi_bin << ", N_tracks_sample_d = " << N_tracks_sample_d << endl;
                    //if(z_bin == ez_bin && mult_bin == emult_bin && Psi_bin == ePsi_bin)
                    {
                        F_mixed_event[ez_bin][emult_bin][ePsi_bin] ->cd();

                        Int_t Fill_flag = 1; // if all events in the loop have some entries left over then fill the event later
                        Double_t Et_total = 0.0;
                        for(Int_t ME_loop_counter = 0; ME_loop_counter < N_tracks_sample; ME_loop_counter++)
                        {
                            //input_SE_ME[i_SE_ME]->GetEntry( counter - N_max_events + ME_loop_counter + 1 );
                            //JetTrackEvent_ME = *JetTrackEvent;

                            //-----------------------------------------------------------------------------
                            // Event information
                            Float_t prim_vertex_x   = JetTrackEvent_ME[ME_loop_counter]->getx();
                            Float_t prim_vertex_y   = JetTrackEvent_ME[ME_loop_counter]->gety();
                            Float_t prim_vertex_z   = JetTrackEvent_ME[ME_loop_counter]->getz();
                            Int_t   RunId           = JetTrackEvent_ME[ME_loop_counter]->getid();
                            Float_t refMult         = JetTrackEvent_ME[ME_loop_counter]->getmult();
                            Float_t n_prim          = JetTrackEvent_ME[ME_loop_counter]->getn_prim();
                            Float_t n_non_prim      = JetTrackEvent_ME[ME_loop_counter]->getn_non_prim();
                            Int_t   n_tofmatch_prim = JetTrackEvent_ME[ME_loop_counter]->getn_tof_prim();
                            Int_t   SE_ME_flag      = JetTrackEvent_ME[ME_loop_counter]->getSE_ME_flag();
                            Float_t ZDCx            = JetTrackEvent_ME[ME_loop_counter]->getZDCx();
                            Float_t BBCx            = JetTrackEvent_ME[ME_loop_counter]->getBBCx();
                            Float_t vzVPD           = JetTrackEvent_ME[ME_loop_counter]->getvzVpd();
                            Int_t   N_Particles     = JetTrackEvent_ME[ME_loop_counter]->getNumParticle();
                            TVector2 QvecEtaPos     = JetTrackEvent_ME[ME_loop_counter]->getQvecEtaPos();
                            TVector2 QvecEtaNeg     = JetTrackEvent_ME[ME_loop_counter]->getQvecEtaNeg();
                            Int_t    cent9          = JetTrackEvent_ME[ME_loop_counter]->getcent9();

                            if((N_Particles-mix_loop) < 0) // at least one event in the loop is out of tracks, stop and don't fill the tree
                            {
                                Fill_flag = 0;
                                break;
                            }

                            if(ME_loop_counter == 0)
                            {
                                // Fill event information for jet
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].clearParticleList();
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setx(prim_vertex_x);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].sety(prim_vertex_y);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setz(prim_vertex_z);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setid(RunId);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setmult(refMult);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setn_prim(n_prim);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setn_non_prim(n_non_prim);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setn_tof_prim(n_tofmatch_prim);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setSE_ME_flag(SE_ME_flag);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setZDCx(ZDCx);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setBBCx(BBCx);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setvzVpd(vzVPD);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setQvecEtaPos(QvecEtaPos);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setQvecEtaNeg(QvecEtaNeg);
                                JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].setcent9(cent9);

                                h_tracks_vs_z_vertex[ez_bin][emult_bin][ePsi_bin] ->Fill(N_z_vertex_sample,N_tracks_sample_d);
                                h_Psi2[ez_bin][emult_bin][ePsi_bin]               ->Fill(Psi2);
                            }
                            //-----------------------------------------------------------------------------

                            // select random track from this event
                            //cout << "ME_loop_counter = " << ME_loop_counter << ", N_Particles = " << N_Particles << ", mix_loop = " << mix_loop << endl;
                            Int_t track_num_random    = (Int_t)ran.Integer(N_Particles-mix_loop);
                            //if(ME_loop_counter == 0) cout << "track_num_random = " << track_num_random << ", array size = " << (Int_t)ME_track_number_vector[ME_loop_counter].size() << endl;
                            Int_t track_id_num_random = ME_track_number_vector[ME_loop_counter][track_num_random];
                            ME_track_number_vector[ME_loop_counter][track_num_random] = ME_track_number_vector[ME_loop_counter][N_Particles-mix_loop-1];

                            //cout << "Test C, N_Particles = " << N_Particles << ", track_id_num_random = " << track_id_num_random << endl;

                            //JetTrackEvent_ME.resize(N_max_events);

                            //mult_start_stop_delta[eBeamTimeNum][1]
                            //file_entries_SE_ME[0]

                            // Particle information
                            JetTrackParticle_ME         = JetTrackEvent_ME[ME_loop_counter]->getParticle(track_id_num_random);
                            Float_t dca                 = JetTrackParticle_ME->get_dca_to_prim();
                            Float_t m2                  = JetTrackParticle_ME->get_Particle_m2 ();
                            Float_t nSPi                = JetTrackParticle_ME->get_Particle_nSigmaPi();
                            Float_t nSK                 = JetTrackParticle_ME->get_Particle_nSigmaK();
                            Float_t nSP                 = JetTrackParticle_ME->get_Particle_nSigmaP();
                            Float_t qp                  = JetTrackParticle_ME->get_Particle_qp();
                            Float_t nhitsfit            = JetTrackParticle_ME->get_Particle_hits_fit();
                            TLorentzVector TLV_Particle_prim = JetTrackParticle_ME->get_TLV_Particle_prim();
                            TLorentzVector TLV_Particle_glob = JetTrackParticle_ME->get_TLV_Particle_glob();

                            TLorentzVector TLV_Particle_use = TLV_Particle_prim;
                            if(eflab_prim_glob == 1) TLV_Particle_use = TLV_Particle_glob;

                            Et_total += TLV_Particle_use.Et();

                            JetTrackParticle_Fill = JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin].createParticle();
                            JetTrackParticle_Fill ->set_dca_to_prim(dca);
                            JetTrackParticle_Fill ->set_Particle_m2(m2);
                            JetTrackParticle_Fill ->set_Particle_nSigmaPi(nSPi);
                            JetTrackParticle_Fill ->set_Particle_nSigmaK(nSK);
                            JetTrackParticle_Fill ->set_Particle_nSigmaP(nSP);
                            JetTrackParticle_Fill ->set_Particle_qp(qp);
                            JetTrackParticle_Fill ->set_TLV_Particle_prim(TLV_Particle_prim);
                            JetTrackParticle_Fill ->set_TLV_Particle_glob(TLV_Particle_glob);
                            JetTrackParticle_Fill ->set_Particle_hits_fit(nhitsfit);

                            Double_t track_pT  = TLV_Particle_use.Pt();
                            if(track_pT != track_pT) continue; // that is a NaN test. It always fails if track_pT = nan.

                            Float_t Phi = TLV_Particle_use.Phi();
                            Float_t Eta = TLV_Particle_use.PseudoRapidity();
                            h_Phi_vs_eta[ez_bin][emult_bin][ePsi_bin] ->Fill(Eta,Phi);

                            if(0)
                            {
                                cout << "mix_loop = " << mix_loop << ", ME_loop_counter = " << ME_loop_counter << ", track_id_num_random = " << track_id_num_random
                                    <<  ", track_num_random = " << track_num_random << ", Phi = " << Phi << ", Eta = " << Eta << ", refMult = " << refMult << endl;
                            }
                        }
                        h_Et[ez_bin][emult_bin][ePsi_bin]                 ->Fill(Et_total);

                        if(Fill_flag == 1)
                        {
                            Tree_JetTrackEvent_Fill[ez_bin][emult_bin][ePsi_bin]  ->Fill();
                        }
                    }
                }
                ME_event_counter = 0;
            }

        }
    }
    return 1;
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
Int_t StJetAnalysis::Finish()
{
    cout << "Finish started" << endl;

    Outputfile->cd();
    Outputfile->mkdir("Test");
    Outputfile->cd("Test");
    for(Int_t i = 0; i < 4; i++)
    {
        h_delta_phi_trigger_recoil[i] ->Write();
        h_eta_recoil[i]               ->Write();
        h_phi_recoil[i]               ->Write();
        h_eta_trigger[i]              ->Write();
        h_phi_trigger[i]              ->Write();
    }
    Outputfile->cd();
    h_Delta_phi_trig_jet ->Write();

    if(eMode == 0)
    {
        Outputfile->cd();

        p_parameters  ->Write();
        p_jet_area_values ->Write();

        //h_tracks_vs_z_vertex ->Write();
        //h_Psi2               ->Write();

        if(eMode == 0 || eMode == 1)
        {
            for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
            {
                //h_Phi_vs_eta[i_z]            ->Write();
                h_Phi_vs_eta_random_phi[i_z] ->Write();
            }
        }

        Outputfile->Close();
    }

    if(eMode == 42)
    {

        p_parameters                           ->Write();
        p_jet_area_values                      ->Write();
        p_v2_vs_pt_jet                         ->Write();
        p_v2_vs_pt                             ->Write();
        h_Psi_Full                             ->Write();
        h_Psi_etapos                           ->Write();
        h_Psi_etaneg                           ->Write();
        h_phi                                  ->Write();

        Outputfile->cd();
        Outputfile->mkdir("Delta_pt");
        Outputfile->cd("Delta_pt");

        for(Int_t i_orig_smear = 0; i_orig_smear < N_orig_smear; i_orig_smear++)
        {
            for(Int_t i = 0; i < N_jet_areas; i++)
            {
                for(Int_t i_Psi = ePsi_bin; i_Psi < ePsi_bin+1; i_Psi++)
                {
                    h_Delta_pt_vs_embed_pt[i][i_orig_smear][i_Psi]        ->Write();
                    h_Delta_pt_vs_embed_pt_weight[i][i_orig_smear][i_Psi] ->Write();
                }
            }
        }


        Outputfile->cd();

        Outputfile->Close();
    }

    if(eMode == 310)
    {
        for(Int_t i_hist = 0; i_hist < N_m2_radial_dist; i_hist++)
        {
            for(Int_t i_charge = 0; i_charge < 2; i_charge++)
            {
                h2D_m2_vs_pT[i_hist][i_charge] ->Write();
            }
        }
    }


    if(eMode == 3 || eMode == 31 || eMode == 32 || eMode == 312 || eMode == 311)
    {
        Outputfile->cd();

        for(Int_t i_charge = 0; i_charge < 2; i_charge++)
        {
            HistName = "p_pt_Ach_";
            HistName += i_charge;
            p_pt_Ach[i_charge] -> Write();
        }
        h_Ach -> Write();



        if(eMode == 312)
        {
            h_PYTHIA_hard_bin_weigh_factors        ->Write();
            h_PYTHIA_hard_bin_high_pT_N_events     ->Write();

            Outputfile->cd();
            Outputfile->mkdir("Delta_pt");
            Outputfile->cd("Delta_pt");

            for(Int_t i_orig_smear = 0; i_orig_smear < N_orig_smear; i_orig_smear++)
            {
                for(Int_t i = 0; i < N_jet_areas; i++)
                {
                    for(Int_t i_Psi = ePsi_bin; i_Psi < ePsi_bin+1; i_Psi++)
                    {
                        h_Delta_pt_vs_embed_pt[i][i_orig_smear][i_Psi] ->Write();
                    }
                }
            }


            Outputfile->cd();
        }
        NT_ReCoil_Jet                          ->Write();
        h_area_vs_jet_pt                       ->Write();
        h_area                                 ->Write();
        h_track_pT                             ->Write();
        h_track_pT_piplus                      ->Write();
        h_track_pT_cut                         ->Write();
        h_tracks_above_threshold_per_event     ->Write();

        for(Int_t i_hist = 0; i_hist < N2D_tracks_above_threshold; i_hist++)
        {
            h2D_tracks_above_threshold[i_hist] ->Write();
        }
        //c_3D            ->Write();
        p_parameters                           ->Write();
        p_jet_area_values                      ->Write();
        h_ratio_sub_lead_to_lead_pt_vs_lead_pt ->Write();
        h_sub_lead_vs_lead_pt                  ->Write();
        h_jet_rho_vs_Et                        ->Write();

        if(eMode == 311 || eMode == 312)
        {
            h_1D_det_diff_part_lvl_area        ->Write();
            h_1D_det_diff_part_lvl_area_lt5GeV ->Write();
            h_2D_det_vs_part_lvl_area          ->Write();
            h_2D_det_vs_part_lvl_area_lt5GeV   ->Write();

            for(Int_t i = 0; i < N_jet_areas; i++)
            {
                h_Sim_original_pT_cut[i]                      ->Write();
                h_Sim_original_pT_all[i]                      ->Write();
                h_Sim_smeared_pT_not_matched[i]               ->Write();
                h2D_Sim_matched_pT_vs_original_pT[i]          ->Write();
                h2D_Sim_original_pT_vs_matched_pT[i]          ->Write();
                h_matched_tracks_fraction[i]                  ->Write();
                h2D_matched_tracks_fraction_vs_original_pT[i] ->Write();
            }
        }

        for(Int_t i = 0; i < 2; i++)
        {
            p_Array_leading_pt_bins[i] ->Write();
        }

        for(Int_t i_orig_smear = 0; i_orig_smear < 4; i_orig_smear++)
        {
            h_trigger_track[i_orig_smear]               ->Write();
            h_trigger_track_vs_global_bin[i_orig_smear] ->Write();
        }

        Outputfile->cd();
        Outputfile->mkdir("Jet_vs_area");
        Outputfile->cd("Jet_vs_area");

        for(Int_t ipt_lead = 0; ipt_lead < (N_leading_pt_bins+1); ipt_lead++)
        {
            h_area_vs_recoil_jet_pt[ipt_lead] ->Write();
            h_rho[ipt_lead]                   ->Write();
        }

        Outputfile->cd();
        Outputfile->mkdir("Track_distributions");
        Outputfile->cd("Track_distributions");

        for(Int_t i_orig_smear = 0; i_orig_smear < N_orig_smear; i_orig_smear++)
        {
            for(Int_t ipt_lead = 0; ipt_lead < (N_leading_pt_bins+1); ipt_lead++)
            {
                h_N_tracks_dijet[ipt_lead][i_orig_smear]       ->Write();
                h2D_mult_vs_global_bin[ipt_lead][i_orig_smear] ->Write();
            }
        }

        Outputfile->cd();
        Outputfile->mkdir("Jet_QA");
        Outputfile->cd("Jet_QA");

        for(Int_t i_orig_smear = 0; i_orig_smear < N_orig_smear; i_orig_smear++)
        {
            for(Int_t i = 0; i < N_jet_areas; i++)
            {
                h_jet_area[i][i_orig_smear]        ->Write();
                h_jet_rho[i][i_orig_smear]         ->Write();
                h_jet_per_event[i][i_orig_smear]   ->Write();
                h_dijet_per_event[i][i_orig_smear] ->Write();
            }
        }

        Outputfile->cd();
        h2D_dijet_pt_sub_test->Write();
        h2D_dijet_pt_sub_weight_test->Write();
        Outputfile->mkdir("Jet_spectra");
        Outputfile->cd("Jet_spectra");

        for(Int_t i_orig_smear = 0; i_orig_smear < 4; i_orig_smear++)
        {
            for(Int_t i = 0; i < N_jet_areas; i++)
            {
                for(Int_t ipt_lead = 0; ipt_lead < (N_leading_pt_bins+1); ipt_lead++)
                {
                    h_jet_pt[i][ipt_lead][i_orig_smear]         ->Write();
                    h_jet_pt_sub[i][ipt_lead][i_orig_smear]     ->Write();

                    h2D_dijet_pt[i][ipt_lead][i_orig_smear]     ->Write();
                    h2D_dijet_pt_sub[i][ipt_lead][i_orig_smear] ->Write();
                }
            }
        }

        Outputfile->cd();
        Outputfile->mkdir("QA_hist_arrays");
        Outputfile->cd("QA_hist_arrays");

        for(Int_t i_z = ez_bin; i_z < ez_bin+1; i_z++)
        {
            for(Int_t i_mult = emult_bin; i_mult < emult_bin+1; i_mult++)
            {
                for(Int_t i_Psi = ePsi_bin; i_Psi < ePsi_bin+1; i_Psi++)
                {
                    //h_jet_area_array[i_z][i_mult][i_Psi]           ->Write();
                    //h_jet_rhoarea_array[i_z][i_mult][i_Psi]        ->Write();
                    //h_jet_rho_array[i_z][i_mult][i_Psi]            ->Write();
                    //h_jet_Et_array[i_z][i_mult][i_Psi]             ->Write();
                    //h_jet_Et_array_weight[i_z][i_mult][i_Psi]      ->Write();
                    h_jet_per_event_array[i_z][i_mult][i_Psi]      ->Write();
                    //h_tracks_per_event_array[i_z][i_mult][i_Psi]   ->Write();
                    h2D_jet_rho_vs_mult_array[i_z][i_mult][i_Psi]  ->Write();
                }
            }
        }

        Outputfile->cd();
        Outputfile->mkdir("Acc_Recoil");
        Outputfile->cd("Acc_Recoil");

        for(Int_t i_area_acc = 0; i_area_acc < N_jet_areas; i_area_acc++)
        {
            for(Int_t ipt_lead = 0; ipt_lead < (N_leading_pt_bins+1); ipt_lead++)
            {
                h_N_accepted_recoil_jets[i_area_acc][ipt_lead] ->Write();
            }
        }

        Outputfile->cd();

        Outputfile->Close();
    }


    if(eMode == 1 || eMode == 2 || eMode == 4)
    {


        Int_t start_z    = 0;
        Int_t start_mult = 0;
        Int_t start_Psi  = 0;
        //Int_t start_Psi  = ePsi_bin;
        Int_t stop_z     = N_z_vertex_bins;
        Int_t stop_mult  = N_mult_bins;
        Int_t stop_Psi   = N_Psi_bins;
        //Int_t stop_Psi   = ePsi_bin+1;
        if(eMode == 2
           || eMode == 1
          )
        {
            start_z    = ez_bin;
            start_mult = emult_bin;
            start_Psi  = ePsi_bin;
            stop_z     = ez_bin+1;
            stop_mult  = emult_bin+1;
            stop_Psi   = ePsi_bin+1;
        }

        for(Int_t i_z = start_z; i_z < stop_z; i_z++)
        {
            for(Int_t i_mult = start_mult; i_mult < stop_mult; i_mult++)
            {
                for(Int_t i_Psi = start_Psi; i_Psi < stop_Psi; i_Psi++)
                {
                    F_mixed_event[i_z][i_mult][i_Psi] ->cd();
                    p_parameters     ->Write();
                    p_jet_area_values ->Write();

                    h_PsiA_vs_PsiB  ->Write();
                    h_PsiA          ->Write();
                    h_PsiB          ->Write();
                    h_Psi_Full      ->Write();

                    h_tracks_vs_z_vertex[i_z][i_mult][i_Psi] ->Write();
                    h_Psi2[i_z][i_mult][i_Psi]               ->Write();
                    h_Et[i_z][i_mult][i_Psi]                 ->Write();
                    h_Momentum[i_z][i_mult][i_Psi]           ->Write();
                    h_Phi_vs_eta[i_z][i_mult][i_Psi]         ->Write();

                    Tree_JetTrackEvent_Fill[i_z][i_mult][i_Psi] ->Write("",TObject::kOverwrite);
                    if(eMode == 2) Inputfile                    ->Close();
                    F_mixed_event[i_z][i_mult][i_Psi]           ->Close();
                }
            }
        }
    }

    if(eMode == 11)
    {
        cout << "Save data to output file" << endl;
        Outputfile->cd();
        h_pT_eta_phi                           ->Write();
        p_parameters                           ->Write();
        p_jet_area_values                      ->Write();
        h_tracks_above_threshold_per_event     ->Write();
        for(Int_t i_hist = 0; i_hist < N2D_tracks_above_threshold; i_hist++)
        {
            h2D_tracks_above_threshold[i_hist] ->Write();
        }

        for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
        {
            for(Int_t i_mult = 0; i_mult < N_mult_bins; i_mult++)
            {
                h_tracks_vs_z_vertex_array[i_z][i_mult] ->Write();

                for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
                {
                    for(Int_t i_pT = 0; i_pT < N_track_pt_bins_eta_phi; i_pT++)
                    {
                        h2D_track_eta_vs_phi[i_z][i_mult][i_Psi][i_pT] ->Write();
                    }
                }
            }
        }
        for(Int_t i_z = 0; i_z < N_z_vertex_bins; i_z++)
        {
            for(Int_t i_Psi = 0; i_Psi < N_Psi_bins; i_Psi++)
            {
                h_Psi_vs_z_vertex_array[i_z][i_Psi] ->Write();
            }
        }

        Outputfile->Close();
    }

    return 1;
}
//------------------------------------------------------------------------------------------------------------------



