#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

//------------------------
#include "AliHelix.h"
#include "TLorentzVector.h"
#include "TSystem.h"
//------------------------

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliKalmanTrack.h"

#include "AliTRDpadPlane.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliESDfriend.h"

#include "AliTRDarrayADC.h"

#include "AliPIDResponse.h"
#include "AliPID.h"

#include "AliESDtrackCuts.h"

#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliESDRun.h"

#include "AliMultSelection.h"

#include "AliCDBEntry.h"
#include "TClonesArray.h"
#include "TGeoMatrix.h"
#include "AliAlignObjParams.h"

#include "AliRunTag.h"
#include "TObjString.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliTRDCalPad.h"
#include "AliTRDCalDet.h"
#include "AliTRDCalOnlineGainTable.h"
#include "AliTRDCalROC.h"
#include "TPolyMarker.h"

#include "AliTRDCommonParam.h"

#include "AliESDTrdTracklet.h"
#include "TProfile.h"
#include "TGraph.h"

//------------------------
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
//------------------------

#include "Ali_AS_analysis_TRD_digits.h"


#include <iostream>
#include <iomanip>
using namespace std;

static Int_t flag_plot_event = 0;
static TString HistName;

static TFile* dfile;
static TFile* TRD_alignment_file;
static TFile* TRD_calibration_file_AA;
static TGraph* tg_v_fit_vs_det;
static TGraph* tg_LA_factor_fit_vs_det;
static TH1D* h_v_fit_vs_det;
static TH1D* h_LA_factor_fit_vs_det;
static TGeoHMatrix TM_TRD_rotation_det[540];
static TGeoHMatrix TM_TRD_rotation_sector[18];
static TVector3    TV3_TRD_translation[540];
static const Double_t TRD_lorentz_angle = TMath::DegToRad()*8.8;
static const char *pathdatabase="alien://folder=/alice/data/2016/OCDB"; // for pPb
//static const char *pathdatabase="alien://folder=/alice/data/2015/OCDB"; // for PbPb
static AliTRDCalDet *ChamberVdrift;
static AliTRDCalDet *ChamberT0;
static AliTRDCalDet *ChamberExB;
static AliTRDCalDet *chambergain;
static AliTRDCalOnlineGainTable *KryptoGain;
static AliTRDCalPad *PadNoise;
static AliTRDCalPad *LocalT0_pad;
static AliTRDCalROC* CalROC;
static AliTRDCalROC  *LocalT0;            //  Pad wise T0 calibration object
static const Int_t N_pT_bins = 5;
static const Double_t pT_ranges[N_pT_bins+1] = {0.2,0.5,1.0,2.0,3.0,5.0};
static const Double_t TRD_Impact_distance_in_drift = 3.35;

static AliTRDCommonParam* fParam;
static AliESDfriend *esdFr = NULL;
static AliESDInputHandler *esdH = NULL;
static TString esdFriendTreeFName;
//static Class_peak_finder my_class_peak_finder;

ClassImp(Ali_AS_analysis_TRD_digits)
ClassImp(StJetTrackEvent)
ClassImp(StJetTrackParticle)
ClassImp(StEMCal)

    //________________________________________________________________________
    Ali_AS_analysis_TRD_digits::Ali_AS_analysis_TRD_digits(const char *name)
    : AliAnalysisTaskSE(name),
    JetTrackEvent(0),JetTrackParticle(0),JetEMCal(0),Tree_AS_Event(0), fEventNoInFile(-2), N_good_events(0), fDigitsLoadedFlag(kFALSE),
    fListOfHistos(0x0),fTree(0x0),fCorrTaskSetting(""),fInputEvent(NULL),h_dca(0x0),h_dca_xyz(0x0), h2D_TPC_dEdx_vs_momentum(0x0), h_ADC_tracklet(0x0), h_ADC_vs_time(0x0), fPIDResponse(0), EsdTrackCuts(0)
{
    // Constructor

    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());

    // Output slot #0 id reserved by the base class for AOD
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());

}



//_______________________________________________________________________
Bool_t Ali_AS_analysis_TRD_digits::UserNotify()
{
    cout << "" << endl;
    cout << "In UserNotify" << endl;


    cout << "All pointers deleted" << endl;

    //AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>
    //    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if ( ! esdH ) return kFALSE;
    if ( ! esdH->GetTree() ) return kFALSE;
    if ( ! esdH->GetTree()->GetCurrentFile() ) return kFALSE;


    //-----------------------------------
    TList* list = esdH->GetUserInfo();
    //list->Print();
    //list->Dump();

    //cout << "List: " << endl;
    //list->ls();
    TList* list_cdblist = (TList*)list->FindObject("cdbList");  // list contains paths of calibration data
    Int_t list_size = list_cdblist->Capacity();
    //cout << "cdblist, with size: " << list_size <<  endl;
    //list_cdblist->ls();
    TObjString* obj_string = (TObjString*)list_cdblist->At(0);
    TString string = obj_string->GetString();
    //cout << "String: " << string.Data() << endl;
    Int_t index_A = string.Index("[");
    string.Remove(0,index_A+1);
    //cout << "String: " << string.Data() << endl;
    Int_t index_B = string.Index(",");
    string.Remove(index_B,string.Sizeof());
    Int_t run_number_from_list = string.Atoi();
    cout << "String: " << string.Data() << ", run_number_from_list: " << run_number_from_list << endl;
    //-----------------------------------


    // AliCDBEntry->GetObject()->IsA()->GetName()

    TString fname = esdH->GetTree()->GetCurrentFile()->GetName();
    TString Tree_name = esdH->GetTree()->GetName();
    FileStat_t file_stat;
    Int_t PathInfo = gSystem->GetPathInfo(fname.Data(),file_stat);
    cout << "PathInfo: " << PathInfo << ", fname: " << fname << endl;
    //TFile* file = TFile::Open(fname.Data());
    //cout << "Zombie: " << file->IsZombie() << ", header size: " << file->Sizeof() << ", FileBytesRead: " << file->GetFileBytesRead() << endl;

    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!inputHandler)
    {
	printf("WARNING: Inputhandler not available \n");
    }
    else
    {
	printf("Inputhandler available \n");

	fPIDResponse = inputHandler->GetPIDResponse();

        cout << "Got PID response" << endl;
    }

    fEventNoInFile = -1;
    N_good_events  = 0;

    EsdTrackCuts = new AliESDtrackCuts();

    EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
    EsdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.52);
    EsdTrackCuts->AliESDtrackCuts::SetMinNClustersTPC(50); // 60, Automatically requires TPC refitted tracks?
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXY(10.0);
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(10.0);
    EsdTrackCuts->AliESDtrackCuts::SetPtRange(0.15,200.0); // 0.15, 200.0
    EsdTrackCuts->AliESDtrackCuts::SetEtaRange(-1.0,1.0); // 0.85


    if(fname.Contains("/home/"))
    {
        TRD_alignment_file = TFile::Open("/home/ceres/schmah/ALICE/Database/TRD_Align_2016.root");
	cout << "Local alignment file loaded" << endl;
    }
    else
    {
	TRD_alignment_file = TFile::Open("alien:///alice/data/2016/OCDB/TRD/Align/Data/Run0_999999999_v1_s0.root");
        cout << "Alignment file from database loaded" << endl;
    }

    cout << "End of UserNotify" << endl;
    return kTRUE;
}


//________________________________________________________________________
void Ali_AS_analysis_TRD_digits::UserCreateOutputObjects()
{
    cout << "" << endl;
    cout << "In UserCreateOutputObjects" << endl;


    OpenFile(1);
    cout << "File opened" << endl;

    fListOfHistos = new TList();
    fListOfHistos ->SetOwner();


    OpenFile(2);
    cout << "File opened" << endl;

   

    JetTrackEvent       = new StJetTrackEvent();
    JetTrackParticle    = new StJetTrackParticle();
    JetEMCal            = new StEMCal();
    Tree_AS_Event  = NULL;
    Tree_AS_Event  = new TTree("JetTrackEvent" , "JetTrackEvent" );
    Tree_AS_Event  ->Branch("Events"  , "StJetTrackEvent",  JetTrackEvent);

    PostData(1,fListOfHistos);
    PostData(2,Tree_AS_Event);

    cout << "PostData called" << endl;
}



//________________________________________________________________________
Bool_t Ali_AS_analysis_TRD_digits::NextEvent(Bool_t preload)
{
    fEventNoInFile++;
    //cout << "fEventNoInFile: " << fEventNoInFile << endl;

    return kTRUE;
}

//________________________________________________________________________
void Ali_AS_analysis_TRD_digits::UserExec(Option_t *)
{
    //cout << "" << endl;
    //cout << "Analysis started" << endl;
    //cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;

    Int_t flag_calibrated = 0; // 0 = standard fixed precalibration used, 1 = use pre calibration from root input file
    //-----------------------------------------------------------------
    // IMPORTANT: call NextEvent() for book-keeping
    NextEvent();
    fInputEvent           = InputEvent();
    //if(fEventNoInFile > 50) return;
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------
    // prepare event data structures
    AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD)
    {
	printf("ERROR: fESD not available\n");
	return;
    }



    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if(man)
    {
        //Int_t run_id = man->GetRunFromPath(); // doesn't work
	//cout << "Got AliAnalysisManager, run_id: " << run_id << endl;
	AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
	if(inputHandler)
	{
	    //cout << "Got AliInputEventHandler" << endl;
	    fPIDResponse = inputHandler->GetPIDResponse();
	}
    }
    //cout << "cent: " << fPIDResponse->GetCurrentCentrality() << endl;

    Int_t 	   eventNumber      = fESD ->GetEventNumberInFile();
    Int_t          N_tracks         = fESD ->GetNumberOfTracks();
    Int_t          N_TRD_tracks     = fESD ->GetNumberOfTrdTracks();
    Int_t          N_TRD_tracklets  = fESD ->GetNumberOfTrdTracklets(); // online
    Float_t        magF             = fESD ->GetMagneticField();
    const AliESDVertex* PrimVertex  = fESD ->GetPrimaryVertex();
    Int_t          RunNum           = fESD ->GetRunNumber();
    Double_t       T0zVertex        = fESD ->GetT0zVertex();
    AliCentrality* Centrality       = fESD ->GetCentrality();
    Double_t       MeanBeamIntAA    = fESD ->GetESDRun()->GetMeanIntensity(0,0);

    //printf("RunNum: %d \n",RunNum);

    Double_t Sign_magnetic_field = (magF/fabs(magF));
    //cout << "Trigger: " <<  fESD->GetFiredTriggerClasses() << endl;


    if(N_tracks == 0)
    {
        //printf("   ---------> Event number: %d, N_tracks: %d, N_TRD_tracklets: %d \n",fEventNoInFile,N_tracks,N_TRD_tracklets);

	// Skip empty event
	return;
    }

    



    //-----------------------------------------------------------------
    TClonesArray* TC_tofHits = fESD->GetESDTOFHits();
    //printf(" \n");
    //printf("   ------------> TC_tofHits: \n");
    // AliCDBEntry->GetObject()->IsA()->GetName()   only root 5
    Int_t N_entries_TOF = TC_tofHits ->GetEntries();
    Double_t TOF_R = ((AliESDTOFHit*)TC_tofHits ->First())->GetR();
    //cout << "N_tracks: " << N_tracks <<  ", TC_tofHits: " << TC_tofHits << ", IsA: " << TC_tofHits ->First()->IsA()->GetName() << ", TOF_R: " << TOF_R << ", N_entries_TOF: " << N_entries_TOF << endl;
    for(Int_t i_tof = 0; i_tof < N_entries_TOF; i_tof++)
    {
        Double_t TOF_Z       = ((AliESDTOFHit*)TC_tofHits->At(i_tof))->GetZ();
        Int_t    TOF_channel = ((AliESDTOFHit*)TC_tofHits->At(i_tof))->GetTOFchannel();
        Double_t TOF_time    = ((AliESDTOFHit*)TC_tofHits->At(i_tof))->GetTime();
        //printf("i_tof: %d, TOF_Z: %4.3f, TOF_channel: %d, TOF_time: %d \n",i_tof,TOF_Z,TOF_channel,TOF_time);
    }
    //printf(" \n");
    //-----------------------------------------------------------------



    Int_t SE_ME_Flag = 0;

    TVector2 TV2_Qvec(0.0,0.0);


    // Fill event information for jets
    JetTrackEvent->clearParticleList();
    JetTrackEvent->setx(PrimVertex->GetX());
    JetTrackEvent->sety(PrimVertex->GetY());
    JetTrackEvent->setz(PrimVertex->GetZ());
    JetTrackEvent->setid(RunNum);
    JetTrackEvent->setmult(N_tracks);
    JetTrackEvent->setn_prim(N_tracks);
    JetTrackEvent->setn_TRD_tracklets(N_TRD_tracklets);
    JetTrackEvent->setn_tof_hits(N_entries_TOF);
    JetTrackEvent->setSE_ME_flag(SE_ME_Flag);

    JetTrackEvent ->setBeamIntAA(MeanBeamIntAA);
    JetTrackEvent ->setT0zVertex(T0zVertex);

    AliMultSelection *MultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
    if(MultSelection)
    {
	// V0MEq, V0AEq, V0CEq, SPDTracklets

	JetTrackEvent ->setcent_class_ZNA(MultSelection->GetMultiplicityPercentile("ZNA"));
	JetTrackEvent ->setcent_class_ZNC(MultSelection->GetMultiplicityPercentile("ZNC"));
	JetTrackEvent ->setcent_class_V0A(MultSelection->GetMultiplicityPercentile("V0A"));
	JetTrackEvent ->setcent_class_V0C(MultSelection->GetMultiplicityPercentile("V0C"));
	JetTrackEvent ->setcent_class_V0M(MultSelection->GetMultiplicityPercentile("V0M"));
	JetTrackEvent ->setcent_class_CL0(MultSelection->GetMultiplicityPercentile("CL0"));
	JetTrackEvent ->setcent_class_CL1(MultSelection->GetMultiplicityPercentile("CL1"));
	JetTrackEvent ->setcent_class_SPD(MultSelection->GetMultiplicityPercentile("SPDTracklets"));
	JetTrackEvent ->setcent_class_V0MEq(MultSelection->GetMultiplicityPercentile("V0MEq"));
	JetTrackEvent ->setcent_class_V0AEq(MultSelection->GetMultiplicityPercentile("V0AEq"));
	JetTrackEvent ->setcent_class_V0CEq(MultSelection->GetMultiplicityPercentile("V0CEq"));
    }


    JetTrackEvent->setZDCx(0);
    JetTrackEvent->setBBCx(0);
    JetTrackEvent->setvzVpd(0);

    JetTrackEvent->setQvecEtaPos(TV2_Qvec);
    JetTrackEvent->setQvecEtaNeg(TV2_Qvec);
    JetTrackEvent->setcent9(0);

    //-----------------------------------------------------------------
#if 0
    cout << "" << endl;
    cout << "----------------------------------------------------------------------------------------" << endl;
    printf("Event number: %d, N_tracks: %d, N_TRD_tracklets: %d \n",fEventNoInFile,N_tracks,N_TRD_tracklets);
    printf("cent(ZNA): %4.3f \n",MultSelection->GetMultiplicityPercentile("ZNA"));
    printf("cent(ZNC): %4.3f \n",MultSelection->GetMultiplicityPercentile("ZNC"));
    printf("cent(V0A): %4.3f \n",MultSelection->GetMultiplicityPercentile("V0A"));
    printf("cent(V0C): %4.3f \n",MultSelection->GetMultiplicityPercentile("V0C"));
    printf("cent(V0M): %4.3f \n",MultSelection->GetMultiplicityPercentile("V0M"));
    printf("cent(CL0): %4.3f \n",MultSelection->GetMultiplicityPercentile("CL0"));
    printf("cent(CL1): %4.3f \n",MultSelection->GetMultiplicityPercentile("CL1"));
    printf("cent(SPDTracklets): %4.3f \n",MultSelection->GetMultiplicityPercentile("SPDTracklets"));
    printf("cent(V0MEq): %4.3f \n",MultSelection->GetMultiplicityPercentile("V0MEq"));
    printf("cent(V0AEq): %4.3f \n",MultSelection->GetMultiplicityPercentile("V0AEq"));
    printf("cent(V0CEq): %4.3f \n",MultSelection->GetMultiplicityPercentile("V0CEq"));
#endif
    //-----------------------------------------------------------------




    //-----------------------------------------------------------------------------------------------------
    // EMCal
    JetTrackEvent->clearEMCalList();

    //printf("Start of EMCal clusters \n");
    Int_t nclus;
    TClonesArray * arrClustersProcess = NULL;
    if(!fCorrTaskSetting.CompareTo(""))
    { // special version of correction task needed?
        nclus = fInputEvent->GetNumberOfCaloClusters();
    }
    else
    {
        arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
        if(!arrClustersProcess)
        {
            AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
        }
        nclus = arrClustersProcess->GetEntries();
    }
    //printf("End of getting EMCal array \n");

    // loop over clusters
    for(Long_t i_clus = 0; i_clus < nclus; i_clus++)
    {
        AliVCluster* clus = NULL;
        if(fInputEvent->IsA()==AliESDEvent::Class())
        {
            if(arrClustersProcess)
                clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(i_clus));
            else
                clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i_clus));
        }
        else
        {
            printf("WARNING: No EMCAl ESD information! \n");
        }
        if(!clus) continue;

        Double_t cluster_energy = clus ->E();
        Bool_t   isExotic       = clus ->GetIsExotic();
        Double_t Dx_to_track    = clus ->GetTrackDx();
        Double_t Dz_to_track    = clus ->GetTrackDz();
        Bool_t   IsEMCal        = clus ->IsEMCAL();
        Double_t EMCalChi2      = clus ->Chi2();
        UShort_t NCells         = (UShort_t)clus ->GetNCells();
        Float_t  clusPos[3]={0,0,0};
        clus->GetPosition(clusPos);
        TLorentzVector clusterVector;
        Double_t clusPrimVertex[3] = {PrimVertex->GetX(),PrimVertex->GetY(),PrimVertex->GetZ()};
        clus->GetMomentum(clusterVector,clusPrimVertex);

        JetEMCal = JetTrackEvent->createEMCal();
        JetEMCal ->set_TLV_EMCal(clusterVector);
        JetEMCal ->set_cluster_pos(clusPos[0],clusPos[1],clusPos[2]);
        JetEMCal ->set_Dx_to_track(Dx_to_track);
        JetEMCal ->set_Dz_to_track(Dz_to_track);
        JetEMCal ->set_N_cells_in_cluster(NCells);
        JetEMCal ->set_isExotic(isExotic);

        delete clus;
        //printf("i_clus: %d, energy: %4.3f, NCells: %d \n",i_clus,cluster_energy,(Int_t)NCells);
    }
    //printf("End of EMCal clusters \n");
    //-----------------------------------------------------------------------------------------------------





    Int_t N_good_tracks = 0;


    //-----------------------------------------------------------------
    // Track loop
    //cout << "" << endl;
    //cout << "-----------------------------------------------------------------" << endl;
    //cout << "Start matching " << N_tracks << " TPC tracks with " << TV3_TRD_hits_middle.size() << " TRD pads" << endl;
    N_good_tracks = 0;
    Int_t N_matched_TRD_hits_total = 0;
    for(Int_t iTracks = 0; iTracks < N_tracks; iTracks++)
    {
	//---------------------------------------------------------------
	// Gather track information

	// We always want the ESD track
	AliESDtrack* track = fESD->GetTrack(iTracks);
	if(!track)
	{
	    printf("ERROR: Could not receive track %d\n", iTracks);
	    continue;
	}


        if(!EsdTrackCuts->AcceptTrack(track)) continue;


        Double_t TRD_signal   = track ->GetTRDsignal(); // truncated mean signal?
        Double_t Track_pT     = track ->Pt();
        Double_t Track_p      = track ->P();
        Double_t p_vec[3];
        track->GetPxPyPz(p_vec);
        Int_t    charge       = track ->Charge();
        Double_t Track_phi    = track ->Phi();
	Double_t Track_theta  = track ->Theta();
	Double_t Track_eta    = track ->Eta();
	Double_t TPC_chi2     = track ->GetTPCchi2();
	Double_t TPC_signal   = track ->GetTPCsignal(); // dE/dx?
	Double_t TOF_signal   = track ->GetTOFsignal(); // time-of-flight?
        Double_t Track_length = track ->GetIntegratedLength();
	UShort_t N_TPC_cls    = track ->GetTPCNcls();

        Int_t pT_bin;
	for(Int_t i_pT = 0; i_pT < N_pT_bins; i_pT++)
	{
	    if(Track_pT >= pT_ranges[i_pT] && Track_pT < pT_ranges[i_pT+1])
	    {
                pT_bin = i_pT;
                break;
	    }
	}

	TLorentzVector TLV_p_vec;
	Double_t p_vec_energy = TMath::Sqrt(p_vec[0]*p_vec[0] + p_vec[1]*p_vec[1] + p_vec[2]*p_vec[2] + 0.938*0.938);
	TLV_p_vec.SetPxPyPzE(p_vec[0],p_vec[1],p_vec[2],p_vec_energy);
	//cout << "TLV_p_vec.P: " << TLV_p_vec.P() << ", P: " << Track_p << ", TLV_p_vec.Theta: " << TLV_p_vec.Theta() << ", Theta: " << Track_theta
	//<< ", TLV_p_vec.Phi: " << TLV_p_vec.Phi() << ", phi: " << Track_phi  << endl;


	ULong_t status = track->GetStatus();
	Int_t ITS_refit = 0;
        Int_t TPC_refit = 0;
        Int_t track_status = 0;
	if(((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit))
	{
	    ITS_refit = 1;
	    track_status |= 1 << 0; // setting bit 0 to 1
	}
	if(((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit))
	{
	    TPC_refit = 1;
	    track_status |= 1 << 1; // setting bit 1 to 1
	}


          // e = 0, muon = 1, pion = 2, kaon = 3, proton = 4
	Double_t Track_PID[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

        // nSigma TPC
	Track_PID[0] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
	Track_PID[1] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kMuon);
	Track_PID[2] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
	Track_PID[3] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
	Track_PID[4] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);

        // nSigma TOF, -999 in case there is no TOF hit
	Track_PID[5] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron);
	Track_PID[6] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kMuon);
	Track_PID[7] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
	Track_PID[8] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
	Track_PID[9] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);


	Float_t track_xy_impact,track_z_impact;
	track->GetImpactParameters(track_xy_impact,track_z_impact);
	Double_t track_total_impact = TMath::Sqrt(track_xy_impact*track_xy_impact + track_z_impact*track_z_impact);



	//-------------------
	Int_t N_ITS_cls = 0;
	for(Int_t i_ITS_layer = 0; i_ITS_layer < 6; ++i_ITS_layer)
	{
	    if(track ->HasPointOnITSLayer(i_ITS_layer))
	    {
		N_ITS_cls |= 1 << i_ITS_layer; // setting bit i_ITS_layer to 1
	    }
	}
	//-------------------

	TLorentzVector TL_vec;
	TL_vec.SetPtEtaPhiM(Track_pT,Track_eta,Track_phi,0.1349766);
       

        JetTrackParticle = JetTrackEvent->createParticle();
        JetTrackParticle ->set_dca_to_prim(((Double_t)charge)*track_total_impact);
        JetTrackParticle ->set_Particle_m2(0.0);
        JetTrackParticle ->set_Particle_nSigmaPi(Track_PID[2]);
        JetTrackParticle ->set_Particle_nSigmaK(Track_PID[3]);
        JetTrackParticle ->set_Particle_nSigmaP(Track_PID[4]);
        JetTrackParticle ->set_Particle_qp(0.0);
        JetTrackParticle ->set_Particle_hits_fit(0);
        JetTrackParticle ->set_TLV_Particle_prim(TL_vec);

	N_good_tracks++;

    } // End of TPC track loop
    //cout << "Tracks matched" << endl;


    Tree_AS_Event ->Fill();


    N_good_events++;
}



//________________________________________________________________________
void Ali_AS_analysis_TRD_digits::Terminate(Option_t *)
{
    cout << "In terminate" << endl;
}

