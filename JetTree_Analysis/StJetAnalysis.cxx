#include "StJetAnalysis.h"

ClassImp(StJetAnalysis);



const char JETTRACK_EVENT_TREE[]   = "JetTrackEvent";
const char JETTRACK_EVENT_BRANCH[] = "Events";



static const Int_t N_max_jets = 120;
static Color_t jet_color[N_max_jets] = {2,kAzure-2,kGreen+1,kMagenta+2,kCyan+2,kOrange+7,kTeal-8,kYellow-2,kBlue-9,kPink+7,kCyan-3
,kRed+2,kViolet-2,kSpring-1,kYellow-6,kRed-2,kGreen-10,kBlue-3,kOrange-5,kCyan+4,kRed,kAzure-3,kGreen+2,kMagenta+3,kCyan+3,kOrange+8,kTeal-9,kYellow-3,
kBlue-8,kPink+8,kCyan-4,kRed+3,kViolet-3,kSpring-2,kYellow-7,kRed-3,kGreen-9,kBlue-4,kOrange-5,kCyan+5,
2,kAzure-2,kGreen+1,kMagenta+2,kCyan+2,kOrange+7,kTeal-8,kYellow-2,kBlue-9,kPink+7,kCyan-3
,kRed+2,kViolet-2,kSpring-1,kYellow-6,kRed-2,kGreen-10,kBlue-3,kOrange-5,kCyan+4,kRed,kAzure-3,kGreen+2,kMagenta+3,kCyan+3,kOrange+8,kTeal-9,kYellow-3,
kBlue-8,kPink+8,kCyan-4,kRed+3,kViolet-3,kSpring-2,kYellow-7,kRed-3,kGreen-9,kBlue-4,kOrange-5,kCyan+5,
2,kAzure-2,kGreen+1,kMagenta+2,kCyan+2,kOrange+7,kTeal-8,kYellow-2,kBlue-9,kPink+7,kCyan-3
,kRed+2,kViolet-2,kSpring-1,kYellow-6,kRed-2,kGreen-10,kBlue-3,kOrange-5,kCyan+4,kRed,kAzure-3,kGreen+2,kMagenta+3,kCyan+3,kOrange+8,kTeal-9,kYellow-3,
kBlue-8,kPink+8,kCyan-4,kRed+3,kViolet-3,kSpring-2,kYellow-7,kRed-3,kGreen-9,kBlue-4,kOrange-5,kCyan+5};



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



//------------------------------------------------------------------------------------------------------------------
vector<Double_t>  StJetAnalysis::Get_Helix_params_from_kine(TLorentzVector TLV_particle, TVector3 TV3_vertex, Double_t charge)
{
    // Returns helix parameters from kinematic input, Lorentz vector and vertex + charge

    Double_t x[3];
    Double_t p[3];
    x[0] = TV3_vertex.X();
    x[1] = TV3_vertex.Y();
    x[2] = TV3_vertex.Z();
    p[0] = TLV_particle.Px();
    p[1] = TLV_particle.Py();
    p[2] = TLV_particle.Pz();
    //calculation of Helixparameter taken from http://alidoc.cern.ch/AliRoot/v5-09-36/_ali_helix_8cxx_source.html
    //AliHelix::AliHelix(Double_t x[3], Double_t p[3], Double_t charge, Double_t conversion)
    vector<Double_t> fHelix;
    fHelix.resize(8);
    Double_t pt = TMath::Sqrt(p[0] * p[0] + p[1] * p[1]);
    //

    Double_t b_field	=	0.5;
    Double_t b_fak = b_field * 3. / 1000.;

    Double_t curvature = ((charge/pt) * b_fak);

    fHelix[4] = curvature; // C
    fHelix[3] = p[2] / pt; // tgl
    //
    Double_t xc, yc, rc;
    rc = 1 / fHelix[4];
    xc = x[0] - rc * p[1] / pt;
    yc = x[1] + rc * p[0] / pt;
    //
    fHelix[5] = x[0]; // x0
    fHelix[0] = x[1]; // y0
    fHelix[1] = x[2]; // z0
    //
    //fHelix[6] = xc;
    //fHelix[7] = yc;
    //fHelix[8] = TMath::Abs(rc);
    //
    fHelix[5] = xc;
    fHelix[0] = yc;

    fHelix[6] = pt;
    fHelix[7] = p[2];
    //
    if(TMath::Abs(p[1]) < TMath::Abs(p[0]))
    {
        fHelix[2] = TMath::ASin(p[1] / pt);
        //Helix[2]=asinf(p[1]/pt);
        if (charge * yc < charge * x[1])
            fHelix[2] = TMath::Pi() - fHelix[2];
    } else
    {
        fHelix[2] = TMath::ACos(p[0] / pt);
        //fHelix[2]=acosf(p[0]/pt);
        if (charge * xc > charge * x[0])
            fHelix[2] = -fHelix[2];
    }

    return fHelix;
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void  StJetAnalysis::Init3DGraphics()
{
    TEveManager::Create();

    TEveP_primary_vertex = new TEvePointSet();
    TEveLine_beam_axis = new TEveLine();
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,-650.0);
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,650.0);
    TEveLine_beam_axis ->SetName("beam axis");
    TEveLine_beam_axis ->SetLineStyle(1);
    TEveLine_beam_axis ->SetLineWidth(4);
    TEveLine_beam_axis ->SetMainAlpha(0.7);
    TEveLine_beam_axis ->SetMainColor(kBlue);

    TPL3D_helix = new TEveLine();
    gEve->AddElement(TEveLine_beam_axis);



    //--------------------------
    // Open histogram which defines good and bad chambers
    //TFile* file_TRD_QA = TFile::Open("./Data/chamber_QC.root");
    //h_good_bad_TRD_chambers = (TH1D*)file_TRD_QA ->Get("all_defects_hist");

    TH1I* h_good_bad_TRD_chambers = new TH1I("h_good_bad_TRD_chambers","h_good_bad_TRD_chambers",540,0,540);
    TFile* file_TRD_QA_flags = TFile::Open(GeomDir+"/chamber_QC_flags.root");
    vector<int> *t_flags;
    file_TRD_QA_flags ->GetObject("QC_flags", t_flags);

    // Its a 3 digit binary number. LSB is ADC good = 0 or bad = 1, next bit is anode HV good = 0, or bad = 1, and last bit is drift HV
    // so a 3 means that the ADC and the anode HV was bad, but the drift HV was okay

    // LSB = official QA, bit 1 = no fit, bit 2 = anode HV defect, bit 3 = drift HV defect, bit 4 = adc defect

    // number   adc defect   drift HV defect   anode HD defect    no fit   official QA
    //   0          0               0                0               0          0         --> all good
    //   1          0               0                0               0          1         --> official QA bad, rest good
    //  ...
    //   31         1               1                1               1          1         --> all bad

    Int_t i_chamber = 0;
    for(vector<int>::iterator it = t_flags->begin(); it != t_flags->end(); ++it)
    {
        //cout << "chamber: " << i_chamber << ", it: "  << *it << ", " << t_flags->at(i_chamber) << endl;
        h_good_bad_TRD_chambers ->SetBinContent(i_chamber+1,t_flags->at(i_chamber));
        i_chamber++;
    }
    //--------------------------

    //--------------------------
    // Load TRD geometry
    TFile* file_TRD_geom = TFile::Open(GeomDir+"TRD_Geom.root");
    vec_TH1D_TRD_geometry.resize(3); // x,y,z
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        vec_TH1D_TRD_geometry[i_xyz].resize(8); // 8 vertices
        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            HistName = "vec_TH1D_TRD_geometry_xyz_";
            HistName += i_xyz;
            HistName += "_V";
            HistName += i_vertex;
            vec_TH1D_TRD_geometry[i_xyz][i_vertex] = (TH1D*)file_TRD_geom->Get(HistName.Data());
        }

    }
#if defined(USEEVE)
    vec_eve_TRD_detector_box.resize(540);
#endif
    Int_t color_flag_QC[32];
    for(Int_t i_QC_flag = 0; i_QC_flag < 32; i_QC_flag++)
    {
        color_flag_QC[i_QC_flag] = kCyan;

        Int_t k_bit = 1; // fit
        Int_t bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // no fit
        {
            color_flag_QC[i_QC_flag] = kPink;
        }

        k_bit = 4; // ADC value
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // ADC low
        {
            color_flag_QC[i_QC_flag] = kMagenta;
        }

        k_bit = 2; // anode HV
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // anode HV low
        {
            color_flag_QC[i_QC_flag] = kYellow;
        }

        k_bit = 3; // drift HV bit
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // drift HV defect
        {
            color_flag_QC[i_QC_flag] = kOrange;
        }

        k_bit = 0; // official QA
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // official QA bad
        {
            color_flag_QC[i_QC_flag] = kRed;
        }
    }
    color_flag_QC[31] = kRed;
    //= {kCyan,kPink,kMagenta,kMagenta+2,kOrange,kOrange+2,kRed,kRed+2};

#if defined(USEEVE)
    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        vec_eve_TRD_detector_box[TRD_detector] = new TEveBox;

        HistName = "TRD_box_";
        HistName += TRD_detector;
        vec_eve_TRD_detector_box[TRD_detector] ->SetName(HistName.Data());
        Int_t flag_QC = h_good_bad_TRD_chambers ->GetBinContent(TRD_detector+1);
        if(!flag_QC) // chamber is OK flagged by QA
        {
            vec_eve_TRD_detector_box[TRD_detector]->SetMainColor(kCyan);
            vec_eve_TRD_detector_box[TRD_detector]->SetMainTransparency(95); // the higher the value the more transparent
        }
        else // bad chamber
        {
            vec_eve_TRD_detector_box[TRD_detector]->SetMainColor(color_flag_QC[flag_QC]);
            vec_eve_TRD_detector_box[TRD_detector]->SetMainTransparency(85); // the higher the value the more transparent
        }
        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            Double_t arr_pos_glb[3] = {vec_TH1D_TRD_geometry[0][i_vertex]->GetBinContent(TRD_detector),vec_TH1D_TRD_geometry[1][i_vertex]->GetBinContent(TRD_detector),vec_TH1D_TRD_geometry[2][i_vertex]->GetBinContent(TRD_detector)};
            vec_eve_TRD_detector_box[TRD_detector]->SetVertex(i_vertex,arr_pos_glb[0],arr_pos_glb[1],arr_pos_glb[2]);
        }

        //gEve->AddElement(vec_eve_TRD_detector_box[TRD_detector]);
    }

    gEve->Redraw3D(kTRUE);
#endif
    //--------------------------

}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
Long64_t StJetAnalysis::ReadData()
{
    printf("StJetAnalysis::ReadData() started \n");
    Int_t i_SE_ME = 0;


    //------------------------------------
    // Read input data
    TString SE_ME_List = InListDir+SEList;
    if (!SE_ME_List.IsNull())   // if input file is ok
    {
        cout << "Open file list " << SE_ME_List << endl;
        ifstream in(SE_ME_List);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_SE_ME[i_SE_ME]  = new TChain(JETTRACK_EVENT_TREE , JETTRACK_EVENT_TREE );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile = str;
                    addfile = pinputdir+addfile;
                    Long64_t file_entries;
                    input_SE_ME[i_SE_ME] ->AddFile(addfile.Data(),-1, JETTRACK_EVENT_TREE );
                    file_entries = input_SE_ME[i_SE_ME]->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }

            file_entries_SE_ME[0] = input_SE_ME[0]->GetEntries();
            cout << "Number of entries in chain: " << file_entries_SE_ME[0] << endl;

            JetTrackEvent = new StJetTrackEvent();
            input_SE_ME[i_SE_ME]  ->SetBranchAddress( JETTRACK_EVENT_BRANCH, &JetTrackEvent );
        }
    }
    else return 0;
    //------------------------------------

    printf("StJetAnalysis::ReadData() finished \n");
    return file_entries_SE_ME[0];
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void StJetAnalysis::InitJet()
{
    printf("StJetAnalysis::InitJet() started \n");
    Int_t i_SE_ME = 0;

    outputfile_name = eOutdir;
    outputfile_name += "F_jet_";
    outputfile_name += SEList;
    outputfile_name += ".root";
    outputfile_jet = new TFile(outputfile_name.Data(),"RECREATE");

    //------------------------------------
    h_jet_sub = new TH1D("h_jet_sub","h_jet_sub",1500,-50,100);
    //------------------------------------

    printf("StJetAnalysis::InitJet() finished \n");
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void StJetAnalysis::InitSplit()
{
    printf("StJetAnalysis::InitSplit() started \n");
    Int_t i_SE_ME = 0;
    r3.SetSeed(0);

    printf("N_z_vertex_bins: %d, N_mult_bins: %d, N_Psi_bins: %d \n",N_z_vertex_bins,N_mult_bins,N_Psi_bins);



    //------------------------------------
    // Define output trees
    JetTrackEvent_Fill.resize(N_z_vertex_bins);
    Tree_JetTrackEvent_Fill.resize(N_z_vertex_bins);
    vec_Outputfiles.resize(N_z_vertex_bins);
    for(Int_t i_z_vertex_bin = 0; i_z_vertex_bin < N_z_vertex_bins; i_z_vertex_bin++)
    {
        JetTrackEvent_Fill[i_z_vertex_bin].resize(N_mult_bins);
        Tree_JetTrackEvent_Fill[i_z_vertex_bin].resize(N_mult_bins);
        vec_Outputfiles[i_z_vertex_bin].resize(N_mult_bins);
        for(Int_t i_mult_bin = 0; i_mult_bin < N_mult_bins; i_mult_bin++)
        {
            JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin].resize(N_Psi_bins);
            Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin].resize(N_Psi_bins);
            vec_Outputfiles[i_z_vertex_bin][i_mult_bin].resize(N_Psi_bins);
            for(Int_t i_Psi_bin = 0; i_Psi_bin < N_Psi_bins; i_Psi_bin++)
            {

                //--------------
                // Define output files
                outputfile_name = eOutdir;
                outputfile_name += "F_z_";
                outputfile_name += i_z_vertex_bin;
                outputfile_name += "_mult_";
                outputfile_name += i_mult_bin;
                outputfile_name += "_Psi_";
                outputfile_name += i_Psi_bin;
                outputfile_name += SEList;
                outputfile_name += ".root";
                vec_Outputfiles[i_z_vertex_bin][i_mult_bin][i_Psi_bin] = new TFile(outputfile_name.Data(),"RECREATE");
                //--------------


                JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] = new StJetTrackEvent();
                Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] = NULL;

                //JetTrackEvent_ptr_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] = &JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin]; // old

                Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] = NULL;
                Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] = new TTree(JETTRACK_EVENT_TREE ,  JETTRACK_EVENT_TREE);
                Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin]->Branch(JETTRACK_EVENT_BRANCH , "StJetTrackEvent", JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin]);
                Long64_t maxtreesize = 2000000000;
                Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin]->SetMaxTreeSize(5*Long64_t(maxtreesize));
                Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin]->SetAutoSave( 10000000 );
                //Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin]->SetBasketSize("*",128000);
                Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin]->SetBasketSize("*",128000*10);
            }
        }
    }
    //------------------------------------



    printf("InitSplit finished \n");
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void StJetAnalysis::SplitTrees()
{
    printf("StJetAnalysis::SplitTrees started \n");

    Int_t i_SE_ME = 0;

    //--------------------
    for(Long64_t counter = 0; counter < file_entries_SE_ME[0]; counter++)
    {
        if (counter != 0  &&  counter % 100 == 0)
            cout << "." << flush;
        if (counter != 0  &&  counter % 1000 == 0)
        {
            if((file_entries_SE_ME[0]-0) > 0)
            {
                Double_t event_percent = 100.0*((Double_t)(counter-0))/((Double_t)(file_entries_SE_ME[0]-0));
                cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data " << flush;
            }
        }

        if (!input_SE_ME[i_SE_ME]->GetEntry( counter )) // take the event -> information is stored in event
            break;  // end of data chunk


        //--------------------
        // Event information
        Float_t prim_vertex_x   = JetTrackEvent->getx();
        Float_t prim_vertex_y   = JetTrackEvent->gety();
        Float_t prim_vertex_z   = JetTrackEvent->getz();
        Int_t   RunId           = JetTrackEvent->getid();
        Float_t refMult         = JetTrackEvent->getmult();
        Float_t n_prim          = JetTrackEvent->getn_prim();
        Int_t   SE_ME_flag      = JetTrackEvent->getSE_ME_flag();
        Float_t ZDCx            = JetTrackEvent->getZDCx();
        Float_t BBCx            = JetTrackEvent->getBBCx();
        Float_t vzVPD           = JetTrackEvent->getvzVpd();
        Int_t   N_Particles     = JetTrackEvent->getNumParticle();
        TVector2 QvecEtaPos     = JetTrackEvent->getQvecEtaPos();
        TVector2 QvecEtaNeg     = JetTrackEvent->getQvecEtaNeg();
        Int_t    cent9          = JetTrackEvent->getcent9();

        //printf("event: %lld, N_Particles: %d \n",counter,N_Particles);

        // Fill event information for jet in split categories
        Int_t i_z_vertex_bin = 0;
        Int_t i_mult_bin     = 0;
        Int_t i_Psi_bin      = 0;

        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->clearParticleList();
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setx(prim_vertex_x);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->sety(prim_vertex_y);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setz(prim_vertex_z);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setid(RunId);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setmult(refMult);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setn_prim(n_prim);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setSE_ME_flag(SE_ME_flag);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setZDCx(ZDCx);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setBBCx(BBCx);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setvzVpd(vzVPD);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setQvecEtaPos(QvecEtaPos);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setQvecEtaNeg(QvecEtaNeg);
        JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->setcent9(cent9);
        //--------------------



        //--------------------
        // Loop over track information
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


            JetTrackParticle_Fill = JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->createParticle();
            JetTrackParticle_Fill ->set_dca_to_prim(dca);
            JetTrackParticle_Fill ->set_Particle_m2(m2);
            JetTrackParticle_Fill ->set_Particle_nSigmaPi(nSPi);
            JetTrackParticle_Fill ->set_Particle_nSigmaK(nSK);
            JetTrackParticle_Fill ->set_Particle_nSigmaP(nSP);
            JetTrackParticle_Fill ->set_Particle_qp(qp);
            JetTrackParticle_Fill ->set_Particle_hits_fit(nhitsfit);
            JetTrackParticle_Fill ->set_TLV_Particle_prim(TLV_Particle_prim);
            JetTrackParticle_Fill ->set_TLV_Particle_glob(TLV_Particle_glob);
        }
        //--------------------

        Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin]  ->Fill();

    } // end of event loopo
    //--------------------


    //--------------------
    for(Int_t i_z_vertex_bin = 0; i_z_vertex_bin < N_z_vertex_bins; i_z_vertex_bin++)
    {
        for(Int_t i_mult_bin = 0; i_mult_bin < N_mult_bins; i_mult_bin++)
        {
            for(Int_t i_Psi_bin = 0; i_Psi_bin < N_Psi_bins; i_Psi_bin++)
            {
                vec_Outputfiles[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->cd();
                Tree_JetTrackEvent_Fill[i_z_vertex_bin][i_mult_bin][i_Psi_bin] ->Write("",TObject::kOverwrite);
            }
        }
    }
    //--------------------

}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void StJetAnalysis::MakeJets(Int_t graphics)
{
    //printf("StJetAnalysis::MakeJets() started \n");

    TVector3 TV3_beam(0.0,0.0,1.0);

    //--------------------
    vector<PseudoJet> jets_fiducial[2]; // [not smeared, smeared] only for PYTHIA
    Double_t jet_rho_array[2];

    // choose a jet definition
    JetDefinition jet_def(antikt_algorithm, eJet_R);

    // jet area definition
    Double_t ghost_maxrap = 1.0; // Fiducial cut for background estimation
    GhostedAreaSpec area_spec(ghost_maxrap);
    AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
    ClusterSequenceArea clust_seq_hard(vec_PJ_particles, jet_def, area_def);

    // run the clustering, extract the jets
    double ptmin = 0.2;
    vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));
    Selector Fiducial_cut_selector = SelectorAbsEtaMax(1.0 - eJet_R); // Fiducial cut for jets
    //Selector Fiducial_cut_selector = SelectorAbsEtaMax(1.0); // Fiducial cut for jets
    jets_fiducial[0] = Fiducial_cut_selector(jets_all);

    // print out some info
    //cout << "Clustered with " << jet_def.description() << endl;
    //--------------------



    //--------------------
    // background estimation
    //cout << "Define JetDefinition" << endl;
    JetDefinition jet_def_bkgd(kt_algorithm, eBkg_R); // <--
    //JetDefinition jet_def_bkgd(antikt_algorithm, jet_R); // test
    //cout << "Define AreaDefinition" << endl;
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
    //AreaDefinition area_def_bkgd(active_area,GhostedAreaSpec(ghost_maxrap,1,0.005));
    //cout << "Define selector" << endl;
    Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(eRem_n_hardest)); // <--
    //Selector selector = SelectorAbsEtaMax(1.0 - jet_R); // test

    //cout << "Define JetMedianBackgroundEstimator" << endl;
    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd); // <--
    //JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def, area_def); // test
    //cout << "Define Subtractor" << endl;

    Subtractor subtractor(&bkgd_estimator);
    //cout << "Define bkgd_estimator" << endl;
    bkgd_estimator.set_particles(vec_PJ_particles);

    //cout << "Calculate jet_rho and jet_sigma" << endl;
    //Double_t jet_rho   = bkgd_estimator.rho() - 0.212;   // TO BE CHANGED, equivallent to a 60 MeV shift for R = 0.3 (rho*pi*R*2 = 60 MeV)
    Double_t jet_rho   = bkgd_estimator.rho();
    jet_rho_array[0] = jet_rho;
    //cout << "jet_sigma" << endl;
    Double_t jet_sigma = bkgd_estimator.sigma();

    //printf("rho: %4.3f, sigma: %4.3f \n",jet_rho,jet_sigma);
    //--------------------



    //--------------------
    if(graphics)
    {
        vec_eve_Jets.clear();
        vec_eve_Jets.resize((Int_t)jets_fiducial[0].size());
        vec_eve_tracks_jets.clear();
    }

    Double_t Radius_jet_plot = 380.0;

    // Jet loop
    Int_t N_total_constituents = 0;
    if(jet_rho >= 0.0)
    {
        TLorentzVector TLV_jet;
        for(Int_t i_jet = 0; i_jet < (Int_t)jets_fiducial[0].size(); i_jet++)
        {
            Float_t jet_pt     = jets_fiducial[0][i_jet].perp();
            if(jet_pt <= 0.0) continue;
            Float_t jet_area   = jets_fiducial[0][i_jet].area();
            Float_t jet_pt_sub = jets_fiducial[0][i_jet].perp() - jet_rho*jet_area;
            Float_t jet_eta    = jets_fiducial[0][i_jet].eta();
            Float_t jet_phi    = jets_fiducial[0][i_jet].phi();

            TLV_jet.SetPtEtaPhiE(jet_pt,jet_eta,jet_phi,jet_pt);
            TVector3 TV3_jet(TLV_jet.Px(),TLV_jet.Py(),TLV_jet.Pz());
            TVector3 TV3_jet_dir(TLV_jet.Px(),TLV_jet.Py(),TLV_jet.Pz());
            TV3_jet_dir *= 1.0/TV3_jet_dir.Mag();
            Double_t perp_TV3_jet = TV3_jet.Perp();
            Double_t scaling_radius_plot = Radius_jet_plot/perp_TV3_jet;
            TV3_jet *= scaling_radius_plot;

            TVector3 TV3_jet_perpA = TV3_beam.Cross(TV3_jet_dir);
            TV3_jet_perpA *= 1.0/TV3_jet_perpA.Mag();
            TVector3 TV3_jet_perpB = TV3_jet_dir.Cross(TV3_jet_perpA);
            TV3_jet_perpB *= 1.0/TV3_jet_perpB.Mag();

            printf("i_jet: %d, jet_pt: %4.3f, jet_area: %4.3f, jet_eta: %4.3f, jet_phi: %4.3f \n",i_jet,jet_pt,jet_area,jet_eta,jet_phi);
            h_jet_sub ->Fill(jet_pt_sub);


            if(graphics)
            {
                vec_eve_Jets[i_jet] = new TEveBox;

                HistName = "Jet_";
                HistName += i_jet;
                vec_eve_Jets[i_jet] ->SetName(HistName.Data());
                vec_eve_Jets[i_jet] ->SetMainColor(kAzure+2);
                vec_eve_Jets[i_jet] ->SetMainTransparency(65); // the higher the value the more transparent

                TVector3 TV3_vertex[8];
                //Double_t cluster_energy_factor = TMath::Power(TLV_EMCal.E(),0.5)*10.0;
                Double_t cluster_energy_factor = jet_pt*3.0;
                Double_t cluster_width_factor = 6.0;
                TV3_vertex[0] = TV3_jet + cluster_width_factor*TV3_jet_perpA;
                TV3_vertex[1] = TV3_jet + cluster_width_factor*TV3_jet_perpA + cluster_width_factor*TV3_jet_perpB;
                TV3_vertex[2] = TV3_jet - cluster_width_factor*TV3_jet_perpA + cluster_width_factor*TV3_jet_perpB;
                TV3_vertex[3] = TV3_jet - cluster_width_factor*TV3_jet_perpA;
                TV3_vertex[4] = TV3_jet + cluster_width_factor*TV3_jet_perpA + TV3_jet_dir*cluster_energy_factor;
                TV3_vertex[5] = TV3_jet + cluster_width_factor*TV3_jet_perpA + cluster_width_factor*TV3_jet_perpB   + TV3_jet_dir*cluster_energy_factor;
                TV3_vertex[6] = TV3_jet - cluster_width_factor*TV3_jet_perpA + cluster_width_factor*TV3_jet_perpB   + TV3_jet_dir*cluster_energy_factor;
                TV3_vertex[7] = TV3_jet - cluster_width_factor*TV3_jet_perpA + TV3_jet_dir*cluster_energy_factor;

                for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
                {
                    Double_t arr_pos_glb[3] = {TV3_vertex[i_vertex][0],TV3_vertex[i_vertex][1],TV3_vertex[i_vertex][2]};
                    //printf("   i_vertex: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_vertex,TV3_vertex[i_vertex][0],TV3_vertex[i_vertex][1],TV3_vertex[i_vertex][2]);
                    vec_eve_Jets[i_jet]->SetVertex(i_vertex,arr_pos_glb[0],arr_pos_glb[1],arr_pos_glb[2]);
                }

                gEve->AddElement(vec_eve_Jets[i_jet]);

            }


            //--------------------
            // constituent track loop
            vector<PseudoJet> jet_constituents = jets_fiducial[0][i_jet].constituents();
            if(graphics)
            {
                N_total_constituents += (Int_t)jet_constituents.size();
                vec_eve_tracks_jets.resize(N_total_constituents);
            }

            TLorentzVector TLV_jet_const;
            for(Int_t i_constituent = 0; i_constituent < (Int_t)jet_constituents.size(); i_constituent++)
            {
                Float_t jet_const_pt  = jet_constituents[i_constituent].perp();
                if(jet_const_pt <= 0.0) continue;
                Float_t jet_const_phi = jet_constituents[i_constituent].phi();
                Float_t jet_const_eta = jet_constituents[i_constituent].eta();
                Int_t   user_index    = jet_constituents[i_constituent].user_index();

                //--------------------
                if(graphics)
                {
                    TLV_jet_const.SetPtEtaPhiE(jet_const_pt,jet_const_eta,jet_const_phi,jet_const_pt);
                    TVector3 TV3_jet_const(TLV_jet_const.Px(),TLV_jet_const.Py(),TLV_jet_const.Pz());
                    TVector3 TV3_jet_const_dir(TLV_jet_const.Px(),TLV_jet_const.Py(),TLV_jet_const.Pz());
                    TV3_jet_const_dir *= 1.0/TV3_jet_const_dir.Mag();
                    Double_t perp_TV3_jet_const = TV3_jet_const.Perp();
                    Double_t scaling_radius_plot_const = Radius_jet_plot/perp_TV3_jet_const;
                    TV3_jet_const *= scaling_radius_plot_const;

                    TVector3 TV3_jet_const_perpA = TV3_beam.Cross(TV3_jet_const_dir);
                    TV3_jet_const_perpA *= 1.0/TV3_jet_const_perpA.Mag();
                    TVector3 TV3_jet_const_perpB = TV3_jet_const_dir.Cross(TV3_jet_const_perpA);
                    TV3_jet_const_perpB *= 1.0/TV3_jet_const_perpB.Mag();

                    vec_eve_tracks_jets[i_constituent] = new TEveBox;

                    HistName = "Jet_const_";
                    HistName += i_constituent;


                    vec_eve_tracks_jets[i_constituent] ->SetName(HistName.Data());
                    if(i_jet < N_max_jets)
                    {
                        vec_eve_tracks_jets[i_constituent] ->SetMainColor(jet_color[i_jet]);
                    }
                    else vec_eve_tracks_jets[i_constituent] ->SetMainColor(kRed);
                    vec_eve_tracks_jets[i_constituent] ->SetMainTransparency(65); // the higher the value the more transparent

                    TVector3 TV3_vertex[8];
                    //Double_t cluster_energy_factor = TMath::Power(TLV_EMCal.E(),0.5)*10.0;
                    Double_t cluster_energy_factor = jet_const_pt*20.0;
                    Double_t cluster_width_factor = 5.0;
                    TV3_vertex[0] = TV3_jet_const + cluster_width_factor*TV3_jet_const_perpA;
                    TV3_vertex[1] = TV3_jet_const + cluster_width_factor*TV3_jet_const_perpA + cluster_width_factor*TV3_jet_const_perpB;
                    TV3_vertex[2] = TV3_jet_const - cluster_width_factor*TV3_jet_const_perpA + cluster_width_factor*TV3_jet_const_perpB;
                    TV3_vertex[3] = TV3_jet_const - cluster_width_factor*TV3_jet_const_perpA;
                    TV3_vertex[4] = TV3_jet_const + cluster_width_factor*TV3_jet_const_perpA + TV3_jet_const_dir*cluster_energy_factor;
                    TV3_vertex[5] = TV3_jet_const + cluster_width_factor*TV3_jet_const_perpA + cluster_width_factor*TV3_jet_const_perpB   + TV3_jet_const_dir*cluster_energy_factor;
                    TV3_vertex[6] = TV3_jet_const - cluster_width_factor*TV3_jet_const_perpA + cluster_width_factor*TV3_jet_const_perpB   + TV3_jet_const_dir*cluster_energy_factor;
                    TV3_vertex[7] = TV3_jet_const - cluster_width_factor*TV3_jet_const_perpA + TV3_jet_const_dir*cluster_energy_factor;

                    for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
                    {
                        Double_t arr_pos_glb[3] = {TV3_vertex[i_vertex][0],TV3_vertex[i_vertex][1],TV3_vertex[i_vertex][2]};
                        //printf("   i_vertex: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_vertex,TV3_vertex[i_vertex][0],TV3_vertex[i_vertex][1],TV3_vertex[i_vertex][2]);
                        vec_eve_tracks_jets[i_constituent]->SetVertex(i_vertex,arr_pos_glb[0],arr_pos_glb[1],arr_pos_glb[2]);
                    }

                    gEve->AddElement(vec_eve_tracks_jets[i_constituent]);
                }
                //--------------------
            }
            //--------------------


        } // end of jet loop
    }
    //--------------------

    if(graphics) gEve->Redraw3D(kTRUE);

}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
Int_t StJetAnalysis::LoopEvent(Int_t iEvent, Int_t graphics)
{
    Int_t i_SE_ME = 0;

    if (!input_SE_ME[i_SE_ME]->GetEntry( iEvent )) // take the event -> information is stored in event
        return 0;

    vec_PJ_particles.clear();

    //--------------------
    // Event information
    Float_t prim_vertex_x   = JetTrackEvent->getx();
    Float_t prim_vertex_y   = JetTrackEvent->gety();
    Float_t prim_vertex_z   = JetTrackEvent->getz();
    Int_t   RunId           = JetTrackEvent->getid();
    Float_t refMult         = JetTrackEvent->getmult();
    Float_t n_prim          = JetTrackEvent->getn_prim();
    Int_t   SE_ME_flag      = JetTrackEvent->getSE_ME_flag();
    Float_t ZDCx            = JetTrackEvent->getZDCx();
    Float_t BBCx            = JetTrackEvent->getBBCx();
    Float_t vzVPD           = JetTrackEvent->getvzVpd();
    Int_t   N_Particles     = JetTrackEvent->getNumParticle();
    UShort_t N_EMCal        = JetTrackEvent->getNumEMCal();
    TVector2 QvecEtaPos     = JetTrackEvent->getQvecEtaPos();
    TVector2 QvecEtaNeg     = JetTrackEvent->getQvecEtaNeg();
    Int_t    cent9          = JetTrackEvent->getcent9();

    TVector3 TV3_prim_vertex(prim_vertex_x,prim_vertex_y,prim_vertex_z);
    TVector3 TV3_beam(0.0,0.0,1.0);
    printf("N_Particles: %d, N_EMCal: %d \n",N_Particles,N_EMCal);
    //--------------------



    //--------------------
    // Loop over track information
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
        //TLorentzVector TLV_Particle_glob = JetTrackParticle->get_TLV_Particle_glob();

        PseudoJet Fill_PseudoJet(TLV_Particle_prim.Px(),TLV_Particle_prim.Py(),TLV_Particle_prim.Pz(),TLV_Particle_prim.E());
        vec_PJ_particles.push_back(Fill_PseudoJet);
    }
    //--------------------



    //--------------------
    // Loop over EMCal information
    Float_t pos_EMCal_clus[3] = {0.0};
    vec_eve_EMCal_cluster.clear();
    vec_eve_EMCal_cluster.resize(N_EMCal);
    for(Int_t i_EMCal = 0; i_EMCal < N_EMCal; i_EMCal++)
    {
        // Particle information
        JetEMCalParticle            = JetTrackEvent->getEMCal(i_EMCal);
        TLorentzVector TLV_EMCal    = JetEMCalParticle->get_TLV_EMCal();
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            pos_EMCal_clus[i_xyz] = JetEMCalParticle->get_cluster_pos(i_xyz);
        }

        TVector3 TV3_EMCal_cluster(pos_EMCal_clus[0],pos_EMCal_clus[1],pos_EMCal_clus[2]);
        TVector3 TV3_EMCal_cluster_dir(pos_EMCal_clus[0],pos_EMCal_clus[1],pos_EMCal_clus[2]);
        TV3_EMCal_cluster_dir *= 1.0/TV3_EMCal_cluster_dir.Mag();
        TVector3 TV3_EMCal_perpA = TV3_beam.Cross(TV3_EMCal_cluster_dir);
        TV3_EMCal_perpA *= 1.0/TV3_EMCal_perpA.Mag();
        TVector3 TV3_EMCal_perpB = TV3_EMCal_cluster_dir.Cross(TV3_EMCal_perpA);
        TV3_EMCal_perpB *= 1.0/TV3_EMCal_perpB.Mag();

        Float_t    Dx_to_track        = JetEMCalParticle->get_Dx_to_track();
        Float_t    Dz_to_track        = JetEMCalParticle->get_Dz_to_track();
        UShort_t   N_cells_in_cluster = JetEMCalParticle->get_N_cells_in_cluster();
        Bool_t     isExotic           = JetEMCalParticle->get_isExotic();
        Bool_t     isPHOS             = JetEMCalParticle->get_isPHOS();
        Bool_t     isEMCal            = JetEMCalParticle->get_isEMCal();

        //if(TLV_EMCal.E() > 5.0)
        {
            if(isEMCal) printf("EMCal, i_EMCal: %d, N_cells_in_cluster: %d, energy: %4.3f, iEvent: %d \n",i_EMCal,(Int_t)N_cells_in_cluster,TLV_EMCal.E(),iEvent);
        }
        //if(isPHOS)  printf("PHOS, i_EMCal: %d, N_cells_in_cluster: %d, energy: %4.3f, iEvent: %d \n",i_EMCal,(Int_t)N_cells_in_cluster,TLV_EMCal.E(),iEvent);


        if(graphics)
        {
            vec_eve_EMCal_cluster[i_EMCal] = new TEveBox;

            HistName = "EMCal_cluster_";
            HistName += i_EMCal;
            vec_eve_EMCal_cluster[i_EMCal] ->SetName(HistName.Data());
            if(isEMCal) vec_eve_EMCal_cluster[i_EMCal] ->SetMainColor(kRed);
            if(isPHOS)  vec_eve_EMCal_cluster[i_EMCal] ->SetMainColor(kGreen);
            vec_eve_EMCal_cluster[i_EMCal] ->SetMainTransparency(65); // the higher the value the more transparent

            TVector3 TV3_vertex[8];
            //Double_t cluster_energy_factor = TMath::Power(TLV_EMCal.E(),0.5)*10.0;
            Double_t cluster_energy_factor = TLV_EMCal.E()*10.0;
            if(isPHOS) cluster_energy_factor = TMath::Power(TLV_EMCal.E(),0.5)*1.0;
            Double_t cluster_width_factor = 5.0;
            TV3_vertex[0] = TV3_EMCal_cluster + cluster_width_factor*TV3_EMCal_perpA;
            TV3_vertex[1] = TV3_EMCal_cluster + cluster_width_factor*TV3_EMCal_perpA + cluster_width_factor*TV3_EMCal_perpB;
            TV3_vertex[2] = TV3_EMCal_cluster - cluster_width_factor*TV3_EMCal_perpA + cluster_width_factor*TV3_EMCal_perpB;
            TV3_vertex[3] = TV3_EMCal_cluster - cluster_width_factor*TV3_EMCal_perpA;
            TV3_vertex[4] = TV3_EMCal_cluster + cluster_width_factor*TV3_EMCal_perpA + TV3_EMCal_cluster_dir*cluster_energy_factor;
            TV3_vertex[5] = TV3_EMCal_cluster + cluster_width_factor*TV3_EMCal_perpA + cluster_width_factor*TV3_EMCal_perpB   + TV3_EMCal_cluster_dir*cluster_energy_factor;
            TV3_vertex[6] = TV3_EMCal_cluster - cluster_width_factor*TV3_EMCal_perpA + cluster_width_factor*TV3_EMCal_perpB   + TV3_EMCal_cluster_dir*cluster_energy_factor;
            TV3_vertex[7] = TV3_EMCal_cluster - cluster_width_factor*TV3_EMCal_perpA + TV3_EMCal_cluster_dir*cluster_energy_factor;

            for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
            {
                Double_t arr_pos_glb[3] = {TV3_vertex[i_vertex][0],TV3_vertex[i_vertex][1],TV3_vertex[i_vertex][2]};
                //printf("   i_vertex: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_vertex,TV3_vertex[i_vertex][0],TV3_vertex[i_vertex][1],TV3_vertex[i_vertex][2]);
                vec_eve_EMCal_cluster[i_EMCal]->SetVertex(i_vertex,arr_pos_glb[0],arr_pos_glb[1],arr_pos_glb[2]);
            }

            gEve->AddElement(vec_eve_EMCal_cluster[i_EMCal]);
        }
    }
    //--------------------



    //printf("N_Particles: %d \n",N_Particles);
    MakeJets(graphics);

    if(graphics) gEve->Redraw3D(kTRUE);

    return 1;
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void StJetAnalysis::WriteJet()
{
    outputfile_jet ->cd();
    h_jet_sub ->Write();
    outputfile_jet ->Close();
}
//------------------------------------------------------------------------------------------------------------------



