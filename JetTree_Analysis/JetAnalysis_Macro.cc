

R__LOAD_LIBRARY(StJetAnalysis_cxx.so);

void JetAnalysis_Macro(Int_t event_plot, TString suffix, Int_t BeamTime_num, Int_t Random_par,
                       Int_t Mode, Int_t In_Mode, Int_t z_bin, Int_t mult_bin, Int_t Psi_bin,
                       TString infilelist, Double_t PYTHIA_eff_factor, Int_t Centrality,
                       Int_t flab_prim_glob, Double_t Jet_R, Double_t Bkg_R, Int_t Rem_n_hardest, Double_t max_pt_down_scale
                      )
{

    Int_t iSystem = 0; // 0 = local, 1 = HD

    // First type:
    // gSystem ->Load("libfastjet.so");
    // gSystem ->Load("libfastjettools.0.so");
    // .L StJetAnalysis.cxx++

    // New tree structure, new parameters
    // root4star -b -q JetAnalysis_Macro.cc\(0,20000,\"_test\",6,-1,31,1,5,3,0,\"Jet_tree_V2_split_list_0_10_1651_1700\",1.0,0,0,0.3,0.3,3,3.0\)    // 0%-10%
    // root4star -b -q JetAnalysis_Macro.cc\(0,20000,\"_test\",6,-1,31,1,5,3,0,\"Jet_tree_V2_split_list_60_80_1751_1800\",1.0,1,0,0.3,0.3,3,3.0\)   // 60%-80%
    // root4star -b -q JetAnalysis_Macro.cc\(0,20000,\"_test\",6,-1,42,1,5,3,0,\"Jet_tree_V2_split_list_0_10_1651_1700\",1.0,0,0,0.3,0.3,3,3.0\)    // 0%-10%

    // root4star -b -q JetAnalysis_Macro.cc\(0,20000,\"_test\",6,-1,311,1,5,3,0,\"PYTHIA_list_hard_bin_55PtHard65\",1.0,0,0,0.3,0.3,0,3.0\)    // PYTHIA hard bin

    // root4star -b -q JetAnalysis_Macro.cc\(0,20000,\"_test\",6,-1,312,2,5,3,0,\"PYTHIA_list_mode312\",1.0,0,0,0.3,0.3,3,3.0\)    // PYTHIA hard bin
    // root4star -b -q JetAnalysis_Macro.cc\(0,20000,\"_test\",6,-1,312,1,5,3,0,\"PYTHIA_list_mode312\",1.0,0,0,0.3,0.3,3,3.0\)    // PYTHIA hard bin
    // root4star -b -q JetAnalysis_Macro.cc\(0,20000,\"_test\",6,-1,312,2,5,3,0,\"PYTHIA_list_mode312\",1.0,0,0,0.5,0.3,3,3.0\)    // PYTHIA hard bin

    // root -b -q JetAnalysis_Macro.cc\(0,20000,\"_test\",6,-1,1,2,5,3,0,\"Test_list.txt\",1.0,0,0,0.5,0.3,3,3.0\)
    // .x JetAnalysis_Macro.cc(0,20000,"_test",6,-1,1,2,5,3,0,"Test_list.txt",1.0,0,0,0.5,0.3,3,3.0)

    // New
    // .x JetAnalysis_Macro.cc(0,"_test",6,-1,1,2,5,3,0,"PbPb_run0_test.txt",1.0,0,0,0.5,0.3,3,3.0)

    // Mode
    // 0   = not important
    // 1   = create new trees subdivided into the different subcategories (EP, z-vertex,...) - First step
    // 4   = create new same event files with split tracks
    // 11  = create sample histos (do before mode 2) - Second step
    // 3   = jet finder (use Random_par = -1 for the first loop (SE and ME) to create the Et histograms)
    //       use In_Mode = 24 to load mixed event files with split tracks
    // 2   = create mixed event files (In_Mode 1 or 4)
    // 31  = recoil jet analysis
    // 32  = recoil jet analysis with high pT track suppression
    // 311 = recoil jet analysis - PYTHIA + calculate matching efficiency matrix
    // 312 = recoil jet analysis - embed PYTHIA jets into ME -> closure test
    // 42  = Delta pt calculation
    // 310 = m2 spectra from recoil jets

    // BeamTime_num:
    // 0 = 7.7
    // 1 = 11.5
    // 2 = 39
    // 3 = 62.4
    // 4 = 19.6
    // 5 = 27
    // 6 = 200
    // 7 = 14.5

    // Random_par
    // -1 = do not use weighting factors
    // 1 = use weighting factors

    // Centrality
    // 0 = 0-10%
    // 1 = 60-80%

    // flab_prim_glob
    // 0 = use primary tracks
    // 1 = use global tracks

    // Jet_R
    // Jet radius for anti-kt, usually 0.2, 0.3, 0.4, 0.5

    // Bkg_R
    // Jet radius for kt -> background, usually similar to Jet_R

    // Rem_n_hardest
    // how many of the hardest jet candidates should be removed from background rho calculation, usually 2-3 for SE and 0 for ME

    // max_pt_down_scale
    // at which track pT should the downscaling start, usually 3-4 GeV/c

    cout << "StJetAnalysis_Macro started" << endl;

    gSystem ->Load("libfastjet.so");
    gSystem ->Load("libfastjettools.0.so");
    gSystem ->Load("libTree");
    if(iSystem == 0)
    {
        gSystem ->Load("StJetAnalysis_cxx.so");
    }
    else
    {
        gSystem ->Load("/home/ceres/schmah/ALICE/Jet_Analysis/JetTree_Analysis/StJetAnalysis_cxx.so");
    }

    /*
    //gSystem ->Load("/global/homes/a/aschmah/local/lib/libfastjet.so");
    //gSystem ->Load("/global/homes/a/aschmah/local/lib/libfastjettools.so");
    gSystem ->Load("/misc/alisoft/alibuild/sw/ubuntu1604_x86-64/fastjet/v3.2.1_1.024-alice3-7/lib/libfastjet.so");
    gSystem ->Load("/misc/alisoft/alibuild/sw/ubuntu1604_x86-64/fastjet/v3.2.1_1.024-alice3-7/lib/libfastjettools.so");
    //gSystem ->Load("/u/aschmah/local/lib/libfastjet.so");
    //gSystem ->Load("/u/aschmah/local/lib/libfastjettools.so");
    gSystem ->Load("../../Utils/.sl64_gcc482/lib/StJetTrackEvent.so");
    //gSystem ->Load("StRefMultCorr"); // ALEX
    gSystem ->Load("StJetAnalysis");
    */

    //**************************** Set graphic style ***************************************
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,
                                                 greens, blues, NCont);
    gStyle->SetNumberContours(NCont);
    //**************************************************************************************

    cout << "Start of analysis" << endl;

    StJetAnalysis *StJetAnalysis_Ana = new StJetAnalysis();

    //---------------------------------------------------------------
    Int_t beamtime   = 1;
    Int_t graphics   = 1;
    //Int_t event_plot = 0; // -1 = all events
    StJetAnalysis_Ana->set_N_vertex_mult_Psi_bins(2,2,2);
    StJetAnalysis_Ana->setInListDir("File_lists/");
    StJetAnalysis_Ana->setSEList(infilelist.Data());
    if(iSystem == 0) // local
    {
        if(beamtime == 0)
        {
            StJetAnalysis_Ana->setInputDir("./Data/");
            StJetAnalysis_Ana->setOutdir("./Data/Jet_histos/");
            StJetAnalysis_Ana->setGeomDir("./Data/");
        }
        if(beamtime == 1)
        {
            StJetAnalysis_Ana->setInputDir("./Data/");
            StJetAnalysis_Ana->setOutdir("./Data/Jet_histos/");
            StJetAnalysis_Ana->setGeomDir("./Data/");
        }
    }
    else // HD
    {
        StJetAnalysis_Ana->setGeomDir("/home/ceres/schmah/ALICE/TRD_self_tracking/Data/");
        if(beamtime == 0)
        {
            StJetAnalysis_Ana->setInputDir("/misc/alidata120/alice_u/schmah/Jet_Analysis/pPb_2016/Track_trees/jet_trees_V2/");
            StJetAnalysis_Ana->setOutdir("/misc/alidata120/alice_u/schmah/Jet_Analysis/pPb_2016/Jet_histos/");
        }
        if(beamtime == 1)
        {
            StJetAnalysis_Ana->setInputDir("/misc/alidata120/alice_u/schmah/Jet_Analysis/PbPb_2018/Track_trees/jet_trees_V2/");
            StJetAnalysis_Ana->setOutdir("/misc/alidata120/alice_u/schmah/Jet_Analysis/PbPb_2018/Jet_histos/");
        }
    }


    //StJetAnalysis_Ana->InitSplit();
    //StJetAnalysis_Ana->SplitTrees(); // split the input trees into z-vertex, multiplicity and EP angle bins

    if(graphics) StJetAnalysis_Ana->Init3DGraphics();
    StJetAnalysis_Ana->InitJet();
    Long64_t Number_of_events = StJetAnalysis_Ana->ReadData();

    StJetAnalysis_Ana->setJet_R(0.4);
    StJetAnalysis_Ana->setBkg_R(0.3);
    StJetAnalysis_Ana->setRem_n_hardest(Rem_n_hardest);
    //StJetAnalysis_Ana->LoopEvent(11);

    Int_t return_LoopEvent = 0;

    Long64_t start_event = 0;
    Long64_t stop_event  = (Long64_t)Number_of_events;
    //Int_t stop_event  = 1;
    if(event_plot != -1)
    {
        start_event = event_plot;
        stop_event  = (event_plot+1);
    }


    for(Int_t i_event = start_event; i_event < stop_event; i_event++)
    {
        if (i_event != 0  &&  i_event % 100 == 0)
            cout << "." << flush;
        if (i_event != 0  &&  i_event % 1000 == 0)
        {
            if((Number_of_events-0) > 0)
            {
                Double_t event_percent = 100.0*((Double_t)(i_event-0))/((Double_t)(Number_of_events-0));
                cout << " " << i_event << " (" << event_percent << "%) " << "\n" << "==> Processing data " << flush;
            }
        }

        return_LoopEvent = StJetAnalysis_Ana->LoopEvent(i_event,graphics);
        if(!return_LoopEvent) break;
    }
    StJetAnalysis_Ana->WriteJet();
    //---------------------------------------------------------------



#if 0

    if(BeamTime_num == 6)
    {
        //StJetAnalysis_Ana->setEff_track_rec_parameter_file("/project/projectdirs/star/pwg/starjetc/aschmah/Jet/AuAu200_run11/Eff_parameters/EffParameters2011_New.root");
        //StJetAnalysis_Ana->Read_ABC_params_eff_track_rec_function();

        StJetAnalysis_Ana->setSEList(InListName.Data());
        StJetAnalysis_Ana->setInputDir("/misc/alidata120/alice_u/schmah/Jet_Analysis/pPb_2016/Track_trees/jet_trees_V2/");
        StJetAnalysis_Ana->setSuffix(suffix);
        StJetAnalysis_Ana->setOutputfile(Outfile.Data());

        /*
        StJetAnalysis_Ana->setInputDirPYTHIA("/project/projectdirs/star/pwg/starjetc/aschmah/Jet/AuAu200_run11/Pythia/Trees/");
        if(Mode == 312) StJetAnalysis_Ana->setInputDirPYTHIA("/project/projectdirs/star/pwg/starjetc/aschmah/Jet/AuAu200_run11/");
        StJetAnalysis_Ana->setSampleHisto("/project/projectdirs/star/pwg/starjetc/aschmah/Jet/AuAu200_run11/histo_out_V2/sample_histo/Jet_sample_histos_mode_11_full_run11_w_eta_phi_0_10_V2.root");
        if(Centrality == 1) StJetAnalysis_Ana->setSampleHisto("/project/projectdirs/star/pwg/starjetc/aschmah/Jet/AuAu200_run11/histo_out_V2/sample_histo/Jet_sample_histos_mode_11_full_run11_w_eta_phi_60_80_V2.root");

        StJetAnalysis_Ana->setSE_Et_Histo("/project/projectdirs/star/pwg/starjetc/aschmah/Jet/AuAu200_run11/histo_out_V2/sample_histo/Merge_Jet_200GeV_mode_3_In_mode_1_V20_R_03_0_10_V2.root");
        StJetAnalysis_Ana->setME_Et_Histo("/project/projectdirs/star/pwg/starjetc/aschmah/Jet/AuAu200_run11/histo_out_V2/sample_histo/Merge_Jet_200GeV_mode_3_In_mode_2_V20_R_03_0_10_V2_wo_two_hardest.root");
        if(Centrality == 1)
        {
            StJetAnalysis_Ana->setSE_Et_Histo("/project/projectdirs/star/pwg/starjetc/aschmah/Jet/AuAu200_run11/histo_out_V2/sample_histo/Merge_Jet_200GeV_mode_3_In_mode_1_R_03_60_80_V2.root");
            StJetAnalysis_Ana->setME_Et_Histo("/project/projectdirs/star/pwg/starjetc/aschmah/Jet/AuAu200_run11/histo_out_V2/sample_histo/Merge_Jet_200GeV_mode_3_In_mode_2_R_03_60_80_V2.root");
        }
        StJetAnalysis_Ana->setPYTHIA_eff_factor(PYTHIA_eff_factor);
        StJetAnalysis_Ana->setCentrality(Centrality);
        StJetAnalysis_Ana->setReCenteringFile("/global/homes/a/aschmah/STAR/ReCentering/Data/Recentering_FitParams_AuAu200.root");
        */
    }

    StJetAnalysis_Ana->setEnergy(BeamTime_num);
    StJetAnalysis_Ana->setRandom(Random_par);
    StJetAnalysis_Ana->setStartEvent(start_event);
    StJetAnalysis_Ana->setStopEvent(stop_event);
    StJetAnalysis_Ana->setMode(Mode);
    StJetAnalysis_Ana->setIn_Mode(In_Mode);
    StJetAnalysis_Ana->setOutdir(Outdir);
    StJetAnalysis_Ana->setBins(z_bin,mult_bin,Psi_bin);
    StJetAnalysis_Ana->setPrim_Glob(flab_prim_glob);
    StJetAnalysis_Ana->setJet_R(Jet_R);
    StJetAnalysis_Ana->setBkg_R(Bkg_R);
    StJetAnalysis_Ana->setRem_n_hardest(Rem_n_hardest);
    StJetAnalysis_Ana->set_Max_pt_down_scale(max_pt_down_scale);
    StJetAnalysis_Ana->Init();

#endif

    /*
    StJetAnalysis_Ana->Make();
    StJetAnalysis_Ana->Finish();
    */
    cout << "End of analysis" << endl;
}