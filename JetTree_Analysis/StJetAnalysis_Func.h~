
#ifndef STJETANALYSIS_H
#define STJETANALYSIS_H




//----------------------------------------------------------------------------------------
bool sortFunc( const vector<Double_t>& p1,
              const vector<Double_t>& p2 )
{
    // for 2D vector sorting along column, third element
    return p1[2] < p2[2];
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Eff_track_rec_function(Double_t* x,Double_t* par)
{
    // Track reconstruction efficiency parametrization
    Double_t pt,y;
    Double_t A,B,C;A=par[0];B=par[1];C=par[2];
    pt=x[0];
    y=A*(exp(-pow(B/pt,C)));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// Momentum smearing
Int_t Apply_mom_smearing_and_efficiency(Int_t kflab_prim_glob, Float_t qp_In, Int_t kCentrality, Int_t i_Particle_use_In, Float_t m2_In,
                                        Float_t ePYTHIA_eff_factor_In, TLorentzVector TLV_Particle_In, PseudoJet& Fill_PseudoJet_smear_out,
                                        TF1* f_EfficiencyVsPt_In[][7][6])
{
    Double_t smear_factor = 0.01; // sigma_pT/pT = (1% for primaries, 2% for globals) * pT
    if(kflab_prim_glob == 1) smear_factor = 0.02; // global tracks are used
    Double_t pt_original = TLV_Particle_In.Pt();
    Double_t sigma_smear = smear_factor*pt_original*pt_original;
    Double_t pt_smear    = gRandom->Gaus(pt_original,sigma_smear);

    TLorentzVector TLV_Particle_smear = TLV_Particle_In;
    TLV_Particle_smear *= pt_smear/pt_original;

    PseudoJet Fill_PseudoJet_smear(TLV_Particle_smear.Px(),TLV_Particle_smear.Py(),TLV_Particle_smear.Pz(),TLV_Particle_smear.E());
    Fill_PseudoJet_smear.set_user_index(i_Particle_use_In);

    // New track efficiencies
    Int_t PID_eff = 2; // PID: p, anti-p, pi+, pi-, K+, K-
    if(qp_In > 0.0)
    {
        if(m2_In < 1.0 && m2_In > 0.8 ) PID_eff = 0; // p+
        if(m2_In < 0.3 && m2_In > 0.2 ) PID_eff = 4; // K+
        if(m2_In < 0.1 && m2_In > -0.1) PID_eff = 2; // pi+
    }
    if(qp_In < 0.0)
    {
        if(m2_In < 1.0 && m2_In > 0.8 ) PID_eff = 1; // p-
        if(m2_In < 0.3 && m2_In > 0.2 ) PID_eff = 5; // K-
        if(m2_In < 0.1 && m2_In > -0.1) PID_eff = 3; // pi-
    }


    // Centrality: 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80
    Int_t cent9_eff = 8;
    if(kCentrality == 0) // 0-10% -> randomly choose between efficiency for 0-5% and 5%-10%
    {
        cent9_eff = 0;
        if(ran_gen.Rndm() < 0.5) cent9_eff = 1;
    }
    if(kCentrality == 1) // 60-80% -> randomly choose between efficiency for 60%-70% and 70%-80%
    {
        cent9_eff = 7;
        if(ran_gen.Rndm() < 0.5) cent9_eff = 8;
    }


    Int_t RunId_In     = 0;
    Double_t epsilon = ePYTHIA_eff_factor_In*f_EfficiencyVsPt_In[cent9_eff][RunId_In][PID_eff]->Eval(TLV_Particle_smear.Pt());
    // cout << "pt = " << TLV_Particle_smear.Pt() << ", espilon = " << epsilon << ", ePYTHIA_eff_factor = " << ePYTHIA_eff_factor << endl;
    Double_t rnd = gRandom->Uniform(0,1);
    // cout << "pT = " << TLV_Particle_smear.Pt() << ", epsilon = " << epsilon << ", rnd = " << rnd << endl;

    if(rnd <= epsilon)
    {
        Fill_PseudoJet_smear_out = Fill_PseudoJet_smear;
        //cout << "Particle was accepted, efficiency: " << epsilon << ", pT: " << TLV_Particle_smear.Pt() << ", PID: " << PID_eff << endl;
        return 1;
    }
    else
    {
        //cout << "Particle was rejected due to track reconstruction efficiency" << endl;
        return 0;
    }
}
//----------------------------------------------------------------------------------------


/*
//----------------------------------------------------------------------------------------
Int_t Calc_Corr_EventPlane_Angle(StJetTrackEvent* JetTrackEvent_in, Double_t &EP_eta_pos, Double_t &EP_eta_neg,
                                Double_t &EP_Qx_eta_pos, Double_t &EP_Qy_eta_pos, Double_t &EP_Qx_eta_neg, Double_t &EP_Qy_eta_neg)
{
    // Calculates the Q-vector


    const Double_t nHitsFitA_EP_cut      = 14;
    const Double_t MomentumA_EP_low_cut  = 0.15;
    const Double_t MomentumA_EP_high_cut = 5.0;
    const Double_t dcaAB_EP_cut          = 2.0;
    const Double_t eta_EP_cut            = 1.0;
    const Double_t Event_radius_cut      = 2.0;
    const Double_t eta_gap               = 0.05; // 0.05

    // Event information
    Double_t  prim_vertex_x    = JetTrackEvent_in->getx();
    Double_t  prim_vertex_y    = JetTrackEvent_in->gety();
    Double_t  prim_vertex_z    = JetTrackEvent_in->getz();
    Int_t     RunId            = JetTrackEvent_in->getid();
    Double_t  refMult          = JetTrackEvent_in->getmult();
    Double_t  n_prim           = JetTrackEvent_in->getn_prim();
    Double_t  n_non_prim       = JetTrackEvent_in->getn_non_prim();
    Int_t     n_tofmatch_prim  = JetTrackEvent_in->getn_tof_prim();
    Double_t  ZDCx             = JetTrackEvent_in->getZDCx();
    Double_t  BBCx             = JetTrackEvent_in->getBBCx();
    Double_t  vzVPD            = JetTrackEvent_in->getvzVpd();
    Int_t     N_Particles      = JetTrackEvent_in->getNumParticle();
    Int_t     cent9            = JetTrackEvent_in->getcent9();

    Int_t z_bin = -1;
    if(prim_vertex_z >= -z_acceptance[gBeamTimeNum] && prim_vertex_z <= z_acceptance[gBeamTimeNum])
    {
        z_bin = (Int_t)((prim_vertex_z-start_z)/delta_z);

    }

    if(prim_vertex_x*prim_vertex_x + prim_vertex_y*prim_vertex_y > Event_radius_cut*Event_radius_cut) return 0;
    if(!(z_bin >= 0 && z_bin < n_z_vertex_bins)) return 0;

    EP_Qx_eta_pos        = 0.0;
    EP_Qy_eta_pos        = 0.0;
    Int_t Qtracks_used_eta_pos    = 0;
    EP_Qx_eta_neg        = 0.0;
    EP_Qy_eta_neg        = 0.0;
    Int_t Qtracks_used_eta_neg    = 0;

    for(Int_t i_Particle = 0; i_Particle < N_Particles; i_Particle++)
    {
        // Particle information
        StJetTrackParticle *JetTrackParticle_in = JetTrackEvent_in   ->getParticle(i_Particle);
        Double_t dca                     = JetTrackParticle_in->get_dca_to_prim();
        Double_t qp                      = JetTrackParticle_in->get_Particle_qp();
        Double_t nhitsfit                = JetTrackParticle_in->get_Particle_hits_fit();
        TLorentzVector TLV_Particle_glob = JetTrackParticle_in->get_TLV_Particle_glob();
        Double_t track_eta = TLV_Particle_glob.PseudoRapidity();
        Double_t track_phi = TLV_Particle_glob.Phi();
        Double_t track_pT  = TLV_Particle_glob.Pt();

        if(
           dca                < dcaAB_EP_cut
           && nhitsfit        > nHitsFitA_EP_cut
           && fabs(qp)        > MomentumA_EP_low_cut
           && fabs(qp)        < MomentumA_EP_high_cut
           && fabs(track_eta) < eta_EP_cut
           && fabs(track_eta) > eta_gap
          )
        {
            Double_t p_t_weight = 1.0;
            if(track_pT < 2.0)  p_t_weight = track_pT;
            if(track_pT >= 2.0) p_t_weight = 2.0;

            Double_t iQx = TMath::Cos(2.0*track_phi);
            Double_t iQy = TMath::Sin(2.0*track_phi);

            if(track_eta >= 0.0)
            {
                EP_Qx_eta_pos += p_t_weight*iQx;
                EP_Qy_eta_pos += p_t_weight*iQy;
                Qtracks_used_eta_pos++;
            }
            if(track_eta < 0.0)
            {
                EP_Qx_eta_neg += p_t_weight*iQx;
                EP_Qy_eta_neg += p_t_weight*iQy;
                Qtracks_used_eta_neg++;
            }

        }
    }


    Int_t file_bin_rc = (Int_t)(h_runId_index_rc_inv->GetBinContent(h_runId_index_rc_inv->FindBin(RunId)));

    // [Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full][z-bin][parameter]
    Double_t par0_eta_pos_Qx = h_rc_QxQy_etapm_z_vs_index_[0][z_bin][0]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[0][z_bin][0]->FindBin(file_bin_rc));
    Double_t par0_eta_pos_Qy = h_rc_QxQy_etapm_z_vs_index_[1][z_bin][0]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[1][z_bin][0]->FindBin(file_bin_rc));
    Double_t par0_eta_neg_Qx = h_rc_QxQy_etapm_z_vs_index_[2][z_bin][0]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[3][z_bin][0]->FindBin(file_bin_rc));
    Double_t par0_eta_neg_Qy = h_rc_QxQy_etapm_z_vs_index_[3][z_bin][0]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[4][z_bin][0]->FindBin(file_bin_rc));

    Double_t par1_eta_pos_Qx = h_rc_QxQy_etapm_z_vs_index_[0][z_bin][1]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[0][z_bin][1]->FindBin(file_bin_rc));
    Double_t par1_eta_pos_Qy = h_rc_QxQy_etapm_z_vs_index_[1][z_bin][1]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[1][z_bin][1]->FindBin(file_bin_rc));
    Double_t par1_eta_neg_Qx = h_rc_QxQy_etapm_z_vs_index_[2][z_bin][1]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[2][z_bin][1]->FindBin(file_bin_rc));
    Double_t par1_eta_neg_Qy = h_rc_QxQy_etapm_z_vs_index_[3][z_bin][1]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[3][z_bin][1]->FindBin(file_bin_rc));

    Double_t par2_eta_pos_Qx = h_rc_QxQy_etapm_z_vs_index_[0][z_bin][2]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[0][z_bin][2]->FindBin(file_bin_rc));
    Double_t par2_eta_pos_Qy = h_rc_QxQy_etapm_z_vs_index_[1][z_bin][2]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[1][z_bin][2]->FindBin(file_bin_rc));
    Double_t par2_eta_neg_Qx = h_rc_QxQy_etapm_z_vs_index_[2][z_bin][2]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[2][z_bin][2]->FindBin(file_bin_rc));
    Double_t par2_eta_neg_Qy = h_rc_QxQy_etapm_z_vs_index_[3][z_bin][2]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[3][z_bin][2]->FindBin(file_bin_rc));

    Double_t eta_pos_corr_Qx = par0_eta_pos_Qx + par1_eta_pos_Qx*refMult + par2_eta_pos_Qx*refMult*refMult;
    Double_t eta_pos_corr_Qy = par0_eta_pos_Qy + par1_eta_pos_Qy*refMult + par2_eta_pos_Qy*refMult*refMult;
    Double_t eta_neg_corr_Qx = par0_eta_neg_Qx + par1_eta_neg_Qx*refMult + par2_eta_neg_Qx*refMult*refMult;
    Double_t eta_neg_corr_Qy = par0_eta_neg_Qy + par1_eta_neg_Qy*refMult + par2_eta_neg_Qy*refMult*refMult;


    EP_Qy_eta_pos = (EP_Qy_eta_pos/Qtracks_used_eta_pos)-eta_pos_corr_Qy;
    EP_Qx_eta_pos = (EP_Qx_eta_pos/Qtracks_used_eta_pos)-eta_pos_corr_Qx;
    EP_Qy_eta_neg = (EP_Qy_eta_neg/Qtracks_used_eta_neg)-eta_neg_corr_Qy;
    EP_Qx_eta_neg = (EP_Qx_eta_neg/Qtracks_used_eta_neg)-eta_neg_corr_Qx;

    EP_eta_pos = TMath::ATan2(EP_Qy_eta_pos,EP_Qx_eta_pos);
    EP_eta_pos /= 2.0;
    EP_eta_neg = TMath::ATan2(EP_Qy_eta_neg,EP_Qx_eta_neg);
    EP_eta_neg /= 2.0;


    return 1;
}
//----------------------------------------------------------------------------------------
*/


//----------------------------------------------------------------------------------------
Double_t LevyFitFunc_pT(Double_t* x_val, Double_t* par)
{
    // One over pT term is removed -> original pT distribution
    // Fit function for d2N/(2pi dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt(pT*pT+m0*m0);
    y = pT*B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    return y;
}
//----------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
Double_t GaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0);
    return y;
}
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color

    if(x<0||y<0)
    {   // defaults
      x=gPad->GetLeftMargin()*1.15;
      y=(1-gPad->GetTopMargin())*1.04;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->Draw();
    return text;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x1_val);
    Zero_line -> SetX2(x2_val);
    Zero_line -> SetY1(y1_val);
    Zero_line -> SetY2(y2_val);
    Zero_line -> SetLineWidth(LineWidth);
    Zero_line -> SetLineStyle(LineStyle);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_STAR_3D()
{

    TGeoManager *geom = new TGeoManager("geom","My 3D Project");
    //------------------Creat materials------------------------------
    TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
    TGeoMaterial *Fe = new TGeoMaterial("Fe",55.84,26.7,7.87);
    Fe->SetTransparency(80); // higher value means more transparent, 100 is maximum


    TGeoMaterial *M_outer_tube = new TGeoMaterial("M_outer_tube",55.84,26.7,7.87);
    M_outer_tube->SetTransparency(93); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_IDS = new TGeoMaterial("M_IDS",55.84,26.7,7.87);
    M_IDS       ->SetTransparency(80); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_beampipe = new TGeoMaterial("M_beampipe",55.84,26.7,7.87);
    M_beampipe       ->SetTransparency(70); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_Pixel_support = new TGeoMaterial("M_Pixel_support",55.84,26.7,7.87);
    M_Pixel_support    ->SetTransparency(70); // higher value means more transparent, 100 is maximum


    //------------------Create media---------------------------------
    TGeoMedium *Air = new TGeoMedium("Air",0,vacuum);
    TGeoMedium *Iron = new TGeoMedium("Iron",1,Fe);
    TGeoMedium *Me_outer_tube = new TGeoMedium("Me_outer_tube",1,M_outer_tube);
    TGeoMedium *Me_IDS        = new TGeoMedium("Me_IDS",1,M_IDS);
    TGeoMedium *Me_beampipe   = new TGeoMedium("Me_beampipe",1,M_beampipe);
    TGeoMedium *Me_Pixel_support   = new TGeoMedium("Me_Pixel_support",1,M_Pixel_support);

    //------------------Create TOP volume----------------------------
    TGeoVolume *top = geom->MakeBox("top",Air,500,500,500);
    geom->SetTopVolume(top);
    geom->SetTopVisible(0);
    // If you want to see the boundary, please input the number, 1 instead of 0.
    // Like this, geom->SetTopVisible(1);


    TGeoVolume *inner_field_tube       = geom->MakeTube("inner_field_tube",Iron,49.5,50.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *outer_field_tube       = geom->MakeTube("outer_field_tube",Me_outer_tube,199.5,200.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_central_part       = geom->MakeTube("IDS_central_part",Me_IDS,42.8/2.0,43.0/2.0,56.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_side_parts         = geom->MakeTube("IDS_side_parts",Me_IDS,79.3/2.0,79.5/2.0,(222.7-64.0)/2.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_connection_parts_R = geom->MakeCone("IDS_connection_parts_R",Me_IDS,(64.0-56.0)/2.0,42.8/2.0,43.0/2.0,79.3/2.0,79.5/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max
    TGeoVolume *IDS_connection_parts_L = geom->MakeCone("IDS_connection_parts_L",Me_IDS,(64.0-56.0)/2.0,79.3/2.0,79.5/2.0,42.8/2.0,43.0/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max

    TGeoVolume *beampipe_central_part       = geom->MakeTube("beampipe_central_part",Me_beampipe,4.05/2.0,4.15/2.0,141.5);  // r_min, r_max, dz (half of total length)
    TGeoVolume *beampipe_side_parts         = geom->MakeTube("beampipe_side_parts",Me_beampipe,9.52/2.0,9.62/2.0,100.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *beampipe_connection_parts_R = geom->MakeCone("beampipe_connection_parts_R",Me_beampipe,(191.5-141.5)/2.0,4.05/2.0,4.15/2.0,9.52/2.0,9.62/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max
    TGeoVolume *beampipe_connection_parts_L = geom->MakeCone("beampipe_connection_parts_L",Me_beampipe,(191.5-141.5)/2.0,9.52/2.0,9.62/2.0,4.05/2.0,4.15/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max

    TGeoVolume *Pixel_support       = geom->MakeTube("Pixel_support",Me_Pixel_support,21.8/2.0,22.0/2.0,56.0);  // r_min, r_max, dz (half of total length)

    inner_field_tube       ->SetLineColor(4);
    outer_field_tube       ->SetLineColor(kRed-8);
    IDS_central_part       ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_side_parts         ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_connection_parts_R ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_connection_parts_L ->SetLineColor(2);  // Inner Detector Support (IDS)

    beampipe_central_part       ->SetLineColor(3);  // (beampipe)
    beampipe_side_parts         ->SetLineColor(3);  // (beampipe)
    beampipe_connection_parts_R ->SetLineColor(3);  // (beampipe)
    beampipe_connection_parts_L ->SetLineColor(3);  // (beampipe)

    Pixel_support ->SetLineColor(kYellow-3);  // (pixel support)

    top->AddNodeOverlap(inner_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(outer_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(IDS_central_part,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(IDS_side_parts,1,new TGeoTranslation(0,0,64.0+(222.7-64.0)/2.0));
    top->AddNodeOverlap(IDS_side_parts,1,new TGeoTranslation(0,0,-(64.0+(222.7-64.0)/2.0)));
    top->AddNodeOverlap(IDS_connection_parts_R,1,new TGeoTranslation(0,0,56.0+(64.0-56.0)/2.0));
    top->AddNodeOverlap(IDS_connection_parts_L,1,new TGeoTranslation(0,0,-(56.0+(64.0-56.0)/2.0)));

    top->AddNodeOverlap(beampipe_central_part,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(beampipe_side_parts,1,new TGeoTranslation(0,0,191.5+100.0));
    top->AddNodeOverlap(beampipe_side_parts,1,new TGeoTranslation(0,0,-(191.5+100.0)));
    top->AddNodeOverlap(beampipe_connection_parts_R,1,new TGeoTranslation(0,0,141.4+(191.5-141.5)/2.0));
    top->AddNodeOverlap(beampipe_connection_parts_L,1,new TGeoTranslation(0,0,-(141.4+(191.5-141.5)/2.0)));

    top->AddNodeOverlap(Pixel_support,1,new TGeoTranslation(0,0,0));
    /*
    top->DrawClone("ogl");



    const Int_t n_TPC_points = 50;
    TPolyLine3D   *TPC_endcaps[4];
    TPolyLine3D   *TPC_tube[4];
    TPolyLine3D   *TPC_tube_lines[n_TPC_points+1];

    Float_t radius_table[4] = {200,200,3.81,3.81};
    Float_t z_val_table[4]  = {200,-200,200,-200};

    Float_t radius_table_tube[4] = {50,50,50,50};
    Float_t z_val_table_tube[4]  = {200,-200,100,-100};

    for(Int_t r = 0; r < 4; r++)
    {
        TPC_endcaps[r] = new TPolyLine3D();
        Float_t radius   = radius_table[r];
        Float_t x_offset = 0.0;
        Float_t y_offset = 0.0;
        Float_t z_tpc_val   = z_val_table[r];
        for(Int_t t = 0; t < n_TPC_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            TPC_endcaps[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
        }
        TPC_endcaps[r]->SetLineStyle(0);
        TPC_endcaps[r]->SetLineColor(28); // 28
        TPC_endcaps[r]->SetLineWidth(2);
        TPC_endcaps[r]->DrawClone("ogl");
    }

    for(Int_t r = 0; r < 4; r++)
    {
        TPC_tube[r] = new TPolyLine3D();
        Float_t radius   = radius_table_tube[r];
        Float_t x_offset = 0.0;
        Float_t y_offset = 0.0;
        Float_t z_tpc_val   = z_val_table_tube[r];
        for(Int_t t = 0; t < n_TPC_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            TPC_tube[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
            if(r == 0 && (t%4 == 0))
            {
                TPC_tube_lines[t] = new TPolyLine3D();
                TPC_tube_lines[t]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
                TPC_tube_lines[t]->SetNextPoint(x_tpc_val,y_tpc_val,z_val_table_tube[r+1]);
                TPC_tube_lines[t]->SetLineStyle(0);
                TPC_tube_lines[t]->SetLineColor(28); // 28
                TPC_tube_lines[t]->SetLineWidth(1);
                //TPC_tube_lines[t]->DrawClone("ogl");
            }
        }
        TPC_tube[r]->SetLineStyle(0);
        TPC_tube[r]->SetLineColor(28); // 28
        TPC_tube[r]->SetLineWidth(2);
        TPC_tube[r]->DrawClone("ogl");
    }

    TPolyLine3D   *BeamLine;
    BeamLine       = new TPolyLine3D(2);
    BeamLine   ->SetPoint(0,0,0,-550);
    BeamLine   ->SetPoint(1,0,0,550);
    BeamLine   ->SetLineStyle(0);
    BeamLine   ->SetLineColor(4);
    BeamLine   ->SetLineWidth(2);
    BeamLine   ->DrawClone("ogl");
    */
}
//----------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------
Double_t efficiency11(Double_t pt, TF1* effLow, TF1* effHigh)
{
    // Efficiency for PYTHIA run11
    Double_t eff;
    if(pt<=1.2)eff = effLow->Eval(pt);
    else eff = effHigh->Eval(pt);
    //eff=eff+eff*0.05;
    return eff;
}
//----------------------------------------------------------------------------------------


#endif /* STJETANALYSIS_H */