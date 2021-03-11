#ifndef ALI_AS_ANALYSIS_TRD_DIGITS_H
#define ALI_AS_ANALYSIS_TRD_DIGITS_H

class AliTRDdigitsManager;

#include "AliAnalysisTaskSE.h"

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "StJetTrackEvent.h"
#include "StJetTrackEventLinkDef.h"
#include "TVector2.h"

ClassImp(StEMCal)
ClassImp(StJetTrackParticle)
ClassImp(StJetTrackEvent)



class Ali_AS_analysis_TRD_digits : public AliAnalysisTaskSE
{
public:
    Ali_AS_analysis_TRD_digits()
	: AliAnalysisTaskSE(),
	JetTrackEvent(0),JetTrackParticle(0),JetEMCal(0),Tree_AS_Event(0), fEventNoInFile(-2), N_good_events(0), fDigitsLoadedFlag(kFALSE),
        fCorrTaskSetting(""),fInputEvent(NULL),h_dca(0x0),h_dca_xyz(0x0),h2D_TPC_dEdx_vs_momentum(0x0),h_ADC_tracklet(0x0),h_ADC_vs_time(0x0)
    {
	cout << "" << endl;
	cout << "***************************************************************************************" << endl;
	cout << "In Ali_AS_analysis_TRD_digits.h constructor" << endl;
	cout << "***************************************************************************************" << endl;
	cout << "" << endl;
    }
	Ali_AS_analysis_TRD_digits(const char *name);
	//virtual ~Ali_AS_analysis_TRD_digits() {}

	virtual void   UserCreateOutputObjects();
	virtual Bool_t UserNotify();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);

	AliHelix aliHelix;

    protected:

	Bool_t NextEvent(Bool_t preload=kFALSE);
	TList           *fListOfHistos;       //! list of output histograms
	TTree           *fTree;               //! output tree
	AliPIDResponse  *fPIDResponse;        //! PID handling
	AliESDtrackCuts *EsdTrackCuts;        //!

    private:

        StJetTrackEvent*    JetTrackEvent;
        StJetTrackParticle* JetTrackParticle;
        StEMCal* JetEMCal;
	TTree       *Tree_AS_Event;

	Int_t fEventNoInFile;
	Int_t N_good_events;
        Int_t fDigitsLoadedFlag;
        TString fCorrTaskSetting; // EMCal
        AliVEvent*            fInputEvent;                                          // current event

	std::vector<TH1D*> h_dca;
	std::vector< std::vector<TH1D*> > h_dca_xyz;
        TH2D* h2D_TPC_dEdx_vs_momentum;
        vector<TH1D*> h_ADC_tracklet;
        vector<TProfile*> h_ADC_vs_time;

	Ali_AS_analysis_TRD_digits(const Ali_AS_analysis_TRD_digits&); // not implemented
	Ali_AS_analysis_TRD_digits& operator=(const Ali_AS_analysis_TRD_digits&); // not implemented

	ClassDef(Ali_AS_analysis_TRD_digits, 1);
};



#endif
