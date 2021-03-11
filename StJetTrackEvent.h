#ifndef __STJETTRACKEVENT_H__
#define __STJETTRACKEVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector2.h"
//#include "StarClassLibrary/StThreeVectorF.hh"

// A. Schmah 05/31/2013
// A. Schmah 11/11/2013
// A. Schmah 02.03/2021 added ALICE EMCal container



//------------------------------------------------------------------------------------
class StEMCal : public TObject
{
    // EMCal data container
private:
    // Track properties
    TLorentzVector TLV_EMCal; // Lorentz vector properties of this EMCal cluster
    Float_t  cluster_pos[3];
    Float_t  Dx_to_track;
    Float_t  Dz_to_track;
    UShort_t N_cells_in_cluster;
    Bool_t   isExotic;

public:
    StEMCal() :
        TLV_EMCal(),cluster_pos(),Dx_to_track(0),Dz_to_track(0),N_cells_in_cluster(0),isExotic(0)
    {

    }
        ~StEMCal()
        {

        }

        // setters
        void set_TLV_EMCal(TLorentzVector tlv)                {TLV_EMCal  = tlv;    }
        void set_cluster_pos(Float_t x, Float_t y, Float_t z) {cluster_pos[0] = x; cluster_pos[1] = y; cluster_pos[2] = z;}
        void set_Dx_to_track(Float_t f)                       {Dx_to_track = f;}
        void set_Dz_to_track(Float_t f)                       {Dz_to_track = f;}
        void set_N_cells_in_cluster(UShort_t s)               {N_cells_in_cluster = s;}
        void set_isExotic(Bool_t b)                           {isExotic = b;}

        // getters
        TLorentzVector get_TLV_EMCal() const                  { return TLV_EMCal;   }
        Float_t        get_cluster_pos(Int_t index) const     { return cluster_pos[index];}
        Float_t        get_Dx_to_track() const                { return Dx_to_track;}
        Float_t        get_Dz_to_track() const                { return Dz_to_track;}
        UShort_t       get_N_cells_in_cluster() const         { return N_cells_in_cluster;}
        Bool_t         get_isExotic() const                   { return isExotic;}


        ClassDef(StEMCal,1)  //
};
//------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------
class StJetTrackParticle : public TObject
{
private:
    // Track properties
    TLorentzVector TLV_Particle_prim; // Lorentz vector properties of this particle (primary track)
    TLorentzVector TLV_Particle_glob; // Lorentz vector properties of this particle (global track)

    Float_t        dca_to_prim; // distance of closest approach to mother particle decay vertex (or primary vertex)
    Float_t        Particle_m2;
    Float_t        Particle_nSigmaPi;
    Float_t        Particle_nSigmaK;
    Float_t        Particle_nSigmaP;
    Float_t        Particle_qp;
    Float_t        Particle_hits_fit;

public:
    StJetTrackParticle() :
        TLV_Particle_prim(),TLV_Particle_glob(),dca_to_prim(-1),Particle_m2(-1),Particle_nSigmaPi(-1),
        Particle_nSigmaK(-1),Particle_nSigmaP(-1),Particle_qp(-1),Particle_hits_fit(-1)
    {

    }
        ~StJetTrackParticle()
        {

        }

        // setters
        void set_dca_to_prim(Float_t f)                     { dca_to_prim = f;            }
        void set_Particle_m2(Float_t f)                     { Particle_m2 = f;            }
        void set_Particle_nSigmaPi(Float_t f)               { Particle_nSigmaPi = f;      }
        void set_Particle_nSigmaK(Float_t f)                { Particle_nSigmaK = f;       }
        void set_Particle_nSigmaP(Float_t f)                { Particle_nSigmaP = f;       }
        void set_Particle_qp(Float_t f)                     { Particle_qp = f;            }
        void set_Particle_hits_fit(Float_t f)               { Particle_hits_fit = f;      }
        void set_TLV_Particle_prim(TLorentzVector tlv)      { TLV_Particle_prim = tlv;    }
        void set_TLV_Particle_glob(TLorentzVector tlv)      { TLV_Particle_glob = tlv;    }

        // getters
        Float_t get_dca_to_prim()              const        { return dca_to_prim;         }
        Float_t get_Particle_m2 ()             const        { return Particle_m2;         }
        Float_t get_Particle_nSigmaPi()        const        { return Particle_nSigmaPi;   }
        Float_t get_Particle_nSigmaK()         const        { return Particle_nSigmaK;    }
        Float_t get_Particle_nSigmaP()         const        { return Particle_nSigmaP;    }
        Float_t get_Particle_qp()              const        { return Particle_qp;         }
        Float_t get_Particle_hits_fit()        const        { return Particle_hits_fit;   }
        TLorentzVector get_TLV_Particle_prim() const        { return TLV_Particle_prim;   }
        TLorentzVector get_TLV_Particle_glob() const        { return TLV_Particle_glob;   }


        ClassDef(StJetTrackParticle,1)  // A simple track of a particle
};
//------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------
class StJetTrackEvent : public TObject
{
private:
    Float_t x;
    Float_t y;
    Float_t z;
    Int_t   id;
    Float_t mult;
    Int_t   n_prim;
    Int_t   n_TRD_tracklets;
    Int_t   n_tof_hits;
    Int_t   SE_ME_flag;

    Float_t ZDCx;
    Float_t BBCx;
    Float_t vzVpd;

    TVector2 QvecEtaPos;
    TVector2 QvecEtaNeg;

    Float_t   cent_class_ZNA; // ZDC neutral A
    Float_t   cent_class_ZNC; // ZDC neutral C
    Float_t   cent_class_V0A; // V0 A
    Float_t   cent_class_V0C; // V0 C
    Float_t   cent_class_V0M; // V0 average
    Float_t   cent_class_CL0; // clusters in layer 0
    Float_t   cent_class_CL1; // clusters in layer 1
    Float_t   cent_class_SPD; // SPD
    Float_t   cent_class_V0MEq; //
    Float_t   cent_class_V0AEq; //
    Float_t   cent_class_V0CEq; //

    Float_t BeamIntAA; // ZDC coincidence rate
    Float_t T0zVertex; // z-vertex position from VPD

    Int_t cent9;
    TString TriggerWord; // Trigger word

    UShort_t      fNumParticle;
    UShort_t      fNumEMCal;

    TClonesArray* fParticle; //->
    TClonesArray* fEMCal; //->

public:
    StJetTrackEvent() :
        x(-1),y(-1),z(-1),id(-1),mult(0),n_prim(0),n_TRD_tracklets(0),
        n_tof_hits(0),SE_ME_flag(-1),ZDCx(-1),BBCx(-1),vzVpd(-1),QvecEtaPos(),QvecEtaNeg(),
        cent_class_ZNA(0),cent_class_ZNC(0),cent_class_V0A(0),cent_class_V0C(0),cent_class_V0M(0),cent_class_CL0(0),cent_class_CL1(0),
        cent_class_SPD(0),cent_class_V0MEq(0),cent_class_V0AEq(0),cent_class_V0CEq(0),BeamIntAA(-1),T0zVertex(-1),cent9(0),
        TriggerWord(),fNumParticle(0),fNumEMCal(0)
    {
        fParticle      = new TClonesArray( "StJetTrackParticle", 10 );
        fEMCal         = new TClonesArray( "StEMCal", 10 );
    }
        ~StJetTrackEvent()
        {
            delete fParticle;
            fParticle = NULL;
        }

        void       setx(Float_t r)                    { x = r;                         }
        Float_t    getx() const                       { return x;                      }

        void       sety(Float_t r)                    { y = r;                         }
        Float_t    gety() const                       { return y;                      }

        void       setz(Float_t r)                    { z = r;                         }
        Float_t    getz() const                       { return z;                      }

        void       setid(Int_t  r)                    { id = r;                        }
        Int_t      getid() const                      { return id;                     }

        void       setmult(Float_t r)                 { mult = r;                      }
        Float_t    getmult() const                    { return mult;                   }

        void       setn_prim(Int_t r)                 { n_prim = r;                    }
        Int_t      getn_prim() const                  { return n_prim;                 }

        void       setn_TRD_tracklets(Int_t r)        { n_TRD_tracklets = r;                }
        Int_t      getn_TRD_tracklets() const         { return n_TRD_tracklets;             }

        void       setn_tof_hits(Int_t r)             { n_tof_hits = r;                }
        Int_t      getn_tof_hits() const              { return n_tof_hits;             }

        void       setSE_ME_flag(Int_t r)             { SE_ME_flag = r;                }
        Int_t      getSE_ME_flag() const              { return SE_ME_flag;             }

        void       setZDCx(Float_t r)                 { ZDCx = r;                      }
        Float_t    getZDCx() const                    { return ZDCx;                   }

        void       setBBCx(Float_t r)                 { BBCx = r;                      }
        Float_t    getBBCx() const                    { return BBCx;                   }

        void       setvzVpd(Float_t r)                { vzVpd = r;                     }
        Float_t    getvzVpd() const                   { return vzVpd;                  }

        void       setQvecEtaPos(TVector2 vec)        { QvecEtaPos = vec;              }
        TVector2   getQvecEtaPos() const              { return QvecEtaPos;             }

        void       setQvecEtaNeg(TVector2 vec)        { QvecEtaNeg = vec;              }
        TVector2   getQvecEtaNeg() const              { return QvecEtaNeg;             }

        void       setcent9(Int_t r)                  { cent9 = r;                     }
        Int_t      getcent9() const                   { return cent9;                  }

        void       setcent_class_ZNA(Float_t r)             { cent_class_ZNA = r;                }
	Float_t      getcent_class_ZNA() const              { return cent_class_ZNA;             }

	void       setcent_class_ZNC(Float_t r)             { cent_class_ZNC = r;                }
	Float_t      getcent_class_ZNC() const              { return cent_class_ZNC;             }

	void       setcent_class_V0A(Float_t r)             { cent_class_V0A = r;                }
	Float_t      getcent_class_V0A() const              { return cent_class_V0A;             }

	void       setcent_class_V0C(Float_t r)             { cent_class_V0C = r;                }
	Float_t      getcent_class_V0C() const              { return cent_class_V0C;             }

	void       setcent_class_V0M(Float_t r)             { cent_class_V0M = r;                }
	Float_t      getcent_class_V0M() const              { return cent_class_V0M;             }

	void       setcent_class_CL0(Float_t r)             { cent_class_CL0 = r;                }
	Float_t      getcent_class_CL0() const              { return cent_class_CL0;             }

	void       setcent_class_CL1(Float_t r)             { cent_class_CL1 = r;                }
	Float_t      getcent_class_CL1() const              { return cent_class_CL1;             }

	void       setcent_class_SPD(Float_t r)             { cent_class_SPD = r;                }
	Float_t      getcent_class_SPD() const              { return cent_class_SPD;             }

	void       setcent_class_V0MEq(Float_t r)             { cent_class_V0MEq = r;                }
	Float_t      getcent_class_V0MEq() const              { return cent_class_V0MEq;             }

	void       setcent_class_V0AEq(Float_t r)             { cent_class_V0AEq = r;                }
	Float_t      getcent_class_V0AEq() const              { return cent_class_V0AEq;             }

	void       setcent_class_V0CEq(Float_t r)             { cent_class_V0CEq = r;                }
	Float_t      getcent_class_V0CEq() const              { return cent_class_V0CEq;             }

	void       setBeamIntAA(Float_t r)                 { BeamIntAA = r;                      }
	Float_t    getBeamIntAA() const                    { return BeamIntAA;                   }

	void       setT0zVertex(Float_t r)            { T0zVertex = r;                     }
	Float_t    getT0zVertex() const               { return T0zVertex;                  }


        void       setTriggerWord(TString s)          { TriggerWord = s;}
        TString    getTriggerWord() const             { return TriggerWord; }

        //--------------------------------------
        // Particle
        StJetTrackParticle* createParticle()
        {
            if (fNumParticle == fParticle->GetSize())
                fParticle->Expand( fNumParticle + 5 );
            if (fNumParticle >= 50000)
            {
                Fatal( "StJetTrackParticle::createParticle()", "ERROR: Too many Particles (>5000)!" );
                exit( 2 );
            }

            new((*fParticle)[fNumParticle++]) StJetTrackParticle;
            return (StJetTrackParticle*)((*fParticle)[fNumParticle - 1]);
        }
        void clearParticleList()
        {
            fNumParticle   = 0;
            fParticle      ->Clear();
        }
        UShort_t getNumParticle() const
        {
            return fNumParticle;
        }
        StJetTrackParticle* getParticle(UShort_t i) const
        {
            return i < fNumParticle ? (StJetTrackParticle*)((*fParticle)[i]) : NULL;
        }
        //--------------------------------------



        //--------------------------------------
        // EMCal
        StEMCal* createEMCal()
        {
            if (fNumEMCal == fEMCal->GetSize())
                fEMCal->Expand( fNumEMCal + 5 );
            if (fNumEMCal >= 50000)
            {
                Fatal( "StEMCal::createEMCal()", "ERROR: Too many EMCals (>5000)!" );
                exit( 2 );
            }

            new((*fEMCal)[fNumEMCal++]) StEMCal;
            return (StEMCal*)((*fEMCal)[fNumEMCal - 1]);
        }
        void clearEMCalList()
        {
            fNumEMCal   = 0;
            fEMCal      ->Clear();
        }
        UShort_t getNumEMCal() const
        {
            return fNumEMCal;
        }
        StEMCal* getEMCal(UShort_t i) const
        {
            return i < fNumEMCal ? (StEMCal*)((*fEMCal)[i]) : NULL;
        }
        //--------------------------------------


        ClassDef(StJetTrackEvent,1)  // A simple event compiled of tracks
};
//------------------------------------------------------------------------------------



#endif // __STJETTRACKEVENT_H__
