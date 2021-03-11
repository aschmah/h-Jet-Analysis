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
    Int_t   n_non_prim;
    Int_t   n_tof_prim;
    Int_t   SE_ME_flag;

    Float_t ZDCx;
    Float_t BBCx;
    Float_t vzVpd;

    TVector2 QvecEtaPos;
    TVector2 QvecEtaNeg;

    Int_t cent9;

    UShort_t      fNumParticle;

    TClonesArray* fParticle; //->

public:
    StJetTrackEvent() :
        x(-1),y(-1),z(-1),id(-1),mult(0),n_prim(0),n_non_prim(0),
        n_tof_prim(0),SE_ME_flag(-1),ZDCx(-1),BBCx(-1),vzVpd(-1),QvecEtaPos(),QvecEtaNeg(),cent9(0),fNumParticle(0)
    {
        fParticle      = new TClonesArray( "StJetTrackParticle", 10 );
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

        void       setn_non_prim(Int_t r)             { n_non_prim = r;                }
        Int_t      getn_non_prim() const              { return n_non_prim;             }

        void       setn_tof_prim(Int_t r)             { n_tof_prim = r;                }
        Int_t      getn_tof_prim() const              { return n_tof_prim;             }

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


        ClassDef(StJetTrackEvent,1)  // A simple event compiled of tracks
};
//------------------------------------------------------------------------------------



#endif // __STJETTRACKEVENT_H__
