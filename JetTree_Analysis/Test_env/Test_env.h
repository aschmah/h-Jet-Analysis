
using namespace std;

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "TString.h"

#include "TObject.h"

#include<TMath.h>

#include "Ali_TRD_Self_Event.h"
#include "Ali_TRD_Self_EventLinkDef.h"

#include "StJetTrackEvent.h"
#include "StJetTrackEventLinkDef.h"

ClassImp(Ali_TRD_Self_Event)
ClassImp(StJetTrackEvent)

//----------------------------------------------------------------------------------------
class Test_env
{
private:
    Ali_TRD_Self_Event*   TRD_Self_Event_out;

public:
    Test_env();
    //~Test_env();

    void Init();

    ClassDef(Test_env, 1)
};
//----------------------------------------------------------------------------------------

