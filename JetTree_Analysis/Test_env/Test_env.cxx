#include "Test_env.h"

Test_env::Test_env()
{

}

void Test_env::Init()
{
    cout << "Init started" << endl;
    Ali_TRD_Self_Event* event  = new Ali_TRD_Self_Event();
    cout << "Test A" << endl;
    StJetTrackEvent*    eventB = new StJetTrackEvent();
    cout << "Init finished" << endl;
}