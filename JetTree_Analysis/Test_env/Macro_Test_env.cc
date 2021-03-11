R__LOAD_LIBRARY(Test_env_cxx.so);

void Macro_Test_env()
{
    gSystem ->Load("Test_env_cxx.so");

    Test_env* blubb = new Test_env();
    blubb ->Init();
}