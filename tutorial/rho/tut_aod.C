// Macro for running Panda digitization, reconstruction and pid tasks
// to run the macro:
// root  full_complete.C  or in root session root>.x  full_complete.C
void tut_aod(Int_t nEvents = 0, TString prefix = "signal")
{
  //-----User Settings:------------------------------------------------------
  TString  parAsciiFile   = "all.par";
  TString  output         = "pid";
  
  // -----   Initial Settings   --------------------------------------------
  PndMasterRunAna *fRun= new PndMasterRunAna();
  fRun->SetInput("dummy");
  fRun->SetOutput(output);
  fRun->SetParamAsciiFile(parAsciiFile);
  fRun->Setup(prefix);

  // -----   Add tasks   ----------------------------------------------------
  fRun->AddDigiTasks(kFALSE);
  fRun->AddRecoTasks(kFALSE);
  fRun->AddPidTasks();

  // -----   Intialise and run   --------------------------------------------
  fRun->Init();
  fRun->Run(0, nEvents);
  fRun->Finish();
}
