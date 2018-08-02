int sim_complete(Int_t nEvents = 5000, TString  SimEngine ="TGeant3", Double_t BeamMomentum = 9.231552)
{
  TString parAsciiFile   = "all.par";
  
  TString prefix         = "ppbar_psi2s2pi_jpsi2pi_2mu_5000";     // prefix string for output files
  
  TString inputGenerator = "psi2spipi_Jpsipipi_mumu.dec";
 
  PndMasterRunSim *fRun = new PndMasterRunSim();
  
  fRun->SetInput(inputGenerator);
  fRun->SetName(SimEngine);
  fRun->SetParamAsciiFile(parAsciiFile);
  fRun->SetNumberOfEvents(nEvents);
  fRun->SetBeamMom(BeamMomentum);
  fRun->SetStoreTraj(kTRUE);
  fRun->Setup(prefix);
  fRun->CreateGeometry();
  fRun->SetGenerator();
  
  FairFilteredPrimaryGenerator *primGen = fRun->GetFilteredPrimaryGenerator();
  primGen->SetVerbose(0);
  
  fRun->AddSimTasks();
  fRun->Init();
  fRun->Run(nEvents); 
  fRun->Finish();
  
  return 0;
}

