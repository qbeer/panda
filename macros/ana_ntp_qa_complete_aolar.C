class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;
class RhoTuple;

void ana_ntp_qa_complete_aolar(int nevts = 5000)
{
  TString parAsciiFile = "all.par";
  TString prefix = "ppbar_psi2s2pi_jpsi2pi_2mu_5000";
  TString input = "psi2spipi_Jpsipipi_mumu.dec";
  TString output = "ana";
  TString friend1 = "pid";

  PndMasterRunAna *fRun = new PndMasterRunAna();
  fRun->SetInput(input);
  fRun->SetOutput(output);
  fRun->AddFriend(friend1);
  fRun->SetParamAsciiFile(parAsciiFile);
  fRun->Setup(prefix);

  gStyle->SetOptFit(1011);

  fRun->Init();

  TFile *out = TFile::Open(prefix + "_output_ana.root", "RECREATE");

  RhoTuple *npbarp = new RhoTuple("npbarp", "pbarp System");

  PndAnalysis *theAnalysis = new PndAnalysis();
  if (nevts == 0) nevts = theAnalysis->GetEntries();

  RhoCandList pbarp, muplus, muminus, piplus, piminus, jpsi, psi2s, jpsiFit, psi2sFit;

  double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass(); // Get nominal PDG mass of the J/psi
  RhoMassParticleSelector *jpsiMassSel = new RhoMassParticleSelector("jpsi", m0_jpsi, 1.0);

  TString pidSelection = "PidAlgoEmcBayes;PidAlgoDrc;PidAlgoDisc;PidAlgoStt;PidAlgoMdtHardCuts";

  // changed it accordingly to sim_complete beam momentum, inital lorentz vector of pbarpSystem
  // (p_x, p_y, p_z, E/c) style
  TLorentzVector ini1(0, 0, 9.231552, 9.279084);
  TLorentzVector ini2(0,0,0,0.938);
  TLorentzVector ini(ini1 + ini2);

  PndRhoTupleQA qa(theAnalysis, ini.P());

  int i = 0;
  while (theAnalysis->GetEvent() && i++ < nevts)
  {
    
    if(i%100 == 0) std::cout << "EVENT: #" << i << std::endl;

    theAnalysis->FillList(muplus, "MuonAllPlus", pidSelection);
    theAnalysis->FillList(muminus, "MuonAllMinus", pidSelection);
    theAnalysis->FillList(piplus, "PionAllPlus", pidSelection);
    theAnalysis->FillList(piminus, "PionAllMinus", pidSelection);

    jpsi.Combine(muplus, muminus);
    jpsi.SetType(443); // checked - j/psi

    for (int j = 0; j < jpsi.GetLength(); ++j)
    {
      
      RhoKinVtxFitter vtxfitter(jpsi[j]);
      vtxfitter.Fit();
      RhoCandidate* fitvtx_jpsi = jpsi[j]->GetFit();
      
      RhoKinFitter mfitter(fitvtx_jpsi);
			mfitter.AddMassConstraint(m0_jpsi);
			mfitter.Fit();
      fitvtx_jpsi = fitvtx_jpsi->GetFit();
      
      if (fitvtx_jpsi) jpsiFit.Add(fitvtx_jpsi);

    }

    psi2s.Combine(jpsiFit, piplus, piminus);
    psi2s.SetType(100443); // checked - psi(2s)

    for (int j = 0; j < psi2s.GetLength(); ++j)
    {

      RhoKinVtxFitter vtxfitter(psi2s[j]); 
      vtxfitter.Fit();
      RhoCandidate* fitted_psi2s = psi2s[j]->GetFit();
      
      if(fitted_psi2s) psi2sFit.Add(fitted_psi2s);

    }

    pbarp.Combine(psi2sFit, piplus, piminus);
    pbarp.SetType(88888); // checked - pbarpSystem

    for (int j = 0; j < pbarp.GetLength(); ++j){
      
      qa.qaComp("pbarp_", pbarp[j], npbarp);

      RhoDecayTreeFitter dtf(pbarp[j],ini);
      dtf.setMassConstraint(pbarp[j]->Daughter(0)->Daughter(0));
      
      int check = dtf.Fit();

      npbarp->Column("fitflag",(Float_t)check, -999.9f);
      npbarp->Column("fitchiq",(Float_t)dtf.chiSquare(), -999.9f);
      qa.qaComp("pbarpfit_", pbarp[j]->GetFit(), npbarp);
      
      npbarp->DumpData();
    }

    psi2sFit.Cleanup();
    jpsiFit.Cleanup();

  } // end while loop

  out->cd();

  npbarp->GetInternalTree()->Write();

  out->Save();

}
