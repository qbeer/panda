class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;


void tut_ana_mcmatch(int nevts = 0, TString prefix = "signal")
{
 	// *** some variables
	int i=0,j=0, k=0, l=0;
	gStyle->SetOptFit(1011);
	
	// *** the output file for FairRunAna
	TString OutFile="out_dummy.root";  
					
	// *** the files coming from the simulation
	TString inPidFile  = prefix+"_pid.root";    // this file contains the PndPidCandidates and McTruth
	TString inParFile  = prefix+"_par.root";
	
	// *** PID table with selection thresholds; can be modified by the user
	TString pidParFile = TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all.par";	
	
	// *** initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	FairRunAna* fRun = new FairRunAna();
	FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
	fRun->SetSource(new FairFileSource(inPidFile));
	
	// *** setup parameter database 	
	FairParRootFileIo* parIO = new FairParRootFileIo();
	parIO->open(inParFile);
	FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
	parIOPid->open(pidParFile.Data(),"in");
	
	rtdb->setFirstInput(parIO);
	rtdb->setSecondInput(parIOPid);
	rtdb->setOutput(parIO);  
	
	fRun->SetOutputFile(OutFile);
	fRun->Init(); 
	
	// *** create an output file for all histograms
	TFile *out = TFile::Open(prefix+"_ana_mcmatch.root","RECREATE");
	
	// *** create some histograms
	TH1F *hjpsim_all = new TH1F("hjpsim_all","J/#psi mass",200,0,4.5);
	TH1F *hpsim_all  = new TH1F("hpsim_all","#psi(2S) mass",200,0,5);
	
	TH1F *hjpsim_nm = new TH1F("hjpsim_nm","J/#psi mass (no truth match)",200,0,4.5);
	TH1F *hpsim_nm  = new TH1F("hpsim_nm","#psi(2S) mass (no truth match)",200,0,5);
	
	TH1F *hjpsim_ftm = new TH1F("hjpsim_ftm","J/#psi mass (full truth match)",200,0,4.5);
	TH1F *hpsim_ftm  = new TH1F("hpsim_ftm","#psi(2S) mass (full truth match)",200,0,5);
	
	TH1F *hjpsim_diff = new TH1F("hjpsim_diff","J/#psi mass diff to truth",100,-2,2);
	TH1F *hpsim_diff  = new TH1F("hpsim_diff","#psi(2S) mass diff to truth",100,-2,2);
	
	//
	// Now the analysis stuff comes...
	//
	
	
	// *** the data reader object
	PndAnalysis* theAnalysis = new PndAnalysis();
	if (nevts==0) nevts= theAnalysis->GetEntries();
	
	// *** RhoCandLists for the analysis
	RhoCandList muplus, muminus, piplus, piminus, jpsi, psi2s;
	
	// *** Mass selector for the jpsi cands
	RhoMassParticleSelector *jpsiMassSel=new RhoMassParticleSelector("jpsi",3.096,1.0);
	
	// ***
	// the event loop
	// ***
	while (theAnalysis->GetEvent() && i++<nevts)
	{
		if ((i%100)==0) cout<<"evt " << i << endl;
		
		// *** Select with no PID info ('All'); type and mass are set 		
		theAnalysis->FillList(muplus,  "MuonAllPlus");
		theAnalysis->FillList(muminus, "MuonAllMinus");
		theAnalysis->FillList(piplus,  "PionAllPlus");
		theAnalysis->FillList(piminus, "PionAllMinus");
		
		// *** combinatorics for J/psi -> mu+ mu-
		jpsi.Combine(muplus, muminus);
		
		// ***
		// *** do the TRUTH MATCH for jpsi
		// ***
		jpsi.SetType(443);
				
		for (j=0;j<jpsi.GetLength();++j) 
		{
			hjpsim_all->Fill( jpsi[j]->M() );
			
			if (theAnalysis->McTruthMatch(jpsi[j])) 
			{
				hjpsim_ftm->Fill( jpsi[j]->M() );
				hjpsim_diff->Fill( jpsi[j]->GetMcTruth()->M() - jpsi[j]->M() );
			}
			else 
				hjpsim_nm->Fill( jpsi[j]->M() );
		}
				
		// *** some rough mass selection
		jpsi.Select(jpsiMassSel);
		
		// *** combinatorics for psi(2S) -> J/psi pi+ pi-
		psi2s.Combine(jpsi, piplus, piminus);
		
		// ***
		// *** do the TRUTH MATCH for psi(2S)
		// ***
		psi2s.SetType(88888);

		for (j=0;j<psi2s.GetLength();++j) 
		{
			hpsim_all->Fill( psi2s[j]->M() );
			
			if (theAnalysis->McTruthMatch(psi2s[j]))
			{ 
			 	hpsim_ftm->Fill( psi2s[j]->M() );
				hpsim_diff->Fill( psi2s[j]->GetMcTruth()->M() - psi2s[j]->M() );
			}
			else 
				hpsim_nm->Fill( psi2s[j]->M() );
		}			
	}
	
	// *** write out all the histos
	out->cd();
	
	hjpsim_all->Write();
	hpsim_all->Write();
	hjpsim_nm->Write();
	hpsim_nm->Write();
	hjpsim_ftm->Write();
	hpsim_ftm->Write();
	hjpsim_diff->Write();
	hpsim_diff->Write();
			
	out->Save();
	
}
