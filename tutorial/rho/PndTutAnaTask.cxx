// ************************************************************************
//
// psi(2S) -> J/psi (-> mu+ mu-) pi+ pi- Analysis Example Task
// 
// for the Rho Tutorial, see
// http://panda-wiki.gsi.de/cgi-bin/viewauth/Computing/PandaRootRhoTutorial
//
// K.Goetzen 7/2013
// ************************************************************************


// The header file
#include "PndTutAnaTask.h"

// C++ headers
#include <string>
#include <iostream>

// FAIR headers
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"

// ROOT headers
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"

// RHO headers
#include "RhoCandidate.h"
#include "RhoHistogram/RhoTuple.h"
#include "RhoFactory.h"
#include "RhoMassParticleSelector.h"

// Analysis headers
#include "PndAnalysis.h"
#include "Pnd4CFitter.h"
#include "PndKinVtxFitter.h"
#include "PndKinFitter.h"
#include "PndVtxPoca.h"
		
		
using std::cout;
using std::endl;


// -----   Default constructor   -------------------------------------------
PndTutAnaTask::PndTutAnaTask() :
  FairTask("Panda Tutorial Analysis Task") { 
}
// -------------------------------------------------------------------------


// -----   Destructor   ----------------------------------------------------
PndTutAnaTask::~PndTutAnaTask() { }
// -------------------------------------------------------------------------


// ----- Method to select true PID candidates
int PndTutAnaTask::SelectTruePid(PndAnalysis *ana, RhoCandList &l)
{
	int removed = 0;
	
	for (int ii=l.GetLength()-1;ii>=0;--ii)
	{
		if ( !(ana->McTruthMatch(l[ii])) )
		{
			l.Remove(l[ii]);
			removed++;
		}
	}
	
	return removed;
}
// -------------------------------------------------------------------------


// -----   Public method Init   --------------------------------------------
InitStatus PndTutAnaTask::Init() 
{		
	// initialize analysis object
	fAnalysis = new PndAnalysis();
	
	// reset the event counter
	fEvtCount = 0;

	// Mass selector for the jpsi cands
	fJpsiMassSel=new RhoMassParticleSelector("jpsi",3.096,1.0);
	
	// create the histograms
	hjpsim_all = new TH1F("hjpsim_all","J/#psi mass (all)",200,0,4.5);
	hpsim_all  = new TH1F("hpsim_all","#psi(2S) mass (all)",200,0,5);
	
	hjpsim_lpid = new TH1F("hjpsim_lpid","J/#psi mass (loose pid)",200,0,4.5);
	hpsim_lpid  = new TH1F("hpsim_lpid","#psi(2S) mass (loose pid)",200,0,5);
	
	hjpsim_tpid = new TH1F("hjpsim_tpid","J/#psi mass (tight pid)",200,0,4.5);
	hpsim_tpid  = new TH1F("hpsim_tpid","#psi(2S) mass (tight pid)",200,0,5);
	
	hjpsim_trpid = new TH1F("hjpsim_trpid","J/#psi mass (true pid)",200,0,4.5);
	hpsim_trpid  = new TH1F("hpsim_trpid","#psi(2S) mass (true pid)",200,0,5);
	
	hjpsim_ftm = new TH1F("hjpsim_ftm","J/#psi mass (full truth match)",200,0,4.5);
	hpsim_ftm  = new TH1F("hpsim_ftm","#psi(2S) mass (full truth match)",200,0,5);
	
	hjpsim_nm = new TH1F("hjpsim_nm","J/#psi mass (no truth match)",200,0,4.5);
	hpsim_nm  = new TH1F("hpsim_nm","#psi(2S) mass (no truth match)",200,0,5);
	
	hjpsim_diff = new TH1F("hjpsim_diff","J/#psi mass diff to truth",100,-2,2);
	hpsim_diff  = new TH1F("hpsim_diff","#psi(2S) mass diff to truth",100,-2,2);
	
	hjpsim_vf   = new TH1F("hjpsim_vf","J/#psi mass (vertex fit)",200,0,4.5);
	hjpsim_4cf  = new TH1F("hjpsim_4cf","J/#psi mass (4C fit)",200,0,4.5);
	hjpsim_mcf  = new TH1F("hjpsim_mcf","J/#psi mass (mass constraint fit)",200,0,4.5);
	
	hjpsi_chi2_vf  = new TH1F("hjpsi_chi2_vf", "J/#psi: #chi^{2} vertex fit",100,0,10);
	hpsi_chi2_4c   = new TH1F("hpsi_chi2_4c",  "#psi(2S): #chi^{2} 4C fit",100,0,250);
	hjpsi_chi2_mf  = new TH1F("hjpsi_chi2_mf", "J/#psi: #chi^{2} mass fit",100,0,10);

	hjpsi_prob_vf  = new TH1F("hjpsi_prob_vf", "J/#psi: Prob vertex fit",100,0,1);
	hpsi_prob_4c   = new TH1F("hpsi_prob_4c",  "#psi(2S): Prob 4C fit",100,0,1);
	hjpsi_prob_mf  = new TH1F("hjpsi_prob_mf", "J/#psi: Prob mass fit",100,0,1);
	
	hvpos = new TH2F("hvpos","(x,y) projection of fitted decay vertex",100,-2,2,100,-2,2);
	
	
	// *** the lorentz vector of the initial psi(2S)
	fIni.SetXYZT(0, 0, 6.231552, 7.240065);
	
	return kSUCCESS;
}

// -------------------------------------------------------------------------
	
void PndTutAnaTask::SetParContainers() 
{
  // Get run and runtime database
  FairRun* run = FairRun::Instance();
  if ( ! run ) Fatal("SetParContainers", "No analysis run");
}

// -------------------------------------------------------------------------


// -----   Public method Exec   --------------------------------------------
void PndTutAnaTask::Exec(Option_t*)
{
	// *** some variables
	int i=0,j=0, k=0, l=0;
	
	// necessary to read the next event
	fAnalysis->GetEventInTask();
	
	if (!(++fEvtCount%100)) cout << "evt "<<fEvtCount<<endl;
	
	// *** RhoCandLists for the analysis
	RhoCandList muplus, muminus, piplus, piminus, jpsi, psi2s;
	
	
	// *** Select with no PID info ('All'); type and mass are set 		
	fAnalysis->FillList(muplus,  "MuonAllPlus");
	fAnalysis->FillList(muminus, "MuonAllMinus");
	fAnalysis->FillList(piplus,  "PionAllPlus");
	fAnalysis->FillList(piminus, "PionAllMinus");
		
	// *** combinatorics for J/psi -> mu+ mu-
	jpsi.Combine(muplus, muminus);
		
		
	// ***
	// *** do the TRUTH MATCH for jpsi
	// ***
	jpsi.SetType(443);
			
	for (j=0;j<jpsi.GetLength();++j) 
	{
		hjpsim_all->Fill( jpsi[j]->M() );
		
		if (fAnalysis->McTruthMatch(jpsi[j]))
		{ 
			hjpsim_ftm->Fill( jpsi[j]->M() );
		 	hjpsim_diff->Fill( jpsi[j]->GetMcTruth()->M() - jpsi[j]->M() );
		}
		else 
			hjpsim_nm->Fill( jpsi[j]->M() );
	}
	
	// ***
	// *** do VERTEX FIT (J/psi)
	// ***
	for (j=0;j<jpsi.GetLength();++j) 
	{
		PndKinVtxFitter vtxfitter(jpsi[j]);	// instantiate a vertex fitter
		vtxfitter.Fit();
		
		double chi2_vtx = vtxfitter.GetChi2();	// access chi2 of fit
		double prob_vtx = vtxfitter.GetProb();	// access probability of fit
		hjpsi_chi2_vf->Fill(chi2_vtx);
		hjpsi_prob_vf->Fill(prob_vtx);			
			
		if ( prob_vtx > 0.01 )				// when good enough, fill some histos
		{
			RhoCandidate *jfit = jpsi[j]->GetFit();	// access the fitted cand
			TVector3 jVtx=jfit->Pos();		// and the decay vertex position
			
			hjpsim_vf->Fill(jfit->M());            
			hvpos->Fill(jVtx.X(),jVtx.Y());
		}
	}
	
	// *** some rough mass selection
	jpsi.Select(fJpsiMassSel);
	
	// *** combinatorics for psi(2S) -> J/psi pi+ pi-
	psi2s.Combine(jpsi, piplus, piminus);
	
	
	// ***
	// *** do the TRUTH MATCH for psi(2S)
	// ***
	psi2s.SetType(88888);
	for (j=0;j<psi2s.GetLength();++j) 
	{
		hpsim_all->Fill( psi2s[j]->M() );
		
		if (fAnalysis->McTruthMatch(psi2s[j])) 
		{
		 	hpsim_ftm->Fill( psi2s[j]->M() );
		 	hpsim_diff->Fill( psi2s[j]->GetMcTruth()->M() - psi2s[j]->M() );
		}
		else 
			hpsim_nm->Fill( psi2s[j]->M() );
	}			
		
	// ***
	// *** do 4C FIT (initial psi(2S) system)
	// ***
	for (j=0;j<psi2s.GetLength();++j) 
	{
		PndKinFitter fitter(psi2s[j]);	// instantiate the kin fitter in psi(2S)
		fitter.Add4MomConstraint(fIni);	// set 4 constraint
		fitter.Fit();		            // do fit
			
		double chi2_4c = fitter.GetChi2();	// get chi2 of fit
		double prob_4c = fitter.GetProb();	// access probability of fit
		hpsi_chi2_4c->Fill(chi2_4c);
		hpsi_prob_4c->Fill(prob_4c);			
			
		if ( prob_4c > 0.01 )			// when good enough, fill some histo
		{
			RhoCandidate *jfit = psi2s[j]->Daughter(0)->GetFit();	// get fitted J/psi
			
			hjpsim_4cf->Fill(jfit->M());
		}
	}		
	
	
	// ***
	// *** do MASS CONSTRAINT FIT (J/psi)
	// ***
	for (j=0;j<jpsi.GetLength();++j) 
	{
		PndKinFitter mfitter(jpsi[j]);		// instantiate the PndKinFitter in psi(2S)
		mfitter.AddMassConstraint(3.0965);	// add the mass constraint
		mfitter.Fit();				// do fit
		
		double chi2_m = mfitter.GetChi2();	// get chi2 of fit
		double prob_m = mfitter.GetProb();	// access probability of fit
		hjpsi_chi2_mf->Fill(chi2_m);
		hjpsi_prob_mf->Fill(prob_m);			
			
		if ( prob_m > 0.01 )				// when good enough, fill some histo
		{
			RhoCandidate *jfit = jpsi[j]->GetFit();	// access the fitted cand
			hjpsim_mcf->Fill(jfit->M());
		}
	}		
	
	
	// ***
	// *** TRUE PID combinatorics
	// ***
	
	// *** do MC truth match for PID type
	SelectTruePid(fAnalysis, muplus);
	SelectTruePid(fAnalysis, muminus);
	SelectTruePid(fAnalysis, piplus);
	SelectTruePid(fAnalysis, piminus);
			
	// *** all combinatorics again with true PID
	jpsi.Combine(muplus, muminus);
	for (j=0;j<jpsi.GetLength();++j) hjpsim_trpid->Fill( jpsi[j]->M() );
	jpsi.Select(fJpsiMassSel);
	
	psi2s.Combine(jpsi, piplus, piminus);
	for (j=0;j<psi2s.GetLength();++j) hpsim_trpid->Fill( psi2s[j]->M() );
	
	
	// ***
	// *** LOOSE PID combinatorics
	// ***
	
	// *** and again with PidAlgoMvd;PidAlgoStt;PidAlgoDrc and loose selection
	fAnalysis->FillList(muplus,  "MuonLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
	fAnalysis->FillList(muminus, "MuonLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
	fAnalysis->FillList(piplus,  "PionLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
	fAnalysis->FillList(piminus, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
	
	jpsi.Combine(muplus, muminus);
	for (j=0;j<jpsi.GetLength();++j) hjpsim_lpid->Fill( jpsi[j]->M() );
	jpsi.Select(fJpsiMassSel);
	
	psi2s.Combine(jpsi, piplus, piminus);
	for (j=0;j<psi2s.GetLength();++j) hpsim_lpid->Fill( psi2s[j]->M() );
	
	
	// ***
	// *** TIGHT PID combinatorics
	// ***
	
	// *** and again with PidAlgoMvd;PidAlgoStt and tight selection
	fAnalysis->FillList(muplus,  "MuonTightPlus",  "PidAlgoMdtHardCuts");
	fAnalysis->FillList(muminus, "MuonTightMinus", "PidAlgoMdtHardCuts");
	fAnalysis->FillList(piplus,  "PionLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
	fAnalysis->FillList(piminus, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
	
	jpsi.Combine(muplus, muminus);
	for (j=0;j<jpsi.GetLength();++j) hjpsim_tpid->Fill( jpsi[j]->M() );
	jpsi.Select(fJpsiMassSel);
	
	psi2s.Combine(jpsi, piplus, piminus);
	for (j=0;j<psi2s.GetLength();++j) hpsim_tpid->Fill( psi2s[j]->M() );
	
}


void PndTutAnaTask::Finish()
{	
	hjpsim_all->Write();
	hpsim_all->Write();
	hjpsim_lpid->Write();
	hpsim_lpid->Write();
	hjpsim_tpid->Write();
	hpsim_tpid->Write();
	hjpsim_trpid->Write();
	hpsim_trpid->Write();
	
	hjpsim_ftm->Write();
	hpsim_ftm->Write();
	hjpsim_nm->Write();
	hpsim_nm->Write();
	
	hpsim_diff->Write();
	hjpsim_diff->Write();
	
	hjpsim_vf->Write();
	hjpsim_4cf->Write();
	hjpsim_mcf->Write();
	
	hjpsi_chi2_vf->Write();
	hpsi_chi2_4c->Write();
	hjpsi_chi2_mf->Write();
			
	hjpsi_prob_vf->Write();
	hpsi_prob_4c->Write();
	hjpsi_prob_mf->Write();
			
	hvpos->Write();
		
}

ClassImp(PndTutAnaTask)
