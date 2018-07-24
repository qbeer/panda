#ifndef PndTutAnaTask_H
#define PndTutAnaTask_H 1


#include "FairTask.h"
#include <map>
#include <string>
#include "TLorentzVector.h"

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class RhoMassParticleSelector;
class PndAnalysis;
class RhoCandList;
class RhoTuple;


class PndTutAnaTask : public FairTask
{

 public:
	
	// ** Default constructor   
	PndTutAnaTask();
	
	// ** Destructor 
	~PndTutAnaTask();	
	
	// ** Virtual method Init 
	virtual InitStatus Init();
	
	// ** Virtual method Exec 
	virtual void Exec(Option_t* opt);
	
	
	virtual void Finish();

 protected:
	
	
 private: 
	// *** event counter
	int fEvtCount;	
	
	// *** mass selector for the J/psi
	RhoMassParticleSelector *fJpsiMassSel;
	
	// *** a method 
	int  SelectTruePid(PndAnalysis *ana, RhoCandList &l);
	
	// *** declare some histograms
	TH1F *hjpsim_all;
	TH1F *hpsim_all;
	
	TH1F *hjpsim_lpid;
	TH1F *hpsim_lpid;
	
	TH1F *hjpsim_tpid;
	TH1F *hpsim_tpid;
	
	TH1F *hjpsim_trpid;
	TH1F *hpsim_trpid;
	
	TH1F *hjpsim_ftm;
	TH1F *hpsim_ftm;
	
	TH1F *hjpsim_nm;
	TH1F *hpsim_nm;
	
	TH1F *hjpsim_diff;
	TH1F *hpsim_diff;
	
	TH1F *hjpsim_vf;
	TH1F *hjpsim_4cf;
	TH1F *hjpsim_mcf;
	
	TH1F *hjpsi_chi2_vf;
	TH1F *hpsi_chi2_4c;
	TH1F *hjpsi_chi2_mf;
	
	TH1F *hjpsi_prob_vf;
	TH1F *hpsi_prob_4c;
	TH1F *hjpsi_prob_mf;
	
	TH2F *hvpos;
	
	// *** the initial 4-vector
	TLorentzVector fIni;
	
	// *** the PndAnalysis object
	PndAnalysis *fAnalysis;
	
	// *** Get parameter containers
	virtual void SetParContainers();
	
	
	ClassDef(PndTutAnaTask,1);
  
};

#endif
