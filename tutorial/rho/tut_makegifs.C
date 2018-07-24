tut_makegifs(TString fn="signal_ana.root", TString type="gif")
{
  TFile *f=new TFile(fn,"READ");
    
  TKey *key;
  
  TCanvas *c1=new TCanvas("c1","c1",10,10,700,600);
  c1->cd();
  
  int i=1;
  
  TIter next(f->GetListOfKeys());
  
  while (key = (TKey*)next())
  {
	TObject *obj = key->ReadObj();
	if (!obj->InheritsFrom("TH1")) continue;
	TString name=TString(obj->GetName())+"."+type;
	cout <<"Creating "<<name<<endl;
	TH1* h=(TH1*) obj;
	h->Draw();
	c1->SaveAs(name);
  }
  //f->Close();
}
