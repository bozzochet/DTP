#include"TH1F.h"
#include"TCanvas.h"
#include"TTree.h"
#include"TFile.h"
#include"iostream"
#include"TStyle.h"
#include"TF1.h"
#include"TrCluster.hh"
#include<vector>
#include"TClonesArray.h"
#include<iostream>
using namespace std;

int main(int argc,char** argv){

auto inFile = TFile::Open(argv[1]);
TTree *events;
inFile->GetObject("Tree", events);
events->Print();
TClonesArray *a = new TClonesArray("TrCluster", 200);
events->SetBranchAddress("Events", &a);

TCanvas* c=new TCanvas();
TH1F* h=new TH1F("h","h",100,-0.5,0.5);

for (int i=0; i< events->GetEntries(); i++) {
	events->GetEntry(i);
	for (int j=0; j< a->GetEntries(); j++){
		TrCluster *cl = (TrCluster*)a->At(j);
		h->Fill(cl->segm);
	}
}

h->Draw();
c->SaveAs("test.pdf");
}



