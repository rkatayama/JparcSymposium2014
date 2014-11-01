#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>

#include<TROOT.h>

#include<TCut.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include<TMultiGraph.h>
#include<TFile.h>
#include<TMath.h>
#include<TLegend.h>
#include<TH1.h>
#include<THStack.h>
#include<TF1.h>
#include<TF2.h>
#include<TTreeFormula.h>
#include<TTree.h>
#include<TRandom.h>
#include<TRandom3.h>
#include<TCanvas.h>
#include<TStyle.h>
#include<TSystem.h>
#include<TVector3.h>
#include<TNtuple.h>

#include<TApplication.h>

#include<time.h>

void makemac(){

	for(int textid=1;textid<=9;textid++){

		for(int eup=0;eup<=1;eup++){

			for(int clockwise=0;clockwise<=1;clockwise++){

				char name[200];
				sprintf(name,"nohup root -l 'SkeltonChange.C(%d,%d,%d)' &",textid,eup,clockwise);
				printf("%s\n",name);
				system(name);

			}//end for

		}//end for

	}//end for

	////////////////////////////////

}//end job
