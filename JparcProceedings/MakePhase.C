#include <assert.h>
#include <time.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <map>

#include<TGaxis.h>

#include<TROOT.h>
#include<TCut.h>
#include<TFile.h>
#include<TMath.h>
#include<TLegend.h>
#include<THStack.h>
#include<TTreeFormula.h>
#include<TRandom.h>
#include<TCanvas.h>
#include<TStyle.h>
#include<TSystem.h>
#include<TApplication.h>

#include <TRandom.h>/* header file for gRandom */
#include <TCanvas.h>/* header file for TCanvas */
#include <TMath.h>  /* header file for TMath */

#include <TH1.h>    /* header file for 1-d histogram */
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>

#include <TGraph.h>
#include <TGraph2D.h>
#include "TFile.h"
#include "TLegend.h"

#include "TTree.h"
#include "TVector3.h"
#include "TNtuple.h"


//using namespace std;

void MakePhase(){//macro start

	gStyle->SetOptStat(1111111);

	typedef std::map<int,pair<double,pair<double,double> > > TMap;	
	TMap tmap[10];

	TH1D *h[10], *h2[10], *h3[10], *h4[10];
	TH1D *H[10], *H2[10] ;
	double avgavg[10];
	double avgavg2[10];


	double b[10]={0}, w[10]={0} ;

		
	TLegend *leg[9];
	for(int l=0;l<=9;l++){

		char strl1[200], strl2[200], strl3[200], strl4[200];
		sprintf(strl1,"hist-%d",l);
		sprintf(strl2,"hist2-%d",l);
		sprintf(strl3,"hist3-%d",l);
		sprintf(strl4,"hist4-%d",l);

		char Strl1[200], Strl2[200] ;
		sprintf(Strl1,"Hist-%d",l);
		sprintf(Strl2,"Hist2-%d",l);

		printf("Strl1=%s\n",Strl1);
		printf("Strl2=%s\n",Strl2);

		h[l] = new TH1D(strl1,"hist",20,-1e-5,1e-5);
		h2[l] = new TH1D(strl2,"hist2",20,-1e-5,1e-5);
		h3[l] = new TH1D(strl3,"hist3",20,-1e-5,1e-5);
		h4[l] = new TH1D(strl4,"hist4",20,-1e-5,1e-5);

		H[l] = new TH1D(Strl1,"Hist",50,-1e-26,1e-26);
		H2[l] = new TH1D(Strl2,"Hist2",50,-1e-26,1e-26);

		avgavg[l]=0;
		avgavg2[l]=0;

	}


	int id=0;

	for(int s=1;s<=9;s++){

		if(s==2)continue;
		if(s==4)continue;
		if(s==5)continue;
		if(s==6)continue;
		if(s==8)continue;

		for(int k=0;k<=1;k++){//start for-out

		//for(int i=1;i<=4;i++){//start for-out

			char up[2000],down[2000];
			char dir[200];

			strcpy(dir,"./Repository/root/");

			if(k==0){

				printf("---------\n");
				sprintf(up,"%s/TextID_%d:Eup:Counter-Clockwise.root",dir,s);
				sprintf(down,"%s/TextID_%d:Edown:Counter-Clockwise.root",dir,s);
				printf("up=%s\n",up);
				printf("down=%s\n",down);

			}else if(k==1){

				printf("---------\n");
				sprintf(up,"%s/TextID_%d:Eup:Clockwise.root",dir,s);
				sprintf(down,"%s/TextID_%d:Edown:Clockwise.root",dir,s);
				printf("up=%s\n",up);
				printf("down=%s\n",down);

			}

			TFile *file,*file2;

			file = TFile::Open(up);
			file2 = TFile::Open(down);

			TTree *tree = (TTree*)file->Get("katayama");
			TTree *tree2 = (TTree*)file2->Get("katayama");

			double timel, spinx, spiny;
			double timel2, spinx2, spiny2;
			double x, y, z;

			tree->SetBranchAddress("t" ,&timel);
			tree->SetBranchAddress("spinx" ,&spinx);
			tree->SetBranchAddress("spiny" ,&spiny);
			tree2->SetBranchAddress("t" ,&timel2);
			tree2->SetBranchAddress("spinx" ,&spinx2);
			tree2->SetBranchAddress("spiny" ,&spiny2);

			tree->SetBranchAddress("x" ,&x);
			tree->SetBranchAddress("y" ,&y);
			tree->SetBranchAddress("z" ,&z);


			for(int j=0;tree->GetEntries()>j;j++){// for-in star

				id++;

				tree->GetEntry(j);
				tree2->GetEntry(j);

				if(k==0){
					h[s]->Fill(spiny-spiny2);
					h2[s]->Fill(spinx-spinx2);

					double dphi;
					dphi = TMath::ATan2(spiny,spinx) - TMath::ATan2(spiny2,spinx2) ;
					H[s]->Fill(dphi/1e-6*1.64e-28);
					avgavg[s] += dphi/1e-6*1.64e-28/(double)tree->GetEntries();
					//printf("avgavg[%d]=%e\n",s,avgavg[s]);


				}//

				if(k==1){
					h3[s]->Fill(spiny-spiny2);
					h4[s]->Fill(spinx-spinx2);

					double dphi;
					dphi = TMath::ATan2(spiny,spinx) - TMath::ATan2(spiny2,spinx2) ;
					H2[s]->Fill(dphi/1e-6*1.64e-28);
					avgavg2[s] += dphi/1e-6*1.64e-28/(double)tree->GetEntries();
					//printf("avgavg2[%d]=%e\n",s,avgavg2[s]);

				if(j==10)
				printf("i=%d,spiny=%0.15e,time=%0.13lf\n",j,spiny,timel);


				}//


			}//end J for


		//}//end i for

		}//end k for

	}//end s for


	TCanvas *c = new TCanvas("a","a",0,0,800,600);
	gStyle->SetOptStat(1111111);
	c->Divide(2,2);
	//c->Divide(3,3);

	gStyle->SetPalette(1);

	double HH[9],HH2[9];
	double HHH[9],HHH2[9];
	double HHHH[9],HHHH2[9];

	for(int l=1;l<=9;l++ ){

		c->cd(l);

		if(l==2)continue;
		if(l==4)continue;
		if(l==5)continue;
		if(l==6)continue;
		if(l==8)continue;

		if(l==1)c->cd(1);
		if(l==3)c->cd(2);
		if(l==7)c->cd(3);
		if(l==9)c->cd(4);


		gStyle->SetOptStat(0);
		printf("--------------------------\n");
		printf("b[%d]=%e,w[%d]=%e\n",l-1,b[l-1],l-1,w[l-1] );
		H[l]->SetMaximum(100);
		H[l]->SetTitle("");
		H[l]->SetName("");
		H[l]->GetXaxis()->SetTitle("false EDM(e cm)");
		H[l]->GetYaxis()->SetTitle("Events");
		H[l]->Draw("EH");
		printf("H[%d]->GetMean()=%e,H[%d]->GetMeanError()=%e\n",l,H[l]->GetMean(),l,H[l]->GetMeanError() );
		printf("avgavg[%d]=%0.2e\n",l,avgavg[l]);
		((TGaxis*)H[l]->GetYaxis())->SetMaxDigits(5);
		H[l]->SetLineStyle(1);
		char bwstr[200];
		sprintf(bwstr,"b=%0.1f nm,w=%0.1f nm",(l==1||l==3)?0.7:2.1,(l==1||l==7)?22.:45.);
		H[l]->SetTitle(bwstr);
		gStyle->SetTitleFontSize(0.08);
		H[l]->SetLabelSize(0.055);
		H[l]->SetLabelSize(0.05,"y");
		H[l]->GetXaxis()->SetTitleSize(0.05);
		H[l]->GetYaxis()->SetTitleSize(0.055);
		H[l]->GetXaxis()->SetTitleOffset(0.9);
		H[l]->GetYaxis()->SetTitleOffset(0.8);
		H2[l]->SetLineColor(2);
		leg[l] = new TLegend(0.1,0.1,0.3,0.3);
		leg[l]->AddEntry(H[l],"Counter-Clockwise","l");
		leg[l]->AddEntry(H2[l],"Clockwise","l");
		leg[l]->SetTextSize(0.045);
		gStyle->SetStatH(0.25); //高さ
		gStyle->SetStatW(0.2); //幅
		H[l]->SetName("C-Clockwise");
		H2[l]->SetName("Clockwise");
		H2[l]->SetLineStyle(2);
		H2[l]->Draw("samesEH");
		leg[l]->Draw();
		printf("H2[%d]->GetMean()=%e,H2[%d]->GetMeanError()=%e\n",l,H2[l]->GetMean(),l,H2[l]->GetMeanError() );
		printf("avgavg2[%d]=%0.2e\n",l,avgavg2[l]);

		HH[l-1]=H[l]->GetMean(); 
		HH2[l-1]=H2[l]->GetMean();

		HHH[l-1]=H[l]->GetMeanError(); 
		HHH2[l-1]=H2[l]->GetMeanError();

		HHHH[l-1]=H[l]->GetRMS(); 
		HHHH2[l-1]=H2[l]->GetRMS();

	}

	printf("----------average----------\n");
	printf("%0.2e\t",HH[0]);
	printf("%0.2e\t",HH[1]);
	printf("%0.2e\n",HH[2]);

	printf("%0.2e\t",HH[3]);
	printf("%0.2e\t",HH[4]);
	printf("%0.2e\n",HH[5]);

	printf("%0.2e\t",HH[6]);
	printf("%0.2e\t",HH[7]);
	printf("%0.2e\n",HH[8]);
	printf("----------Mean Error----------\n");
	printf("%0.2e\t",HHH[0]);
	printf("%0.2e\t",HHH[1]);
	printf("%0.2e\n",HHH[2]);

	printf("%0.2e\t",HHH[3]);
	printf("%0.2e\t",HHH[4]);
	printf("%0.2e\n",HHH[5]);

	printf("%0.2e\t",HHH[6]);
	printf("%0.2e\t",HHH[7]);
	printf("%0.2e\n",HHH[8]);
	printf("----------RMS----------\n");
	printf("%0.2e\t",HHHH[0]);
	printf("%0.2e\t",HHHH[1]);
	printf("%0.2e\n",HHHH[2]);

	printf("%0.2e\t",HHHH[3]);
	printf("%0.2e\t",HHHH[4]);
	printf("%0.2e\n",HHHH[5]);

	printf("%0.2e\t",HHHH[6]);
	printf("%0.2e\t",HHHH[7]);
	printf("%0.2e\n",HHHH[8]);

	printf("------------------------------\n");

	printf("----------average----------\n");
	printf("%0.2e\t",HH2[0]);
	printf("%0.2e\t",HH2[1]);
	printf("%0.2e\n",HH2[2]);

	printf("%0.2e\t",HH2[3]);
	printf("%0.2e\t",HH2[4]);
	printf("%0.2e\n",HH2[5]);

	printf("%0.2e\t",HH2[6]);
	printf("%0.2e\t",HH2[7]);
	printf("%0.2e\n",HH2[8]);
	printf("----------Mean Error----------\n");
	printf("%0.2e\t",HHH2[0]);
	printf("%0.2e\t",HHH2[1]);
	printf("%0.2e\n",HHH2[2]);

	printf("%0.2e\t",HHH2[3]);
	printf("%0.2e\t",HHH2[4]);
	printf("%0.2e\n",HHH2[5]);

	printf("%0.2e\t",HHH2[6]);
	printf("%0.2e\t",HHH2[7]);
	printf("%0.2e\n",HHH2[8]);
	printf("----------RMS----------\n");
	printf("%0.2e\t",HHHH2[0]);
	printf("%0.2e\t",HHHH2[1]);
	printf("%0.2e\n",HHHH2[2]);

	printf("%0.2e\t",HHHH2[3]);
	printf("%0.2e\t",HHHH2[4]);
	printf("%0.2e\n",HHHH2[5]);

	printf("%0.2e\t",HHHH2[6]);
	printf("%0.2e\t",HHHH2[7]);
	printf("%0.2e\n",HHHH2[8]);

	printf("---------------------\n");

}//End MakePhase
