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

void job(int IDSTART=1,int IDEND=1){


	time_t timer;

	time(&timer);
	chrtime = ctime(&timer);
	printf("chrtime=%s\n",chrtime);

	char STR[200];

	//////////////////////////////

	FILE *FP;

	if(( FP=fopen("log.txt","a") )==NULL){
		printf("ファイルのオープンに失敗しました。\n");
		exit(0);
	}

	fprintf(FP,"===================================\n");

	/////////////////////////

	for(int textid=IDSTART;textid<=IDEND;textid++){

		for(int eup=0;eup<=1;eup++){

			for(int clockwise=0;clockwise<=1;clockwise++){

				char name[2000];
				char nohup[2000];

				//////////////////////////
				sprintf(name,"./Repository/mac/TextID=%d:",textid);

				if(eup) strcat(name,"Eup:");
				else strcat(name,"Edown:");

				if(clockwise) strcat(name,"Clockwise.mac");
				else strcat(name,"Counter-Clockwise.mac");
				//////////////////////////

				//////////////////////////
				sprintf(nohup,"./Repository/nohup/TextID=%d:",textid);

				if(eup) strcat(nohup,"Eup:");
				else strcat(nohup,"Edown:");

				if(clockwise) strcat(nohup,"Clockwise.log");
				else strcat(nohup,"Counter-Clockwise.log");
				//////////////////////////


				sprintf(name,"nohup ./../../geant4_workdir/bin/Linux-g++/JparcProceedings %s >%s &",name,nohup);
				printf("%s\n",name);
				system(name);


				time(&timer);
				chrtime = ctime(&timer);
				fprintf(FP,"Jobs ,%s, are throwed at %s",name,chrtime);

			}//end for

		}//end for

	}//end for


	////////////////////////////////


	fclose(FP);

}//end job
