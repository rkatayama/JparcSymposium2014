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

void SkeltonChange(int TextID=0,int Eup=0,int Clockwise=0){

	char name[200];
	sprintf(name,"TextID=%d:",TextID);

	if(Eup) strcat(name,"Eup:");
	else strcat(name,"Edown:");

	if(Clockwise) strcat(name,"Clockwise.addmac");
	else strcat(name,"Counter-Clockwise.addmac");

	FILE *FP,*FP2;

	if(( FP=fopen("skelton.mac","r") )==NULL){
		printf("ファイルのオープンに失敗しました。\n");
	}
	else{

	if(( FP2=fopen(name,"w") )==NULL){
		printf("ファイルのオープンに失敗しました。\n");
	}//end if==NULL
	else{
	char STR[512];

		while(fgets(STR,512,FP)!=NULL){


			//printf("------\n");
			//printf("STR=%s",STR);
			char str1[200],str2[200];
			int aaa=0;
			char STR2[512];

		  sscanf(STR,"%s %s %d",str1,str2,&aaa);

			if(!strcmp(STR,"\n")){
				;
			}else if(!strcmp("/gun/rotation",str1)){

			//	printf("1111111111\n");
				//printf("%s %d\n",str1,Clockwise);
				sprintf(STR2,"%s %d\n",str1,Clockwise);
				fprintf(FP2,"%s",STR2);

			}else if(!strcmp("/field/Edirection",str1)){

			//	printf("2222222222\n");
				//printf("%s %d\n",str1,Eup);
				sprintf(STR2,"%s %d\n",str1,Eup);
				fprintf(FP2,"%s",STR2);

			}else if(!strcmp("/materialboundary/setroughness",str1)){

			//	printf("333333333333\n");
				//printf("%s %s %d\n",str1,str2,TextID);
				sprintf(STR2,"%s %s %d\n",str1,str2,TextID);
				fprintf(FP2,"%s",STR2);

			}else{

			//	printf("44444444444\n");
				fprintf(FP2,"%s",STR);

			}

			//printf("------\n");
		
		}//end while

		}//end else

	fclose(FP2);

	}//end else

	fclose(FP);

}//end job
