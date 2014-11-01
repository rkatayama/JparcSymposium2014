#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(){
	char Eup[200],Edown[200],Clockwise[200],Cclockwise[200];
	strcpy(Eup,"Eup");
	strcpy(Edown,"Edown");
	strcpy(Clockwise,"Clockwise");
	strcpy(Cclockwise,"Counter-Clockwise");

	for(int id=1;id<=9;id++){

		if(id==2)continue; if(id==4)continue;
		if(id==5)continue; if(id==6)continue; if(id==8)continue;

		for(int eupdown=0;eupdown<=1;eupdown++){
			for(int rotation=0;rotation<=1;rotation++){
	
				char str[200];
				char str2[200];

				sprintf(str,"mv TextID=%d:%s:%s.mac",id,(eupdown==0)?Eup:Edown,(rotation==0)?Clockwise:Cclockwise);
				sprintf(str2," TextID=%d:%s:%s.addmac\n",id,(eupdown==0)?Eup:Edown,(rotation==0)?Clockwise:Cclockwise);

				strcat(str,str2);

				printf(str);
				system(str);


			}//end for3
		}//end for2

	}//end for1

}//end main
