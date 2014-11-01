//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: F02SteppingAction.cc,v 1.1.1.1 2004/10/25 12:36:47 kuzniak Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F02SteppingAction.hh"
#include "F02SteppingMessenger.hh"
#include "G4ios.hh"
#include <iomanip>
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
//#include "UCNTimedependentField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include <TTree.h>
#include "TMath.h"
#include <TFile.h>
#include "UCNTimeInformation.hh"

#include "G4RunManager.hh"
#include "UCNDetectorConstruction.hh"
#include "UCNFieldSetup.hh"
#include "UCNPhysicsList.hh"
#include "UCNPrimaryGeneratorAction.hh"

bool init=true;
bool gate_on=false;

bool timekey;
double	next=0e9;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

	F02SteppingAction::F02SteppingAction()
: steppingMessenger(0),CurrentStepNum(0)
{
	steppingMessenger = new F02SteppingMessenger(this);

	tree = new TTree("katayama","katayama");
	tree->Branch("x",&x,"x/D");
	tree->Branch("y",&y,"y/D");
	tree->Branch("z",&z,"z/D");
	tree->Branch("px",&px,"px/D");
	tree->Branch("py",&py,"py/D");
	tree->Branch("pz",&pz,"pz/D");
	tree->Branch("t",&t,"t/D");
	tree->Branch("spinx",&spinx,"spinx/D");
	tree->Branch("spiny",&spiny,"spiny/D");
	tree->Branch("spinz",&spinz,"spinz/D");
	tree->Branch("phase",&phase,"phase/D");
	tree->Branch("dph",&dph,"dph/D");
	tree->Branch("evt",&evt,"t/I");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F02SteppingAction::~F02SteppingAction()
{

	CurrentStepNum++;

	delete steppingMessenger ;
	//delete tree;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02SteppingAction::UserSteppingAction(const G4Step* theStep)
{ 

	G4EventManager* evtManager = G4EventManager::GetEventManager();
	G4Event* Evt = (G4Event*)evtManager->GetConstCurrentEvent();
	evt = Evt->GetEventID();

	G4Track * theTrack = theStep->GetTrack();
	// check if it is alive
	if(theTrack->GetTrackStatus()!=fAlive) { return; }
	// check if it is primary
	if(theTrack->GetParentID()!=0) { return; }


	if(false){
		G4double Time = theStep->GetPostStepPoint()->GetGlobalTime();
		G4double trTime = theTrack->GetGlobalTime();
		printf("F02SteppingAction\n");
		printf("theStep->GetPostStepPoint()->GetGlobalTime=%f\n:theTrack->GetGlobalTime=%f\n",Time,trTime);
	}

	timekey=false;

		G4double time;
		UCNTimeInformation* Info = (UCNTimeInformation*)theTrack->GetUserInformation();
		G4bool precise = (Info!=NULL)?Info->GetPrecise():false;
		//get user variables
		if(theTrack->GetUserInformation()!=0&&precise){	
			//time = (double)Info->GetGlobalTime();
			time = theStep->GetPostStepPoint()->GetGlobalTime();
		}else{
			time = theStep->GetPostStepPoint()->GetGlobalTime();
		}

	if(time>next){
		timekey=true;
	}

	if(timekey){

		next+=10e9;

		G4RunManager* runManager = G4RunManager::GetRunManager();
		UCNPrimaryGeneratorAction* upg = (UCNPrimaryGeneratorAction*)runManager->GetUserPrimaryGeneratorAction();
		UCNDetectorConstruction* udc = (UCNDetectorConstruction*)runManager->GetUserDetectorConstruction();
		UCNFieldSetup *ufs = (UCNFieldSetup*)udc->GetFieldSetup();
		G4bool rot = upg->GetRotation();
	  G4int edir = ufs->GetEFieldDirection();

		//if(upg->GetRotation()) printf("Clockwise\n");else printf("Counter-clockwise\n");
			
		//if(ufs->GetEFieldDirection()) printf("EField is up\n");	else printf("EField is down\n");

		G4ThreeVector pos = theStep->GetPostStepPoint()->GetPosition();
		G4ThreeVector mom = theStep->GetPostStepPoint()->GetMomentum();

		G4int stepnr = theTrack->GetCurrentStepNumber();

		double Pol[3];
		double omega = -183.2471999767;
		double omegar = 3000./60.;// 3000 mm/60 mm
		double omegarel = omega*(2.41/1.91)*(3.33e-11/1e-6);
		omegarel = omega*(2.41/1.91)*(theTrack->GetVelocity())/pow((c_light),2.0)*(10.0*1e3*volt/cm)/(1.0*1e-6*tesla);

		double domega;
			if(rot==true/*which is equivalent to clockwise*/){
				domega = -0.5*omegarel*omegarel/(fabs(omega)+omegar);
				//printf("clockwise\n");
			}else{
				domega = -0.5*omegarel*omegarel/(fabs(omega)-omegar);
				//printf("c-clockwise\n");
			}

		//domega=0; //For Non_rel
		omega+=domega;

		//printf("omegarel=%e\n",omegarel);
		//printf("domega=%e\n",domega);
		//printf("omega=%e\n",omega);

		G4ThreeVector C = theStep->GetPostStepPoint()->GetPolarization();
		Pol[0]=C.x(); Pol[1]=C.y(); Pol[2]=C.z();

		dph = TMath::ATan2(TMath::Cos(omega*time/1e9),TMath::Sin(omega*time/1e9)) - TMath::ATan2(C.y(),C.x()) ;
		double dsy = TMath::Cos(omega*time/1e9)-C.y() ;
		//printf("C.x()=%e:Sim(omega*time/1e9)=%e\n",C.x(),TMath::Sin(omega*time/1e9));

		//fill the results
		x=pos.x(); y=pos.y(); z=pos.z();
		px=mom.x(); py=mom.y(); pz=mom.z();
		t=time/1e9;
		spinx=Pol[0]; spiny=Pol[1]; spinz=Pol[2];
		phase=TMath::ATan2(C.y(),C.x());

		tree->Fill();
	}

	////////////////////////////////////////////////
	if(time>100e9){//it correponds to 100 sec
		TFile *file = TFile::Open("check.root","RECREATE");
		tree->Write();
		file->Close();
		G4cout<<"theTrack->GetVelocity()="<<theTrack->GetVelocity()<<G4endl;
		init=true;
		evtManager->AbortCurrentEvent();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

