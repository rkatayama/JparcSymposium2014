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
// $Id: F02EventAction.cc,v 1.1.1.1 2004/10/25 12:36:44 kuzniak Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "F02EventAction.hh"

#include "F02RunAction.hh"

#include "F02EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

//add by katayama
////////////////
#include "UCNPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "UCNDetectorConstruction.hh"
#include "UCNFieldSetup.hh"
#include "UCNPhysicsList.hh"
////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F02EventAction::F02EventAction(F02RunAction* F02RA)
 : eventMessenger(0),
   runaction(F02RA), verboselevel(0),
   drawFlag("all"), printModulo(10000)
{
  eventMessenger = new F02EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F02EventAction::~F02EventAction()
{
	fclose(fp);
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) 
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
     
  if(verboselevel>1)
    G4cout << "<<< Event  " << evtNb << " started." << G4endl;


	//add by katayama
	if(fp==NULL){

		char str[500];

			G4RunManager* runManager = G4RunManager::GetRunManager();
			UCNPrimaryGeneratorAction* upg = (UCNPrimaryGeneratorAction*)runManager->GetUserPrimaryGeneratorAction();
			UCNDetectorConstruction* udc = (UCNDetectorConstruction*)runManager->GetUserDetectorConstruction();
			UCNFieldSetup *ufs = (UCNFieldSetup*)udc->GetFieldSetup();
			UCNPhysicsList* upl = (UCNPhysicsList*)runManager->GetUserPhysicsList();

			G4int rot = (G4int)upg->GetRotation();
			G4int edir = ufs->GetEFieldDirection();
		  G4int textid = upl->GetUCNMaterialBoundary()->GetTextID();

	sprintf(str,"Repository/log/TextID_%d:%s:%s.log",
			        textid,(edir==1)?"Eup":"Edown",
							(rot==1)?"Clockwise":"Counter-Clockwise");

		if( (fp=fopen(str,"w"))==NULL ){
			G4cout<<"File open is failed at event action !"<<G4endl;
			exit(0);
		}//end if ((fp=fopen(...))==NULL)

		fprintf(fp,"TextID_%d:%s:%s\n",textid,
				(edir==1)?"Eup":"Edown",(rot==1)?"Clockwise":"Counter-Clockwise");
		time(&timer);
		fprintf(fp,"Event Start at %s", ctime(&timer) );

	}//end if fp==NULL

	time(&timer);
		printf("Event Start at %s", ctime(&timer) );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::EndOfEventAction(const G4Event* evt)
{  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
   G4int n_trajectories = 0;
   if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();  
   for(G4int i=0; i<n_trajectories; i++) 
      { G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
        if (drawFlag == "all") trj->DrawTrajectory(50);
        else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                               trj->DrawTrajectory(50); 
      }
  }  

 // if(verboselevel>0)
    G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;
		time(&timer);
		printf("Event %d End at %s",evt->GetEventID(),ctime(&timer) );

  // save rndm status
  if (runaction->GetRndmFreq() == 2)
    { 
     CLHEP::HepRandom::saveEngineStatus("endOfEvent.rndm");   
     G4int evtNb = evt->GetEventID();
     if (evtNb%printModulo == 0)
       { 
        G4cout << "\n---> End of Event: " << evtNb << G4endl;
        CLHEP::HepRandom::showEngineStatus();


       }
    }     


	//add by katayama
	if(evt->GetEventID()%10==0){
	time(&timer);
	fprintf(fp,"Event %d End at %s",evt->GetEventID(),ctime(&timer) );
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int F02EventAction::GetEventno()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID() ;
  return evno ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::setEventVerbose(G4int level)
{
  verboselevel = level ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
