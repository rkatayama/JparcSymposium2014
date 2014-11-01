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
// $Id: F02SteppingAction.hh,v 1.1.1.1 2004/10/25 12:36:46 kuzniak Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F02SteppingAction_h
#define F02SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

//add by katayama 
#include <TTree.h>

class F02SteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F02SteppingAction : public G4UserSteppingAction
{
  public:
    F02SteppingAction();
   ~F02SteppingAction();

		TTree *tree;

		//add by katayama
		double x,y,z,px,py,pz,t;
    double dph;
		double spinx,spiny,spinz;
		double phase;

    void UserSteppingAction(const G4Step*);
		int evt;


		void SetRootFileName(G4String newValue){
			G4cout<<newValue.c_str()<<G4endl;
			RootFileName = newValue;
		}

		///*
		G4int GetCurrentStepNum(){//add by katayama
			return CurrentStepNum;
		}
		//*/

  private:
    F02SteppingMessenger*    steppingMessenger;
		G4String RootFileName;
		G4int CurrentStepNum;
};

#endif
