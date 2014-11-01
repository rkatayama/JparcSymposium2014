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
// $Id: F02SteppingMessenger.cc,v 1.1.1.1 2004/10/25 12:36:47 kuzniak Exp $
// GEANT4 tag $Name:  $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F02SteppingMessenger.hh"

#include "F02SteppingAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F02SteppingMessenger::F02SteppingMessenger(F02SteppingAction* SA)
:steppingAction (SA)
{

  G4cout << "create messenger for the F02SteppingAction class " << G4endl;

  steppingDir = new G4UIdirectory("/stepping/");
  steppingDir->SetGuidance("stepping control");

  RootNameCmd = new G4UIcmdWithAString("/stepping/rootname",this);
  RootNameCmd->SetGuidance("set rootfile name");
  RootNameCmd->SetGuidance("e.g. /stepping/rootname abc.root");
  RootNameCmd->SetParameterName("fileName",true);
  RootNameCmd->SetDefaultValue("test.root");
  RootNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F02SteppingMessenger::~F02SteppingMessenger()
{
  delete steppingDir;
  delete RootNameCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02SteppingMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{

	if(command == RootNameCmd){
		steppingAction->SetRootFileName(newValue);
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

   
