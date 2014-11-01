/////////////////////////////////////////////////////////////////////////////
//
//   Messenger class for UCN shutters, 9.9.04, peter fierlinger
//
/////////////////////////////////////////////////////////////////////////////


#include "UCNMaterialBoundaryMessenger.hh"
#include "UCNMaterialBoundary.hh"
#include "G4UIdirectory.hh"
#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNMaterialBoundaryMessenger::UCNMaterialBoundaryMessenger(UCNMaterialBoundary* SA)
:tPropagator (SA)
{
  
G4cout << "create messenger for the UCN material boundary class " << G4endl;
  materialboundaryDir = new G4UIdirectory("/materialboundary/");
  materialboundaryDir->SetGuidance("savetofile parameters");

  setVerboseCmd = new G4UIcmdWithAnInteger("/materialboundary/setverbose",this);
  setVerboseCmd->SetGuidance("Set verbose level ." );
  setVerboseCmd->SetParameterName("level",true);
  setVerboseCmd->SetDefaultValue(0);

	//add by katayama
  setRoughnessCmd = new G4UIcmdWithAString("/materialboundary/setroughness",this);
  setRoughnessCmd->SetGuidance("Choose surface roughness.");
  setRoughnessCmd->SetGuidance("e.g. /materialboundary/setroughness model textid");
  setRoughnessCmd->SetGuidance("where model must be microroughness or lambert.");
  setRoughnessCmd->SetParameterName("choice",true);
  setRoughnessCmd->SetDefaultValue("");
  setRoughnessCmd->AvailableForStates(G4State_Idle);
	///////////////////////


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNMaterialBoundaryMessenger::~UCNMaterialBoundaryMessenger()
{

  delete setVerboseCmd;
	//add by katayama
	delete setRoughnessCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UCNMaterialBoundaryMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

 if( command == setVerboseCmd )    {    // set shutteropen
   tPropagator->setVerbose(setVerboseCmd->GetNewIntValue(newValue));
 }

	//add by katayama
 if( command == setRoughnessCmd )    {    
   tPropagator->AssignRoughness(newValue);
 }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


