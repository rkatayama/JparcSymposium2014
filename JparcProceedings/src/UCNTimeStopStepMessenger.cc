/////////////////////////////////////////////////////////////////////////////
//
//   Messenger class for UCN TimeStopStep, 9,2014, Ryo Katayama
//
/////////////////////////////////////////////////////////////////////////////

#include "UCNTimeStopStepMessenger.hh"

#include "UCNTimeStopStep.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNTimeStopStepMessenger::UCNTimeStopStepMessenger(UCNTimeStopStep* SA)
:timeStopStep (SA)
{

  G4cout << "create messenger for the UCNTimeStopStep class " << G4endl;

  stopstepDir = new G4UIdirectory("/stopstep/");
  stopstepDir->SetGuidance("stop step control");

  rootNameCmd = new G4UIcmdWithAString("/stopstep/rootname",this);
  rootNameCmd->SetGuidance("set rootfile name");
  rootNameCmd->SetGuidance("e.g. /stopstep/rootname abc.root");
  rootNameCmd->SetParameterName("fileName",true);
  rootNameCmd->SetDefaultValue("test.root");
  rootNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  stoptimeCmd = new G4UIcmdWithADoubleAndUnit("/stopstep/stoptime",this);
  stoptimeCmd->SetGuidance(" Set stoptime.");
  stoptimeCmd->SetParameterName("stopstep",true);
  stoptimeCmd->SetDefaultValue(1.) ; 
  stoptimeCmd->SetDefaultUnit("s") ; 
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNTimeStopStepMessenger::~UCNTimeStopStepMessenger()
{
  delete stopstepDir;
  delete rootNameCmd;

	delete stoptimeCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UCNTimeStopStepMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{

	if(command == rootNameCmd){
		timeStopStep->SetRootFileName(newValue);
	}

	if(command == stoptimeCmd){
    timeStopStep->SetStopTime(stoptimeCmd->GetNewDoubleValue(newValue));
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

   
