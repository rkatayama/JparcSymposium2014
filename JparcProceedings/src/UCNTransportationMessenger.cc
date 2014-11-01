/////////////////////////////////////////////////////////////////////////////
//
//   Messenger class for UCN Transportation, 9,2014, Ryo Katayama
//
/////////////////////////////////////////////////////////////////////////////

#include "UCNTransportationMessenger.hh"

#include "UCNTransportation.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNTransportationMessenger::UCNTransportationMessenger(UCNTransportation* SA)
:transport (SA)
{

  G4cout << "create messenger for the UCNTransportation class " << G4endl;

  transportDir = new G4UIdirectory("/transport/");
  transportDir->SetGuidance("transport control");

	//add by katayama
  precisetimeCmd = new G4UIcmdWithAnInteger("/transport/precisetime",this);//
  precisetimeCmd->SetGuidance("e.g. /transport/precisetime 0, which means the time variable precision is long double.");//
  //precisetimeCmd->SetGuidance(";otherwise, it is double precision.");//
  precisetimeCmd->SetParameterName("choice",true);//
  precisetimeCmd->SetDefaultValue(0);//
  precisetimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNTransportationMessenger::~UCNTransportationMessenger()
{
  delete transportDir;
	delete precisetimeCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UCNTransportationMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{

	if(command == precisetimeCmd){
    transport->SetPreciseTime(precisetimeCmd->GetNewIntValue(newValue));
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

   
