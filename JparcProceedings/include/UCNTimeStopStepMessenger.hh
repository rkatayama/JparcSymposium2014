/////////////////////////////////////////////////////////////////////////////
//
//   Messenger class for UCN TimeStopStep, 9,2014, Ryo Katayama
//
/////////////////////////////////////////////////////////////////////////////

#ifndef UCNTimeStopStepMessenger_h
#define UCNTimeStopStepMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class UCNTimeStopStep;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class UCNTimeStopStepMessenger: public G4UImessenger
{
  public:

   UCNTimeStopStepMessenger(UCNTimeStopStep* );
  ~UCNTimeStopStepMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   UCNTimeStopStep* timeStopStep;

   G4UIdirectory*     stopstepDir;

   G4UIcmdWithAString* rootNameCmd;

   G4UIcmdWithADoubleAndUnit* stoptimeCmd;

};

#endif

