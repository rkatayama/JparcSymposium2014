/////////////////////////////////////////////////////////////////////////////
//
//   Messenger class for UCN Transportation, 9,2014, Ryo Katayama
//
/////////////////////////////////////////////////////////////////////////////

#ifndef UCNTransportationMessenger_h
#define UCNTransportationMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class UCNTransportation;
class G4UIdirectory;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class UCNTransportationMessenger: public G4UImessenger
{
  public:

   UCNTransportationMessenger( UCNTransportation* );
  ~UCNTransportationMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   UCNTransportation* transport;

   G4UIdirectory*     transportDir;

   G4UIcmdWithAnInteger* precisetimeCmd;

};

#endif

