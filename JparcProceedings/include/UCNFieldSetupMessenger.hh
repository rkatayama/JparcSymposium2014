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
// $Id: UCNFieldSetupMessenger.hh,v 1.1.1.1 2004/10/25 12:36:46 kuzniak Exp $
// GEANT4 tag $Name:  $
//
//  modified for UCN, 9.9.04 peter fierlinger

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef UCNFieldSetupMessenger_h
#define UCNFieldSetupMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class UCNFieldSetup;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;


class UCNFieldSetupMessenger: public G4UImessenger
{
  public:
    UCNFieldSetupMessenger(UCNFieldSetup* );
   ~UCNFieldSetupMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    void SetNewValue(G4UIcommand*, G4int);
    
  private:

    UCNFieldSetup*     fEFieldSetup;
    
    G4UIdirectory*             UCNdetDir;

    G4UIcmdWithAnInteger*      StepperCmd;
    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithADoubleAndUnit* MinStepCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;

    G4UIcmdWithAString*        AbsMaterCmd;

  private:
    G4UIcmdWithAnInteger* EdirectionCmd;
    G4UIcmdWithAnInteger* GravityEffectCmd;

};

#endif

