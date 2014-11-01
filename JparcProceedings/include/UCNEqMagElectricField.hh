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
// $Id: UCNEqMagElectricField.hh,v 1.1.1.1 2006/04/26 14:31:54 kuzniak Exp $
// GEANT4 tag $Name:  $
//
//
// class G4EqMagElectricField
//
// Class description:
//
// This is the right-hand side of equation of motion in a combined
// electric and magnetic field.

// History:
// - Created. V.Grichine, 10.11.98
// -------------------------------------------------------------------
//
// 4.9.04 adapted for UCN and gravity by peter fierlinger
// 
//  3.2013 Ryo Katayama, implemented in relativistic spin precession, based on G4EqEMFieldWithSpin

#ifndef UCNEQMAGELECTRICFIELD_hh
#define UCNEQMAGELECTRICFIELD_hh

#include "G4EquationOfMotion.hh"
#include "G4ElectroMagneticField.hh"
#include "TNtuple.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"

#include <fstream>
#include <stdio.h>

class UCNEqMagElectricField : public G4EquationOfMotion
{
public:  // with description
  
  UCNEqMagElectricField(G4ElectroMagneticField *emField );
  //: G4EquationOfMotion( emField ) {;}
  
  ~UCNEqMagElectricField() {;} 
  
  void  SetChargeMomentumMass(G4double particleCharge, // in e+ units
			      G4double MomentumXc,
			      G4double mass);
  
  void EvaluateRhsGivenB(const G4double y[],
			 const G4double Field[],
			 G4double dydx[] ) const;
  // Given the value of the electromagnetic field, this function 
  // calculates the value of the derivative dydx.
  
private:
  
  //G4double        fElectroMagCof ;
  //G4double        fMassCof;
  G4double  fMass;
  G4double  fCharge;
  G4double  fCof_val;
  G4double  anomaly;
  G4double  omegac;
  G4double  omegac2;
  G4double  magN;
  
public:
  
  G4int  num;
  int Number;
  void NumPlus(G4int Num){num = num+1;};
};

#endif
