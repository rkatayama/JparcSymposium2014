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
// $Id: F02StepCut.cc,v 1.1.1.1 2004/10/25 12:36:47 kuzniak Exp $
// GEANT4 tag $Name:  $
//

#include "F02StepCut.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

F02StepCut::F02StepCut(const G4String& aName)
  : G4VDiscreteProcess(aName),MaxChargedStep(DBL_MAX)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

F02StepCut::~F02StepCut()
{
}

F02StepCut::F02StepCut(F02StepCut& right)
    :G4VDiscreteProcess(right)
{}

void F02StepCut::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
}


