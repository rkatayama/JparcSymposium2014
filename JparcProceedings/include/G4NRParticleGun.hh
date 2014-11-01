//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// Original file:
// $Id: G4NRParticleGun.hh,v 1.1 2008/09/17 08:56:03 kuzniak Exp $
// GEANT4 tag $Name:  $
// adapted by Stefan Heule for non-relativistic particles as ultracold neutrons
//

#ifndef G4NRParticleGun_h
#define G4NRParticleGun_h 1


#include "globals.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleMomentum.hh"

class G4Event;
class G4NRParticleGunMessenger;

// class description:
//
//  This is a concrete class of G4VPrimaryGenerator. It shoots a particle of given type
// to a given direction with a given kinetic energy. 
//  The position and time of the primary particle must be set by the corresponding
// set methods of G4VPrimaryGenerator base class, otherwise zero will be set.
//
//  The FAQ to this class is for randomizing position/direction/kinetic energy of primary
// particle. But, G4NRParticleGun does NOT have any way of randomization. Instead, the user's
// concrete implementation of G4VUserPrimaryGeneratorAction which transmits G4Event object
// to this particle gun can randomize these quantities and set to this particle gun before
// invoking GeneratePrimaryVertex() method.
//  Note that, even if the particle gun shoots more than one particles at one invokation of
// GeneratePrimaryVertex() method, all particles have the same physical quantities. If the
// user wants to shoot two particles with different momentum, position, etc., invoke
// GeneratePrimaryVertex() method twice and set quantities on demand to the particle gun.
//

class G4NRParticleGun:public G4VPrimaryGenerator
{
  public: // with description
     G4NRParticleGun();
     G4NRParticleGun(G4int numberofparticles);
     G4NRParticleGun(G4ParticleDefinition * particleDef, 
                   G4int numberofparticles = 1);
     // costructors. "numberofparticles" is number of particles to be shoot at one invokation
     // of GeneratePrimaryVertex() method. All paricles are shoot with the same physical
     // quantities.

  public:
     virtual ~G4NRParticleGun();
     G4NRParticleGun(const G4NRParticleGun &right);

     const G4NRParticleGun & operator=(const G4NRParticleGun &right);
     G4int operator==(const G4NRParticleGun &right) const;
     G4int operator!=(const G4NRParticleGun &right) const;

  public: // with description
     virtual void GeneratePrimaryVertex(G4Event* evt);
     // Creates a primary vertex at the given point and put primary particles to it.
     // Followings are set methods for the particle properties.
     //   SetParticleDefinition should be called first.  
     //   By using SetParticleMomentum(), both particle_momentum_direction and
     //   particle_energy(Kinetic Energy) are set.
     //   
     void SetParticleDefinition
       (G4ParticleDefinition * aParticleDefinition);
     void SetParticleMomentum(G4ParticleMomentum aMomentum);
     inline void SetParticleMomentumDirection
                 (G4ParticleMomentum aMomentumDirection)
     { particle_momentum_direction =  aMomentumDirection.unit(); }
     inline void SetParticleEnergy(G4double aKineticEnergy)
     { particle_energy = aKineticEnergy; }
     inline void SetParticleCharge(G4double aCharge)
     { particle_charge = aCharge; }
     inline void SetParticlePolarization(G4ThreeVector aVal)
     { particle_polarization = aVal; }
     inline void SetNumberOfParticles(G4int i)
     { NumberOfParticlesToBeGenerated = i; }

  public:
     inline G4ParticleDefinition* GetParticleDefinition()
     { return particle_definition; }
     inline G4ParticleMomentum GetParticleMomentumDirection()
     { return particle_momentum_direction; }
     inline G4double GetParticleEnergy()
     { return particle_energy; }
     inline G4double GetParticleCharge()
     { return particle_charge; }
     inline G4ThreeVector GetParticlePolarization()
     { return particle_polarization; }
     inline G4int GetNumberOfParticles()
     { return NumberOfParticlesToBeGenerated; }

  protected:  
     virtual void SetInitialValues();

     G4int                 NumberOfParticlesToBeGenerated;
     G4ParticleDefinition* particle_definition;
     G4ParticleMomentum    particle_momentum_direction;
     G4double	           particle_energy;
     G4double	           particle_charge;
     G4ThreeVector         particle_polarization;

  private:
     G4NRParticleGunMessenger* theMessenger;
};

#endif







