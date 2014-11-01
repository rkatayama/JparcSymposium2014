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
// $Id: UCNBetaDecayChannel.cc,v 1.1.1.1 2006/04/26 14:31:54 kuzniak Exp $
// GEANT4 tag $Name:  $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      18 Sep  2001 H.Kurashige
//
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "UCNBetaDecayChannel.hh"
#include "Randomize.hh"
#include "G4RotationMatrix.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"


UCNBetaDecayChannel::UCNBetaDecayChannel(
					 const G4String& theParentName, 
					 G4double        theBR)
  :G4VDecayChannel("UCN Decay"),
   aENuCorr(-0.102)
{
  // set names for daughter particles
  if (theParentName == "ucn") {
    SetBR(theBR);
    SetParent("ucn");
    SetNumberOfDaughters(0);
    //    SetDaughter(0, "e-");
    //    SetDaughter(1, "anti_nu_e");
    //    SetDaughter(2, "proton");
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "UCNBetaDecayChannel:: constructor :";
      G4cout << " parent particle is not muon but ";
      G4cout << theParentName << G4endl;
    }
#endif
  }
}

UCNBetaDecayChannel::~UCNBetaDecayChannel()
{
}

G4DecayProducts *UCNBetaDecayChannel::DecayIt(G4double) 
{
  //  This class describes free neutron beta decay  kinemtics.
  //  This version neglects neutron/electron polarization  
  //  without Coulomb effect

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "UCNBetaDecayChannel::DecayIt ";
#endif

  if (parent == 0) FillParent();  
  //  if (daughters == 0) FillDaughters();
 
  // parent mass
  // G4double parentmass = parent->GetPDGMass();


  //daughters'mass
  //  G4double daughtermass[3]; 
  //  G4double sumofdaughtermass = 0.0;
  //  for (G4int index=0; index<3; index++){
  //    daughtermass[index] = daughters[index]->GetPDGMass();
  //    sumofdaughtermass += daughtermass[index];
  //  }
  //  G4double xmax = parentmass-sumofdaughtermass;  

  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
  
  //  create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  // calculate daughter momentum
  //  G4double daughtermomentum[3];

  // calcurate electron energy
  //  G4double x;                    // Ee
  //  G4double p;                    // Pe
  //  G4double m = daughtermass[0];  //Me
  //  G4double w;                    // cosine of e-nu angle
  //  G4double r;  
  //  G4double r0;
  //  do {
  //      x = xmax*G4UniformRand();
  //      p = sqrt(x*(x+2.0*m));
  //      w = 1.0-2.0*G4UniformRand();
  //      r = p*(x+m)*(xmax-x)*(xmax-x)*(1.0+aENuCorr*p/(x+m)*w);
  //      r0 = G4UniformRand()*(xmax+m)*(xmax+m)*xmax*xmax*(1.0+aENuCorr);
  //  } while (r < r0);    

  //create daughter G4DynamicParticle 
  // rotation materix to lab frame
  /*
  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double theta = acos(costheta)*rad;
  G4double phi  = 2.0*M_PI*G4UniformRand()*rad;
  G4RotationMatrix rm;
  rm.rotateY(theta);
  rm.rotateZ(phi);

  // daughter 0 (electron) in Z direction
  daughtermomentum[0] = p;
  G4ThreeVector direction0(0.0, 0.0, 1.0);
  direction0 = rm * direction0;
  G4DynamicParticle * daughterparticle0 
         = new G4DynamicParticle( daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle0);

  // daughter 1 (nutrino) in XZ plane
  G4double eNu = xmax-x; 
  G4double cosn = w;
  G4double sinn = sqrt((1.0-cosn)*(1.0+cosn));

  G4ThreeVector direction1(sinn, 0.0, cosn);
  direction1 = rm * direction1;
  G4DynamicParticle * daughterparticle1 
         = new G4DynamicParticle( daughters[1], direction1*eNu);
  products->PushProducts(daughterparticle1);

  // daughter 2 (proton) at REST
  G4ThreeVector direction2(0.0, 0.0, 0.0);
  G4DynamicParticle * daughterparticle2 
         = new G4DynamicParticle( daughters[2], direction2);
  products->PushProducts(daughterparticle2);
  */

 // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "UCNBetaDecayChannel::DecayIt ";
    G4cout << "  create decay products in rest frame " <<G4endl;
    //    products->DumpInfo();
  }
#endif
  return products;
}






