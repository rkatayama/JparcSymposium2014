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
// $Id: UCNEqMagElectricField.cc,v 1.1.1.1 2006/04/26 14:31:54 kuzniak Exp $
// GEANT4 tag $Name:  $
//
//
//  This is the standard right-hand side for equation of motion.
//
//  The only case another is required is when using a moving reference
//  frame ... or extending the class to include additional Forces,
//  eg an electric field
//
//  10.11.98   V.Grichine
//
// -------------------------------------------------------------------
//
//  4.9.04  adapted for gravity by peter fierlinger  
//
//  3.2013 Ryo Katayama, implemented in relativistic spin precession, based on G4EqEMFieldWithSpin


#include "UCNEqMagElectricField.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "TNtuple.h"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"
#include "TFile.h"

#include <fstream>
#include <stdio.h>

int I=0;

bool texton=true;

UCNEqMagElectricField::UCNEqMagElectricField( G4ElectroMagneticField *emField ): G4EquationOfMotion( emField )
{
	G4cout<<"****************"<<G4endl;
	G4cout<<"volt="<<1*volt<<G4endl;
	G4cout<<"volt/cm="<<1*volt/cm<<G4endl;
	G4cout<<"tesla="<<1*tesla<<G4endl;
	printf("c_light=%0.12lf\n",c_light);;
	G4cout<<"electron_chrge="<<electron_charge<<G4endl;
	G4cout<<"MeV="<<1*MeV<<G4endl;
	G4cout<<"rad="<<1*rad<<G4endl;
	G4cout<<"m="<<1*m<<G4endl;
	G4cout<<"cm="<<1*cm<<G4endl;
	G4cout<<"eV="<<eV<<G4endl;
	G4cout<<"****************"<<G4endl;

	anomaly = -1.9130427 - 1.;
	//For neutron, gfactor = -3.826,085,45
	//anomaly is defined as gfactor/2 - 1.

	omegac = 2*((9.27400915e-24)/(1.602176487e-19))/(1836.152701)*(1.602176487e-19/1.054571628e-34/2.997924258e8)*(rad/m/tesla);
	// According to Jackson's text, 
	// BMT equation for CGS unit system is represented as follows:
	// dS/dt(CGS) = (e/mc) S \cross
	//              [ (g/2-1 +1/\gamma) B
	//               -(g/2-1)\gamma/(\gamma+1) (\beta \cdot B)\beta
	//               -(g/2-\gamma/(\gamma+1) \beta \cross E ]
	//
	// Because Geant4 evaluates the equation of motion using a step size,
	// we used the alternative one such as dS/ds = dS/dt*(1/(c*beta)).
	//
	// Therefore, the equation of motion for Geant4UCN corresponds to fllowing formula:
	// dS/ds(CGS) = ((e/c)/mc)/beta* S \cross
	//              [ (g/2-1 +1/\gamma) B
	//               -(g/2-1)\gamma/(\gamma+1) (\beta \cdot B)\beta
	//               -(g/2-\gamma/(\gamma+1) \beta \cross E ]
	//
	// To transform CGS unit system to MKS unit system,
	// we have to divide electric field by light speed and multiply charge with light speed, formally:
	// E->E/c, charge/c -> charge*c/c = charge.
	//
	// Therefore, dS/ds[cgs] have to be changed as follows:
	// dS/ds(MKS) = (e[MKS]/mc)/beta* S \cross
	//              [ (g/2-1 +1/\gamma) B
	//               -(g/2-1)\gamma/(\gamma+1) (\beta \cdot B)\beta
	//               -(g/2-\gamma/(\gamma+1) \beta \cross E/c ]
	//              = (omegac)/beta* S \cross
	//              [ (g/2-1 +1/\gamma) B
	//               -(g/2-1)\gamma/(\gamma+1) (\beta \cdot B)\beta
	//               -(g/2-\gamma/(\gamma+1) \beta \cross E/c ]
	//
	//              = (omegac)/beta* S \cross
	//              [ (g/2-1 +1/\gamma) B
	//               -(g/2-1)\gamma/(\gamma+1) (\beta \cdot B)\beta
	//               -(g/2-\gamma/(\gamma+1) \beta \cross E ]
	//
	// where omegac is defined as e/mc = elementary charge/particle-mass/light-speed.
	//
	// omegac can be represented by using Bohr-magnetron and the fraction of
	// electron mass to particle mass; For example, assuming particle mass is neutron mass,
	// omegac can be represented as follows:
	//
	// omegac = (elementary charge)/(particle-mass*light-velocity)
	//        = e/(mn*c)
	//        = e/(me*c)*(me/mn)
	//        = e*hbar/(2*me)*(me/mn)*(2/(hbar*c))
	//        = (Bohr magnetron)*(me/mn)*(2/(hbar*c))
	//        = (mu_B)*(me/mn)*(2/(hbar*c))
	//        = 2*(mu_B)/(mn/me)/(hbar*c)
	//        = 2*(mu_B)/(mn/me)*(e/e)/(hbar*c)
	//        = 2*(mu_B)/e/(mn/me)*e/(hbar*c)
	//        = 2*(mu_B)/(mn/me)/(hbar*c)
	//
	// Reffering to the above formula, 
	// the omegac for Geant4 with MKS unit assuming neutron mass can thus be represented as follows:
	// 2.*((9.27400915e-24))/(1836.152701)/(1.054571628e-34*2.997924258e8);
	//
	// At the Geant4, mm and nano-second are used, respectively, for the default units.
	// We have to transform MKS units to Geant4 units properly.
	// It is notable that mu_B has the unit of energy/tesla, and mn/me has no dimmension, and also 
	// hbar*c has the unit equivalent to energy*second*(m/second) = energy*m.
	// We should hence multiply 1/m/tesla with the above formular to transform MKS units to Geant4 units properly.
	//
	// As a result, omegac is finally represented as 
	// omegac = 2.*((9.27400915e-24))/(1836.152701)/(1.054571628e-34*2.997924258e8)*(1/m/tesla);
	//

}

	void  
UCNEqMagElectricField::SetChargeMomentumMass(G4double particleCharge, // e+ units
		G4double,
		G4double particleMass)
{
	fMass = particleMass;
	fCharge = particleCharge;
	fCof_val = particleCharge*eplus*c_light;
}



void
UCNEqMagElectricField::EvaluateRhsGivenB(const G4double y[],
		const G4double G[],
		G4double dydx[] ) const
{


	////////////////////////////////////////////////////////////////////
	//   a neutron in a
	//   gravitational field


	if (fCharge == 0 && fMass > 0){
		//  G4cout << "UCNEX3graveqrhs::evaluaterhsgivenb " << G4endl;
		G4double momentum_mag_square = sqr(y[3]) + sqr(y[4]) + sqr(y[5]);
		G4double inv_momentum_magnitude = 1.0 / sqrt( momentum_mag_square );

		G4double Energy = sqrt(momentum_mag_square + fMass*fMass);
		G4double cof2 = Energy/c_light;
		G4double cof1 = inv_momentum_magnitude*fMass;
		G4double inverse_velocity = Energy*inv_momentum_magnitude/c_light;
		dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
		dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
		dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V
		///// here changed for a field like gravity in all directions

		dydx[3] = -G[6]*cof1*cof2;
		dydx[4] = -G[7]*cof1*cof2;
		dydx[5] = -G[8]*cof1*cof2;   	//  m*g

		// Lab Time of flight
		dydx[7] = inverse_velocity;

		dydx[6] = dydx[8] = 0.;//not used

		G4double beta   = sqrt(momentum_mag_square)/Energy ;
		G4double gamma  = Energy/fMass ;

		G4ThreeVector BField(G[0],G[1],G[2]);
		G4ThreeVector EField(G[3],G[4],G[5]);

		G4ThreeVector u(y[3], y[4], y[5]);
		u *= inv_momentum_magnitude;

		G4double ucb = (anomaly+1./gamma)/beta;
		G4double udb = anomaly*beta*gamma/(1.+gamma) * (BField * u);
		G4double uce = (anomaly + 1./(gamma+1.))/c_light; //MKS
		//G4double uce = anomaly + 1./(gamma+1.); //CGSesu

		G4ThreeVector Spin(y[9],y[10],y[11]);

		G4ThreeVector dSpin 
			= 1.0*omegac*( +ucb*( Spin.cross(BField) )-udb*( Spin.cross(u) )-uce*( u*(Spin*EField) - EField*(Spin*u)) ); //Since 2130530


		dydx[ 9] = dSpin.x();
		dydx[10] = dSpin.y();
		dydx[11] = dSpin.z();
		//G4cout << dSpin.x() << " " << dSpin.y() << " " << dSpin.z() << G4endl;



		if(texton){

			if(I%100000000==0){

				G4double Omega = ucb*omegac/inverse_velocity*G[2];
				G4double Speed = 1./inverse_velocity;
				G4cout<<"**************************"<<G4endl;
				G4cout<<"Omega="<<Omega<<G4endl;
				G4cout<<"Speed="<<Speed<<G4endl;
				printf("B=(%0.2e,%0.2e,%0.2e)\n",G[0],G[1],G[2]);
				printf("E=(%0.2e,%0.2e,%0.2e)\n",G[3],G[4],G[5]);
				printf("g=(%0.2e,%0.2e,%0.2e)\n",G[6],G[7],G[8]);
				G4cout<<"dSpin-1st="<<omegac*ucb*(Spin.cross(BField)).mag()<<G4endl;
				G4cout<<"dSpin-2nd="<<-omegac*udb*( Spin.cross(u).mag() )<<G4endl;
				G4cout<<"dSpin-3rd="<<-omegac*uce*( u*(Spin*EField) - EField*(Spin*u)).mag() <<G4endl;
				G4cout<<"**************************"<<G4endl;

			}

			I++;

		}

		return ;
	}
	return;

}
