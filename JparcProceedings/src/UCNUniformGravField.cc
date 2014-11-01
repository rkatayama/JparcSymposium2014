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
// $Id: UCNUniformGravField.cc,v 1.2 2008/09/17 09:07:23 kuzniak Exp $
// GEANT4 tag $Name:  $
//
// 30.1.97 V.Grichine
//
// 9.4.04 adapted for gravity by peter fierlinger
// 
//     
//
// -------------------------------------------------------------------

#include "UCNUniformGravField.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4SteppingControl.hh"
#include "G4VParticleChange.hh"

#include "G4ExceptionSeverity.hh"

#include "UCNFieldSetup.hh"//

#include "TVector3.h"

UCNUniformGravField* UCNUniformGravField::theInstance = 0;
G4int UCNUniformGravField::index = 0;

UCNUniformGravField::UCNUniformGravField(const G4ThreeVector FieldVector )
{
	fFieldComponents[2] = FieldVector.z();
	fFieldComponents[0] = FieldVector.x();
	fFieldComponents[1] = FieldVector.y();
	theInstance = this;
}

UCNUniformGravField::UCNUniformGravField(G4double vField,
		G4double vTheta,
		G4double vPhi    )
{
	if(vField >= 0 &&
			vTheta >= 0 && vTheta <= pi &&
			vPhi >= 0 && vPhi <= twopi)
	{
		fFieldComponents[0] = vField*cos(vTheta) ;
		fFieldComponents[1] = vField*cos(vTheta) ;
		fFieldComponents[2] = vField*cos(vTheta) ;
		theInstance = this;
	}
	else
	{
		G4Exception("Invalid parameters in UCNUniformGravField::UCNUniformGravField", "", (G4ExceptionSeverity)0, "") ;
	}
}

UCNUniformGravField::~UCNUniformGravField()
{}

	UCNUniformGravField::UCNUniformGravField (const UCNUniformGravField &p)
: G4ElectricField(p)
{
	for (G4int i=0; i<6; i++)
		fFieldComponents[i] = p.fFieldComponents[i];
	theInstance = this;
}

	UCNUniformGravField&
UCNUniformGravField::operator = (const UCNUniformGravField &p)
{
	for (G4int i=0; i<6; i++)
		fFieldComponents[i] = p.fFieldComponents[i];
	return *this;
}

// ------------------------------------------------------------------------

void UCNUniformGravField::GetFieldValue (const G4double point[4],
		G4double *G) const 
{

	TVector3 v(point[0],point[1],point[2]+3002.77426) ;
	G[0] = 0.027075/pow(v.Mag(),3.)*(v.x()/v.Mag()) ;
	G[1] = 0.027075/pow(v.Mag(),3.)*(v.y()/v.Mag()) ;
	G[2] = 1.0*1e-6*tesla + 0.027075/pow(v.Mag(),3.)*(v.z()/v.Mag() ) ;

	//G[0]= 0. ;							//default
	//G[1]= 0. ;							//default
	//G[2]= 1.e-6*tesla ;

	G[3]= 0. ;							//default
	G[4]= 0. ;							//default
	G[5]= 0. ;
	if(fUCNFieldSetup->GetEFieldDirection()){
		G[5]=10.0*1e3*volt/cm ;
	}else{
		G[5]=-10.0*1e3*volt/cm ;
	}

	if(fUCNFieldSetup->GetGravityEffect()){
	G[6]= fFieldComponents[0] ;
	G[7]= fFieldComponents[1] ;
	G[8]= fFieldComponents[2] ;
	}else{
	G[6]= 0 ;
	G[7]= 0 ;
	G[8]= 0 ;
	}

	return ;

}

UCNUniformGravField* UCNUniformGravField::GetUniformGravField(){
	//G4cout << "timedependentfield: returning my address " << theInstance << ", " << G4endl;
	return theInstance;
}

