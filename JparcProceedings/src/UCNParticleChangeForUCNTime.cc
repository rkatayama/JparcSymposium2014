/////////////////////////////////////////////////////////////////////////////
//
//   Implemented for the scheme, based on G4ParticleChangeForTransport 10,2014 Ryo Katayama
//
/////////////////////////////////////////////////////////////////////////////

// The class can make time variable updated with long double precision just after AlongStepDoIt.

#include "UCNParticleChangeForUCNTime.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"

#include "UCNTimeInformation.hh"

UCNParticleChangeForUCNTime::UCNParticleChangeForUCNTime()
  : G4ParticleChange(), isMomentumChanged(false), theMaterialChange(0),
    theMaterialCutsCoupleChange(0), theSensitiveDetectorChange(0),
    fpVectorOfAuxiliaryPointsPointer(0),printverbose(false)
{
  if (verboseLevel>2) {
    G4cout << "UCNParticleChangeForUCNTime::UCNParticleChangeForUCNTime() "
           << G4endl;
  }
}

UCNParticleChangeForUCNTime::~UCNParticleChangeForUCNTime() 
{
  if (verboseLevel>2) {
    G4cout << "UCNParticleChangeForUCNTime::~UCNParticleChangeForUCNTime() "
           << G4endl;
  }
}

UCNParticleChangeForUCNTime::
UCNParticleChangeForUCNTime(const UCNParticleChangeForUCNTime &r)
  : G4ParticleChange(r),
    fpVectorOfAuxiliaryPointsPointer(0)
{
  if (verboseLevel>0) {
    G4cout << "UCNParticleChangeForUCNTime::  copy constructor is called "
           << G4endl;
  }
  theTouchableHandle = r.theTouchableHandle;
  isMomentumChanged = r.isMomentumChanged;
  theMaterialChange = r.theMaterialChange;
  theMaterialCutsCoupleChange = r.theMaterialCutsCoupleChange;
  theSensitiveDetectorChange = r.theSensitiveDetectorChange;
}

// assignemnt operator
UCNParticleChangeForUCNTime &
UCNParticleChangeForUCNTime::operator=(const UCNParticleChangeForUCNTime &r)
{
   if (verboseLevel>1) {
    G4cout << "UCNParticleChangeForUCNTime:: assignment operator is called "
           << G4endl;
   }
   if (this != &r)
   {
      theListOfSecondaries = r.theListOfSecondaries;
      theSizeOftheListOfSecondaries = r.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = r.theNumberOfSecondaries;
      theStatusChange = r.theStatusChange;
      theTouchableHandle = r.theTouchableHandle;
      theMaterialChange = r.theMaterialChange;
      theMaterialCutsCoupleChange = r.theMaterialCutsCoupleChange;
      theSensitiveDetectorChange = r.theSensitiveDetectorChange;
      theMomentumDirectionChange = r.theMomentumDirectionChange;
      thePolarizationChange = r.thePolarizationChange;
      thePositionChange = r.thePositionChange;
      theTimeChange = r.theTimeChange;
      theEnergyChange = r.theEnergyChange;
      theVelocityChange        = r.theVelocityChange;
      theTrueStepLength = r.theTrueStepLength;
      theLocalEnergyDeposit = r.theLocalEnergyDeposit;
      theSteppingControlFlag = r.theSteppingControlFlag;
   }
   return *this;
}

//----------------------------------------------------------------
// methods for updating G4Step
//

G4Step* UCNParticleChangeForUCNTime::UpdateStepForAtRest(G4Step* pStep)
{
  // Nothing happens for AtRestDoIt
  if (verboseLevel>0) {
    G4cout << "UCNParticleChangeForUCNTime::UpdateStepForAtRest() is called"
           << G4endl;
    G4cout << " Nothing happens for this method " << G4endl;
  }
  //  Update the G4Step specific attributes
  return UpdateStepInfo(pStep);
}


G4Step* UCNParticleChangeForUCNTime::UpdateStepForAlongStep(G4Step* pStep)
{
  // Smooth curved tajectory representation: let the Step know about
  // the auxiliary trajectory points (jacek 30/10/2002)
  pStep->SetPointerToVectorOfAuxiliaryPoints(fpVectorOfAuxiliaryPointsPointer);

  // copy of G4ParticleChange::UpdateStepForAlongStep
  //  i.e. no effect for touchable

  // A physics process always calculates the final state of the
  // particle relative to the initial state at the beginning
  // of the Step, i.e., based on information of G4Track (or
  // equivalently the PreStepPoint).
  // So, the differences (delta) between these two states have to be
  // calculated and be accumulated in PostStepPoint.

  // Take note that the return type of GetMomentumChange is a
  // pointer to G4ThreeVector. Also it is a normalized
  // momentum vector.

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  G4Track*     aTrack  = pStep->GetTrack();
  G4double     mass = aTrack->GetDynamicParticle()->GetMass();

  // uodate kinetic energy
  //  now assume that no energy change in transportation
  //  However it is not true in electric fields
  //  Case for changing energy will be implemented in future


  // update momentum direction and energy
  if (isMomentumChanged) {
    G4double energy;
    energy= pPostStepPoint->GetKineticEnergy()
                 + (theEnergyChange - pPreStepPoint->GetKineticEnergy());

    // calculate new momentum
    G4ThreeVector pMomentum =  pPostStepPoint->GetMomentum()
                     + ( CalcMomentum(theEnergyChange, theMomentumDirectionChange, mass)
	                  - pPreStepPoint->GetMomentum());
    G4double      tMomentum = pMomentum.mag();
    G4ThreeVector direction(1.0,0.0,0.0);
    if( tMomentum > 0. ){
      G4double  inv_Momentum= 1.0 / tMomentum;
      direction= pMomentum * inv_Momentum;
    }
    pPostStepPoint->SetMomentumDirection(direction);
    pPostStepPoint->SetKineticEnergy( energy );
  }
  if (isVelocityChanged){  pPostStepPoint->SetVelocity(theVelocityChange);
	}

  // stop case should not occur
  //pPostStepPoint->SetMomentumDirection(G4ThreeVector(1., 0., 0.));


  // update polarization
  pPostStepPoint->AddPolarization( thePolarizationChange
  				   - pPreStepPoint->GetPolarization());

  // update position and time
  pPostStepPoint->AddPosition( thePositionChange
			       - pPreStepPoint->GetPosition() );

	///////////////////////////////////////////////////
	
	//add by katayama

	///////////////////////////////////////////////////
	// time variable makes update with logn double precision
	///////////////////////////////////////////////////
  info = (UCNTimeInformation*)(aTrack->GetUserInformation());  
	precise = (info!=NULL)?info->GetPrecise():false;
	if(info!=NULL&&precise){
	pPostStepPoint->SetGlobalTime(info->GetGlobalTime());
	pPostStepPoint->SetLocalTime(info->GetGlobalTime());
	pPostStepPoint->SetProperTime(info->GetProperTime());
	}else{
  pPostStepPoint->AddGlobalTime( theTimeChange - pPreStepPoint->GetLocalTime());
  pPostStepPoint->AddLocalTime( theTimeChange - pPreStepPoint->GetLocalTime());
  pPostStepPoint->AddProperTime( theProperTimeChange - pPreStepPoint->GetProperTime());
	}

	///////////////////////////////////////////////////

#ifdef G4VERBOSE
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes
  //pStep->SetStepLength( theTrueStepLength );
  //  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  pStep->SetControlFlag( theSteppingControlFlag );
  return pStep;
  //  return UpdateStepInfo(pStep);
}

G4Step* UCNParticleChangeForUCNTime::UpdateStepForPostStep(G4Step* pStep)
{
  // A physics process always calculates the final state of the particle

  // Change volume only if some kinetic energy remains
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  if(pPostStepPoint->GetKineticEnergy() > 0.0) {

    // update next touchable
    // (touchable can be changed only at PostStepDoIt)
    pPostStepPoint->SetTouchableHandle( theTouchableHandle );

    pPostStepPoint->SetMaterial( theMaterialChange );
    pPostStepPoint->SetMaterialCutsCouple( theMaterialCutsCoupleChange );
    pPostStepPoint->SetSensitiveDetector( theSensitiveDetectorChange );
  }
  if( this->GetLastStepInVolume() ){
    pStep->SetLastStepFlag();
  }else{
    pStep->ClearLastStepFlag(); 
  }
  // It used to call base class's method
  //   - but this would copy uninitialised data members
  // return G4ParticleChange::UpdateStepForPostStep(pStep);

  // Copying what the base class does would instead
  //   - also not useful
  // return G4VParticleChange::UpdateStepInfo(pStep);

  return pStep;
}


//----------------------------------------------------------------
// methods for printing messages
//

void UCNParticleChangeForUCNTime::DumpInfo() const
{
// use base-class DumpInfo
  G4ParticleChange::DumpInfo();

  G4int oldprc = G4cout.precision(3);
  G4cout << "        Touchable (pointer) : " 
         << std::setw(20) << theTouchableHandle() << G4endl; 
  G4cout.precision(oldprc);
}
//----------------------------------------------------------------


G4bool UCNParticleChangeForUCNTime::CheckIt(const G4Track& aTrack)
{
  G4bool    exitWithError = false;
  G4double  accuracy;
  static G4int nError = 0;
#ifdef G4VERBOSE
  const  G4int maxError = 30;
#endif

  // No check in case of "fStopAndKill" 
  if (GetTrackStatus() ==   fStopAndKill )  return G4VParticleChange::CheckIt(aTrack);

  // MomentumDirection should be unit vector
  G4bool itsOKforMomentum = true;  
  if ( theEnergyChange >0.) {
    accuracy = std::fabs(theMomentumDirectionChange.mag2()-1.0);
    if (accuracy > accuracyForWarning) {
      itsOKforMomentum = false;
      nError += 1;
      exitWithError = exitWithError || (accuracy > accuracyForException);
#ifdef G4VERBOSE
      if (nError < maxError) {
	G4cout << "  G4ParticleChange::CheckIt  : ";
	G4cout << "the Momentum Change is not unit vector !!" 
	       << "  Difference:  " << accuracy << G4endl;
	G4cout << aTrack.GetDefinition()->GetParticleName()
	       << " E=" << aTrack.GetKineticEnergy()/MeV
	       << " pos=" << aTrack.GetPosition().x()/m
	       << ", " << aTrack.GetPosition().y()/m
	       << ", " << aTrack.GetPosition().z()/m
	       <<G4endl;
      }
#endif
    }
  }

  // Both global and proper time should not go back
  G4bool itsOKforGlobalTime = true;  
	UCNTimeInformation* info = (UCNTimeInformation*)aTrack.GetUserInformation();
	precise = (info!=NULL)?info->GetPrecise():false;
	if(info!=NULL&&precise){
  accuracy = (info->GetGlobalTime0()-info->GetGlobalTime())/ns;
	}else{
  accuracy = (aTrack.GetLocalTime()- theTimeChange)/ns;
	}
	if(printverbose){
	printf("UCNParticleChangeForUCNTime::CheckIt()\n");
	printf("aTrack.GetLocalTime()=%f\n",aTrack.GetLocalTime());
	printf("theTimeChange=%f\n",theTimeChange);
	printf("accuracy=%f\n",accuracy);
	printf("accuracyForWarning=%f\n",accuracyForWarning);
	printf("accuracyForException=%f\n",accuracyForException);
	}

  if (accuracy > accuracyForWarning) {
    itsOKforGlobalTime = false;
    nError += 1;
    exitWithError = exitWithError || (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4ParticleChange::CheckIt    : ";
      G4cout << "the local time goes back  !!" 
	     << "  Difference:  " << accuracy  << "[ns] " <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     << " global time=" << aTrack.GetGlobalTime()/ns
	     << " local time=" << aTrack.GetLocalTime()/ns
	     << " proper time=" << aTrack.GetProperTime()/ns
	     << G4endl;
    }
#endif
  }

  G4bool itsOKforProperTime = true;
	UCNTimeInformation* Info = (UCNTimeInformation*)aTrack.GetUserInformation();
	precise = (Info!=NULL)?Info->GetPrecise():false;
	if(Info!=NULL&&precise){
  accuracy = (Info->GetProperTime0()-Info->GetProperTime())/ns;
	}else{
  accuracy = (aTrack.GetProperTime() - theProperTimeChange )/ns;
	}
	if(printverbose){
	printf("UCNParticleChangeForUCNTime::CheckIt()\n");
	printf("aTrack.GetProperTime()=%f\n",aTrack.GetProperTime());
	printf("theTimeChange=%f\n",theProperTimeChange);
	printf("accuracy=%f\n",accuracy);
	printf("accuracyForWarning=%f\n",accuracyForWarning);
	printf("accuracyForException=%f\n",accuracyForException);
	}

  if (accuracy > accuracyForWarning) {
    itsOKforProperTime = false;
    nError += 1;
    exitWithError = exitWithError ||  (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4ParticleChange::CheckIt    : ";
      G4cout << "the proper time goes back  !!" 
	     << "  Difference:  " << accuracy  << "[ns] " <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     << " global time=" << aTrack.GetGlobalTime()/ns
	     << " local time=" << aTrack.GetLocalTime()/ns
	     << " proper time=" << aTrack.GetProperTime()/ns
	     <<G4endl;
    }
#endif
  }

  // Kinetic Energy should not be negative
  G4bool itsOKforEnergy = true;
  accuracy = -1.0*theEnergyChange/MeV;
  if (accuracy > accuracyForWarning) {
    itsOKforEnergy = false;
    nError += 1;
    exitWithError = exitWithError ||   (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4ParticleChange::CheckIt    : ";
      G4cout << "the kinetic energy is negative  !!" 
	     << "  Difference:  " << accuracy  << "[MeV] " <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     <<G4endl;
    }
#endif
  }

  // Velocity  should not be less than c_light
  G4bool itsOKforVelocity = true;
  if (theVelocityChange < 0.) {
    itsOKforVelocity = false;
    nError += 1;
    exitWithError = true;
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4ParticleChange::CheckIt    : ";
      G4cout << "the velocity is negative  !!" 
	     << "  Velocity:  " << theVelocityChange/c_light  <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     <<G4endl;
    }
#endif
  }

  accuracy = theVelocityChange/c_light - 1.0;
  if (accuracy > accuracyForWarning) {
    itsOKforVelocity = false;
    nError += 1;
    exitWithError = exitWithError ||  (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4ParticleChange::CheckIt    : ";
      G4cout << "the velocity is greater than c_light  !!" << G4endl;
      G4cout << "  Velocity:  " << theVelocityChange/c_light  <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     <<G4endl;
    }
#endif
  }

  G4bool itsOK = itsOKforMomentum && itsOKforEnergy && itsOKforVelocity && itsOKforProperTime && itsOKforGlobalTime;
  // dump out information of this particle change
#ifdef G4VERBOSE
  if (!itsOK) { 
    DumpInfo();
  }
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4ParticleChange::CheckIt",
		"TRACK003", EventMustBeAborted,
		"momentum, energy, and/or time was illegal");
  }
  //correction
  if (!itsOKforMomentum) {
    G4double vmag = theMomentumDirectionChange.mag();
    theMomentumDirectionChange = (1./vmag)*theMomentumDirectionChange;
  }
  if (!itsOKforGlobalTime) {
    theTimeChange = aTrack.GetLocalTime();
  }
  if (!itsOKforProperTime) {
    theProperTimeChange = aTrack.GetProperTime();
  }
  if (!itsOKforEnergy) {
    theEnergyChange = 0.0;
  }
  if (!itsOKforVelocity) {
    theVelocityChange = c_light;
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}

UCNTimeInformation* UCNParticleChangeForUCNTime::GetUserInfo(){
	return info;
}
