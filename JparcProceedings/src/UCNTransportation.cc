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
// $Id: UCNTransportation.cc,v 1.5 2008/09/17 09:07:58 kuzniak Exp $
// GEANT4 tag $Name:  $
// 
// ------------------------------------------------------------
//  GEANT 4  include file implementation
//
// ------------------------------------------------------------
//
// This class is a process responsible for the transportation of 
// a particle, ie the geometrical propagation that encounters the 
// geometrical sub-volumes of the detectors.
//
// It is also tasked with part of updating the "safety".
//
// =======================================================================
// Modified:   
//            19 Jan  2006, P.MoraDeFreitas: Fix for suspended tracks (StartTracking)
//            11 Aug  2004, M.Asai: Add G4VSensitiveDetector* for updating stepPoint.
//            21 June 2003, J.Apostolakis: Calling field manager with 
//                            track, to enable it to configure its accuracy
//            13 May  2003, J.Apostolakis: Zero field areas now taken into
//                            account correclty in all cases (thanks to W Pokorski).
//            29 June 2001, J.Apostolakis, D.Cote-Ahern, P.Gumplinger: 
//                          correction for spin tracking   
//            20 Febr 2001, J.Apostolakis:  update for new FieldTrack
//            22 Sept 2000, V.Grichine:     update of Kinetic Energy
// Created:  19 March 1997, J. Apostolakis
// =======================================================================

#include "UCNTransportation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleTable.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManagerStore.hh"
#include "UCNTimeInformation.hh"
#include "UCNTransportationMessenger.hh"
class G4VSensitiveDetector;

//////////////////////////////////////////////////////////////////////////
//
// Constructor

UCNTransportation::UCNTransportation( G4int verboseLevel )
  : G4VProcess( G4String("Transportation"), fTransportation ),
    fParticleIsLooping( false ),
    fPreviousSftOrigin (0.,0.,0.),
    fPreviousSafety    ( 0.0 ),
    fThreshold_Warning_Energy( 100 * MeV ),  
    fThreshold_Important_Energy( 250 * MeV ), 
    fThresholdTrials( 10 ), 
    fUnimportant_Energy( 1 * MeV ), 
    fNoLooperTrials(0),
    fSumEnergyKilled( 0.0 ), fMaxEnergyKilled( 0.0 ), 
    fShortStepOptimisation(false),    // Old default: true (=fast short steps)
    fVerboseLevel( verboseLevel ),
		printverbose(false),
		precise(true)
{
  G4TransportationManager* transportMgr ; 

  transportMgr = G4TransportationManager::GetTransportationManager() ; 

  fLinearNavigator = transportMgr->GetNavigatorForTracking() ; 

  // fGlobalFieldMgr = transportMgr->GetFieldManager() ;

  fFieldPropagator = transportMgr->GetPropagatorInField() ;

  // fpSafetyHelper = fpTransportManager->GetSafetyHelper();  // New 

  // Cannot determine whether a field exists here,
  //  because it would only work if the field manager has informed 
  //  about the detector's field before this transportation process 
  //  is constructed.
  // Instead later the method DoesGlobalFieldExist() is called

  static G4TouchableHandle nullTouchableHandle;  // Points to (G4VTouchable*) 0
  fCurrentTouchableHandle = nullTouchableHandle; 

  fEndGlobalTimeComputed  = false;
  fCandidateEndGlobalTime = 0; 

	theMessenger = new UCNTransportationMessenger(this);
}

//////////////////////////////////////////////////////////////////////////

UCNTransportation::~UCNTransportation()
{

  // --- Alternative code to delete 'junk data' in touchable handle
  //  ** Corresponds to the case that the constructor has
  //       fCurrentTouchableHandle = new G4TouchableHistory();
  //  ** Likely incorrect - as at the end of the simulation 
  //     the handle no longer holds the original 'null' history
  // G4VTouchable*  pTouchable= fCurrentTouchableHandle();  //  Incorrect
  // delete  pTouchable;  // Incorrect

  if( (fVerboseLevel > 0) && (fSumEnergyKilled > 0.0 ) ){ 
    G4cout << " UCNTransportation: Statistics for looping particles " << G4endl;
    G4cout << "   Sum of energy of loopers killed: " <<  fSumEnergyKilled << G4endl;
    G4cout << "   Max energy of loopers killed: " <<  fMaxEnergyKilled << G4endl;
  } 
}

//////////////////////////////////////////////////////////////////////////
//
// Responsibilities:
//    Find whether the geometry limits the Step, and to what length
//    Calculate the new value of the safety and return it.
//    Store the final time, position and momentum.

G4double UCNTransportation::
AlongStepGetPhysicalInteractionLength( const G4Track&  track,
                                             G4double, //  previousStepSize
                                             G4double  currentMinimumStep,
                                             G4double& currentSafety,
                                             G4GPILSelection* selection )
{
  G4double geometryStepLength, newSafety ; 
  fParticleIsLooping = false ;

  // Initial actions moved to  StartTrack()   
  // --------------------------------------
  // Note: in case another process changes touchable handle
  //    it will be necessary to add here (for all steps)   
  // fCurrentTouchableHandle = aTrack->GetTouchableHandle();

  // GPILSelection is set to defaule value of CandidateForSelection
  // It is a return value
  //
  *selection = CandidateForSelection ;

  // Get initial Energy/Momentum of the track
  //
  const G4DynamicParticle*    pParticle  = track.GetDynamicParticle() ;
  const G4ParticleDefinition* pParticleDef   = pParticle->GetDefinition() ;
  G4ThreeVector startMomentumDir       = pParticle->GetMomentumDirection() ;
  G4ThreeVector startPosition          = track.GetPosition() ;

  // G4double   theTime        = track.GetGlobalTime() ;
  // The Step Point safety can be limited by other geometries and/or the 
  // assumptions of any process - it's not always the geometrical safety.
  // We calculate the starting point's isotropic safety here.
  //
  G4ThreeVector OriginShift = startPosition - fPreviousSftOrigin ;
  G4double      MagSqShift  = OriginShift.mag2() ;
  if( MagSqShift >= sqr(fPreviousSafety) )
  {
     currentSafety = 0.0 ;
  }
  else
  {
     currentSafety = fPreviousSafety - std::sqrt(MagSqShift) ;
  }

  // Is the particle charged ?
  //
  G4double              particleCharge = pParticle->GetCharge() ; 

  fGeometryLimitedStep = false ;
  // fEndGlobalTimeComputed = false ;

  // There is no need to locate the current volume. It is Done elsewhere:
  //   On track construction 
  //   By the tracking, after all AlongStepDoIts, in "Relocation"

  // Check whether the particle have an (EM) field force exerting upon it
  //
  G4FieldManager* fieldMgr=0;
  G4bool          fieldExertsForce = false ;

  if( pParticleDef->GetParticleName() == "neutron" ||
      pParticleDef->GetParticleName() == "ucn" || 
      particleCharge != 0.0 )
  {
     fieldMgr= fFieldPropagator->FindAndSetFieldManager( track.GetVolume() ); 
     if (fieldMgr != 0) {
	// Message the field Manager, to configure it for this track
	fieldMgr->ConfigureForTrack( &track );
	// Moved here, in order to allow a transition
	//   from a zero-field  status (with fieldMgr->(field)0
	//   to a finite field  status

        // If the field manager has no field, there is no field !
        fieldExertsForce = (fieldMgr->GetDetectorField() != 0);
     } 
  }

  // G4cout << " G4Transport:  field exerts force= " << fieldExertsForce
  // 	 << "  fieldMgr= " << fieldMgr << G4endl;

  // Choose the calculation of the transportation: Field or not 
  //
  if( !fieldExertsForce ) 
  {
     G4double linearStepLength ;
     if( fShortStepOptimisation && (currentMinimumStep <= currentSafety) )
     {
       // The Step is guaranteed to be taken
       //
       geometryStepLength   = currentMinimumStep ;
       fGeometryLimitedStep = false ;
     }
     else
     {
       //  Find whether the straight path intersects a volume
       //
       linearStepLength = fLinearNavigator->ComputeStep( startPosition, 
                                                         startMomentumDir,
                                                         currentMinimumStep, 
                                                         newSafety) ;
       // Remember last safety origin & value.
       //
       fPreviousSftOrigin = startPosition ;
       fPreviousSafety    = newSafety ; 
       // fSafetyHelper->SetCurrentSafety( newSafety, startPosition);

       // The safety at the initial point has been re-calculated:
       //
       currentSafety = newSafety ;
          
       fGeometryLimitedStep= (linearStepLength <= currentMinimumStep); 
       if( fGeometryLimitedStep )
       {
         // The geometry limits the Step size (an intersection was found.)
         geometryStepLength   = linearStepLength ;
       } 
       else
       {
         // The full Step is taken.
         geometryStepLength   = currentMinimumStep ;
       }
     }
     endpointDistance = geometryStepLength ;

     // Calculate final position
     //
     fTransportEndPosition = startPosition+geometryStepLength*startMomentumDir ;

     // Momentum direction, energy and polarisation are unchanged by transport
     //
     fTransportEndMomentumDir   = startMomentumDir ; 
     fTransportEndKineticEnergy = track.GetKineticEnergy() ;
     fTransportEndSpin          = track.GetPolarization();
     fParticleIsLooping         = false ;
     fMomentumChanged           = false ; 
     fEndGlobalTimeComputed     = false ;
  }
  else   //  A field exerts force
  {
     G4double       momentumMagnitude = pParticle->GetTotalMomentum() ;
     G4ThreeVector  EndUnitMomentum ;
     G4double       lengthAlongCurve ;
     G4double       restMass = pParticleDef->GetPDGMass() ;
 
     fFieldPropagator->SetChargeMomentumMass( particleCharge,    // in e+ units
                                              momentumMagnitude, // in Mev/c 
                                              restMass           ) ;  
     
    
		 if(printverbose){
		 G4cout << "--------AlongStepGPIL-------" << G4endl;
		 G4cout << "check track time: " << track.GetGlobalTime() << G4endl;
		 G4cout << "----------------------------" << G4endl;
		 }

     G4ThreeVector spin        = track.GetPolarization() ;

		 //////////////////////////////////////////////////////////////////////////////
		 //add by kurata and katayama
		 G4double GT = ( (pParticleDef->GetParticleName()=="ucn") && precise )?0:(track.GetGlobalTime());
		 G4double PT = ( (pParticleDef->GetParticleName()=="ucn") && precise )?0:(track.GetProperTime());

     G4FieldTrack  aFieldTrack = G4FieldTrack( startPosition, 
                                               track.GetMomentumDirection(),
                                               0.0, 
                                               track.GetKineticEnergy(),
                                               restMass,
                                               track.GetVelocity(),
                                               GT, //track.GetGlobalTime(), // Lab.
                                               PT, //track.GetProperTime(), // Part.
                                               &spin                  ) ;

		 //////////////////////////////////////////////////////////////////////////////

     if( currentMinimumStep > 0 ) 
     {
        // Do the Transport in the field (non recti-linear)
        //
        lengthAlongCurve = fFieldPropagator->ComputeStep( aFieldTrack,
                                                          currentMinimumStep, 
                                                          currentSafety,
                                                          track.GetVolume() ) ;
        fGeometryLimitedStep= lengthAlongCurve < currentMinimumStep; 
	if( fGeometryLimitedStep ) {
           geometryStepLength   = lengthAlongCurve ;
        } else {
           geometryStepLength   = currentMinimumStep ;
        }
     }
     else
     {
        geometryStepLength   = lengthAlongCurve= 0.0 ;
        fGeometryLimitedStep = false ;
     }

     // Remember last safety origin & value.
     //
     fPreviousSftOrigin = startPosition ;
     fPreviousSafety    = currentSafety ;         
     // fSafetyHelper->SetCurrentSafety( newSafety, startPosition);
       
     // Get the End-Position and End-Momentum (Dir-ection)
     //
     fTransportEndPosition = aFieldTrack.GetPosition() ;

     // Momentum:  Magnitude and direction can be changed too now ...
     //
     fMomentumChanged         = true ; 
     fTransportEndMomentumDir = aFieldTrack.GetMomentumDir() ;

     fTransportEndKineticEnergy  = aFieldTrack.GetKineticEnergy() ; 

     if( fFieldPropagator->GetCurrentFieldManager()->DoesFieldChangeEnergy() )
     {
        // If the field can change energy, then the time must be integrated
        //    - so this should have been updated
        //
        fCandidateEndGlobalTime   = aFieldTrack.GetLabTimeOfFlight();
        fEndGlobalTimeComputed    = true;

        // was ( fCandidateEndGlobalTime != track.GetGlobalTime() );
        // a cleaner way is to have FieldTrack knowing whether time is updated.
     }
     else
     {
        // The energy should be unchanged by field transport,
        //    - so the time changed will be calculated elsewhere
        //
        fEndGlobalTimeComputed = false;

        // Check that the integration preserved the energy 
        //     -  and if not correct this!
        G4double  startEnergy= track.GetKineticEnergy();
        G4double  endEnergy= fTransportEndKineticEnergy; 

        static G4int no_inexact_steps=0, no_large_ediff;
        G4double absEdiff = std::fabs(startEnergy- endEnergy);
        if( absEdiff > perMillion * endEnergy )
        {
          no_inexact_steps++;
          // Possible statistics keeping here ...
        }
        if( fVerboseLevel > 1 )
        {
          if( std::fabs(startEnergy- endEnergy) > perThousand * endEnergy )
          {
            static G4int no_warnings= 0, warnModulo=1,  moduloFactor= 10; 
            no_large_ediff ++;
            if( (no_large_ediff% warnModulo) == 0 )
            {
               no_warnings++;
               G4cout << "WARNING - UCNTransportation::AlongStepGetPIL() " 
	              << "   Energy change in Step is above 1^-3 relative value. " << G4endl
		      << "   Relative change in 'tracking' step = " 
		      << std::setw(15) << (endEnergy-startEnergy)/startEnergy << G4endl
                      << "     Starting E= " << std::setw(12) << startEnergy / MeV << " MeV " << G4endl
                      << "     Ending   E= " << std::setw(12) << endEnergy   / MeV << " MeV " << G4endl;       
               G4cout << " Energy has been corrected -- however, review"
                      << " field propagation parameters for accuracy."  << G4endl;
	       if( (fVerboseLevel > 2 ) || (no_warnings<4) || (no_large_ediff == warnModulo * moduloFactor) ){
		 G4cout << " These include EpsilonStepMax(/Min) in G4FieldManager "
			<< " which determine fractional error per step for integrated quantities. " << G4endl
			<< " Note also the influence of the permitted number of integration steps."
			<< G4endl;
	       }
               G4cerr << "ERROR - UCNTransportation::AlongStepGetPIL()" << G4endl
	              << "        Bad 'endpoint'. Energy change detected"
                      << " and corrected. " 
		      << " Has occurred already "
                      << no_large_ediff << " times." << G4endl;
               if( no_large_ediff == warnModulo * moduloFactor )
               {
                  warnModulo *= moduloFactor;
               }
            }
          }
        }  // end of if (fVerboseLevel)

        // Correct the energy for fields that conserve it
        //  This - hides the integration error
        //       - but gives a better physical answer
        fTransportEndKineticEnergy= track.GetKineticEnergy(); 
     }

     fTransportEndSpin = aFieldTrack.GetSpin();
     fParticleIsLooping = fFieldPropagator->IsParticleLooping() ;
     endpointDistance   = (fTransportEndPosition - startPosition).mag() ;
  }

  // If we are asked to go a step length of 0, and we are on a boundary
  // then a boundary will also limit the step -> we must flag this.
  //
  if( currentMinimumStep == 0.0 ) 
  {
      if( currentSafety == 0.0 )  fGeometryLimitedStep = true ;
  }

  // Update the safety starting from the end-point,
  // if it will become negative at the end-point.
  //
  if( currentSafety < endpointDistance ) 
  {
      G4double endSafety =
               fLinearNavigator->ComputeSafety( fTransportEndPosition) ;
      currentSafety      = endSafety ;
      fPreviousSftOrigin = fTransportEndPosition ;
      fPreviousSafety    = currentSafety ; 
      // fSafetyHelper->SetCurrentSafety( currentSafety, fTransportEndPosition);

      // Because the Stepping Manager assumes it is from the start point, 
      //  add the StepLength
      //
      currentSafety     += endpointDistance ;

#ifdef G4DEBUG_TRANSPORT 
      G4cout.precision(12) ;
      G4cout << "***UCNTransportation::AlongStepGPIL ** " << G4endl  ;
      G4cout << "  Called Navigator->ComputeSafety at " << fTransportEndPosition
             << "    and it returned safety= " << endSafety << G4endl ; 
      G4cout << "  Adding endpoint distance " << endpointDistance 
             << "   to obtain pseudo-safety= " << currentSafety << G4endl ; 
#endif
  }            

  fParticleChange.ProposeTrueStepLength(geometryStepLength) ;

  return geometryStepLength ;
}

//////////////////////////////////////////////////////////////////////////
//
//   Initialize ParticleChange  (by setting all its members equal
//                               to corresponding members in G4Track)

G4VParticleChange* UCNTransportation::AlongStepDoIt( const G4Track& track,
                                                    const G4Step&  stepData )
{

  const G4DynamicParticle*    pParticle  = track.GetDynamicParticle() ;
  const G4ParticleDefinition* pParticleDef   = pParticle->GetDefinition() ;
	if(printverbose){
	printf("-------AlongStepDoIt--------\n");
	printf("track Proper Time=%f\n",pParticle->GetProperTime());
	printf("track Local Time=%f\n",track.GetLocalTime());
	printf("----------------------------\n");
	}

  static G4int noCalls=0;
  static const G4ParticleDefinition* fOpticalPhoton =
           G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

  noCalls++;


  
  fParticleChange.Initialize(track) ;

  //  Code for specific process 
  //
  fParticleChange.ProposePosition(fTransportEndPosition) ;
  fParticleChange.ProposeMomentumDirection(fTransportEndMomentumDir) ;
  fParticleChange.ProposeEnergy(fTransportEndKineticEnergy) ;
  fParticleChange.SetMomentumChanged(fMomentumChanged) ;

  fParticleChange.ProposePolarization(fTransportEndSpin);
  
  G4double deltaTime = 0.0 ;

  // Calculate  Lab Time of Flight (ONLY if field Equations used it!)
     // G4double endTime   = fCandidateEndGlobalTime;
     // G4double delta_time = endTime - startTime;

  G4double startTime = 0.0;    //track.GetGlobalTime() ;  //here needs to change to avoid overflow!

  if (!fEndGlobalTimeComputed)
  {
     // The time was not integrated .. make the best estimate possible
     //
     G4double finalVelocity   = track.GetVelocity() ;
     G4double initialVelocity = stepData.GetPreStepPoint()->GetVelocity() ;
     G4double stepLength      = track.GetStepLength() ;

     deltaTime= 0.0;  // in case initialVelocity = 0 
     const G4DynamicParticle* fpDynamicParticle = track.GetDynamicParticle();
     if (fpDynamicParticle->GetDefinition()== fOpticalPhoton)
     {
        //  A photon is in the medium of the final point
        //  during the step, so it has the final velocity.
        deltaTime = stepLength/finalVelocity ;
     }
     else if (finalVelocity > 0.0)
     {
        G4double meanInverseVelocity ;
        // deltaTime = stepLength/finalVelocity ;
        meanInverseVelocity = 0.5
                            * ( 1.0 / initialVelocity + 1.0 / finalVelocity ) ;
        deltaTime = stepLength * meanInverseVelocity ;
     }
     else if( initialVelocity > 0.0 )
     {
        deltaTime = stepLength/initialVelocity ;
     }
     fCandidateEndGlobalTime   = startTime + deltaTime ;
  }
  else
  {
     deltaTime = fCandidateEndGlobalTime - startTime ;
  }

  //G4cout << "times - stepdoit: " << startTime << " " << deltaTime << " " << fCandidateEndGlobalTime << G4endl;

  fParticleChange.ProposeGlobalTime( fCandidateEndGlobalTime ) ;

  // Now Correct by Lorentz factor to get "proper" deltaTime
  
  G4double  restMass       = track.GetDynamicParticle()->GetMass() ;
  G4double deltaProperTime = deltaTime*( restMass/track.GetTotalEnergy() ) ;

  fParticleChange.ProposeProperTime(deltaProperTime + track.GetProperTime()) ;
  //fParticleChange. ProposeTrueStepLength( track.GetStepLength() ) ;

  // If the particle is caught looping or is stuck (in very difficult
  // boundaries) in a magnetic field (doing many steps) 
  //   THEN this kills it ...
  //
  if ( fParticleIsLooping )
  {
      G4double endEnergy= fTransportEndKineticEnergy;

      if( (endEnergy < fThreshold_Important_Energy) 
	  || (fNoLooperTrials >= fThresholdTrials ) ){
	// Kill the looping particle 
	//
	fParticleChange.ProposeTrackStatus( fStopAndKill )  ;

        // 'Bare' statistics
        fSumEnergyKilled += endEnergy; 
	if( endEnergy > fMaxEnergyKilled) { fMaxEnergyKilled= endEnergy; }

#ifdef G4VERBOSE
	if( (fVerboseLevel > 1) || 
	    ( endEnergy > fThreshold_Warning_Energy )  ) { 
	  G4cout << " UCNTransportation is killing track that is looping or stuck "
		 << G4endl
		 << "   This track has " << track.GetKineticEnergy() / MeV
		 << " MeV energy." << G4endl;
	  G4cout << "   Number of trials = " << fNoLooperTrials 
		 << "   No of calls to AlongStepDoIt = " << noCalls 
		 << G4endl;
	}
#endif
	fNoLooperTrials=0; 
      }
      else{
	fNoLooperTrials ++; 
#ifdef G4VERBOSE
	if( (fVerboseLevel > 2) ){
	  G4cout << "   UCNTransportation::AlongStepDoIt(): Particle looping -  "
		 << "   Number of trials = " << fNoLooperTrials 
		 << "   No of calls to  = " << noCalls 
		 << G4endl;
	}
#endif
      }
  }else{
      fNoLooperTrials=0; 
  }

  // Another (sometimes better way) is to use a user-limit maximum Step size
  // to alleviate this problem .. 

  // Introduce smooth curved trajectories to particle-change
  //
  fParticleChange.SetPointerToVectorOfAuxiliaryPoints
    (fFieldPropagator->GimmeTrajectoryVectorAndForgetIt() );


	//////////////////////////////////////////////////////////////
	//add by kurata and katayama
	if( (pParticleDef->GetParticleName() == "ucn") && precise ){

      G4Track* aTrack = stepData.GetTrack();
			UCNTimeInformation *anInfo = fParticleChange.GetUserInfo();

  //At the first step, associate the UCNTimeInformation with G4VUserInformation that G4track has
  if(track.GetParentID()==0 && track.GetUserInformation()==0)
    {
      anInfo = new UCNTimeInformation(aTrack);
      aTrack->SetUserInformation(anInfo);
			anInfo->SetPrecise(precise);
    }//end if

  UCNTimeInformation* Info = (UCNTimeInformation*)(track.GetUserInformation());  
  Info->UpdateUCNTimeInfo(deltaTime);
  Info->UpdateUCNProperTimeInfo(deltaProperTime);

	}//end if
	//////////////////////////////////////////////////////////////

  return &fParticleChange ;
}

//////////////////////////////////////////////////////////////////////////
//
//  This ensures that the PostStep action is always called,
//  so that it can do the relocation if it is needed.
// 

G4double UCNTransportation::
PostStepGetPhysicalInteractionLength( const G4Track&,
                                            G4double, // previousStepSize
                                            G4ForceCondition* pForceCond )
{ 
  *pForceCond = Forced ; 
  return DBL_MAX ;  // was kInfinity ; but convention now is DBL_MAX
}

/////////////////////////////////////////////////////////////////////////////
//

G4VParticleChange* UCNTransportation::PostStepDoIt( const G4Track& track,
                                                   const G4Step& )
{
  G4TouchableHandle retCurrentTouchable ;   // The one to return

  // Initialize ParticleChange  (by setting all its members equal
  //                             to corresponding members in G4Track)
  // fParticleChange.Initialize(track) ;  // To initialise TouchableChange

  fParticleChange.ProposeTrackStatus(track.GetTrackStatus()) ;

  // If the Step was determined by the volume boundary,
  // logically relocate the particle
  if(fGeometryLimitedStep)
  {  
    // fCurrentTouchable will now become the previous touchable, 
    // and what was the previous will be freed.
    // (Needed because the preStepPoint can point to the previous touchable)

    fLinearNavigator->SetGeometricallyLimitedStep() ;
    fLinearNavigator->
    LocateGlobalPointAndUpdateTouchableHandle( track.GetPosition(),
                                               track.GetMomentumDirection(),
                                               fCurrentTouchableHandle,
                                               true                      ) ;
    // Check whether the particle is out of the world volume 
    // If so it has exited and must be killed.
    //
    if( fCurrentTouchableHandle->GetVolume() == 0 )
    {
       fParticleChange.ProposeTrackStatus( fStopAndKill ) ;
    }
    retCurrentTouchable = fCurrentTouchableHandle ;
    fParticleChange.SetTouchableHandle( fCurrentTouchableHandle ) ;
  }
  else                 // fGeometryLimitedStep  is false
  {                    
    // This serves only to move the Navigator's location
    //
    fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition() ) ;

    // The value of the track's current Touchable is retained. 
    // (and it must be correct because we must use it below to
    // overwrite the (unset) one in particle change)
    //  It must be fCurrentTouchable too ??
    //
    fParticleChange.SetTouchableHandle( track.GetTouchableHandle() ) ;
    retCurrentTouchable = track.GetTouchableHandle() ;
  }         // endif ( fGeometryLimitedStep ) 

  const G4VPhysicalVolume* pNewVol = retCurrentTouchable->GetVolume() ;
  const G4Material* pNewMaterial   = 0 ;
  const G4VSensitiveDetector* pNewSensitiveDetector   = 0 ;
                                                                                       
  if( pNewVol != 0 )
  {
    pNewMaterial= pNewVol->GetLogicalVolume()->GetMaterial();
    pNewSensitiveDetector= pNewVol->GetLogicalVolume()->GetSensitiveDetector();
  }

  // ( <const_cast> pNewMaterial ) ;
  // ( <const_cast> pNewSensitiveDetector) ;

  fParticleChange.SetMaterialInTouchable( (G4Material *) pNewMaterial ) ;
  fParticleChange.SetSensitiveDetectorInTouchable( (G4VSensitiveDetector *) pNewSensitiveDetector ) ;

  const G4MaterialCutsCouple* pNewMaterialCutsCouple = 0;
  if( pNewVol != 0 )
  {
    pNewMaterialCutsCouple=pNewVol->GetLogicalVolume()->GetMaterialCutsCouple();
  }

  if( pNewVol!=0 && pNewMaterialCutsCouple!=0 && pNewMaterialCutsCouple->GetMaterial()!=pNewMaterial )
  {
    // for parametrized volume
    //
    pNewMaterialCutsCouple =
      G4ProductionCutsTable::GetProductionCutsTable()
                             ->GetMaterialCutsCouple(pNewMaterial,
                               pNewMaterialCutsCouple->GetProductionCuts());
  }
  fParticleChange.SetMaterialCutsCoupleInTouchable( pNewMaterialCutsCouple );

  // temporarily until Get/Set Material of ParticleChange, 
  // and StepPoint can be made const. 
  // Set the touchable in ParticleChange
  // this must always be done because the particle change always
  // uses this value to overwrite the current touchable pointer.
  //
  fParticleChange.SetTouchableHandle(retCurrentTouchable) ;

  return &fParticleChange ;
}

// New method takes over the responsibility to reset the state of UCNTransportation
//   object at the start of a new track or the resumption of a suspended track. 

void 
UCNTransportation::StartTracking(G4Track* aTrack)
{
  G4VProcess::StartTracking(aTrack);

// The actions here are those that were taken in AlongStepGPIL
//   when track.GetCurrentStepNumber()==1

  // reset safety value and center
  //
  fPreviousSafety    = 0.0 ; 
  fPreviousSftOrigin = G4ThreeVector(0.,0.,0.) ;
  
  // reset looping counter -- for motion in field
  fNoLooperTrials= 0; 
  // Must clear this state .. else it depends on last track's value
  //  --> a better solution would set this from state of suspended track TODO ? 
  // Was if( aTrack->GetCurrentStepNumber()==1 ) { .. }

  // ChordFinder reset internal state
  //
  if( DoesGlobalFieldExist() ) {
     fFieldPropagator->ClearPropagatorState();   
     // Resets all state of field propagator class (ONLY)
     //  including safety values (in case of overlaps and to wipe for first track).

     // G4ChordFinder* chordF= fFieldPropagator->GetChordFinder();
     // if( chordF ) chordF->ResetStepEstimate();
  }

  // Make sure to clear the chord finders of all fields (ie managers)
  static G4FieldManagerStore* fieldMgrStore= G4FieldManagerStore::GetInstance();
  fieldMgrStore->ClearAllChordFindersState(); 

  // Update the current touchable handle  (from the track's)
  //
  fCurrentTouchableHandle = aTrack->GetTouchableHandle();
}

void UCNTransportation::SetPreciseTime(G4int a){

	precise = (a==0)?false:true;

}
