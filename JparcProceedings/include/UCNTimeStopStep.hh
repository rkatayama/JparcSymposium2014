//////////////////////////////////////////////////////////////////////
//  the class to make an event ended at the time user wants to stop
//  7.2014 Ryo katayama
//  9.2014 changed as it can get stoptime from its messenger by Ryo Katayama
////////////////////////////////////////////////////////////////////////


#ifndef UCNTimeStopStep_h
#define UCNTimeStopStep_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "UCNUCN.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"

#include "UCNTimeStopStepMessenger.hh"

#include <TTree.h>
#include "TMath.h"
#include <TFile.h>
#include "UCNTimeInformation.hh"

class UCNTimeStopStep : public G4VContinuousDiscreteProcess
{

	public: 

		UCNTimeStopStep(const G4String& processName = "UCNTimeStopStep");
		~UCNTimeStopStep();

	public: 

		double x,y,z,px,py,pz,t;
		double spinx,spiny,spinz,phase;
		double ddph,ddsy;
		int evt;
		int rot;
		int edir;
		double b;
		double w;
		double V;
		G4int textid;

		TTree *tree;

		void SetRootFileName(G4String newValue){
			G4cout<<newValue.c_str()<<G4endl;
			RootFileName = newValue;
		}

		void SetStopTime(G4double stoptime);

		inline G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
		G4double GetContinuousStepLimit(const G4Track& aTrack,
				G4double  previousStepSize,
				G4double  currentMinimumStep,
				G4double& currentSafety);

		G4double GetMeanFreePath(const G4Track& aTrack,
				G4double ,
				G4ForceCondition* condition);

		G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
				const G4Step&  aStep);

		G4VParticleChange* AlongStepDoIt(
				const G4Track&,
				const G4Step&);

		static UCNTimeStopStep* GetInstance();

	private:
		UCNTimeStopStepMessenger * theMessenger;
		G4String RootFileName;
		static UCNTimeStopStep* theInstance;
		G4double stoptime;
};

////////////////////
// Inline methods
////////////////////

inline
G4bool UCNTimeStopStep::IsApplicable(const G4ParticleDefinition& 
		aParticleType)
{
	return ( &aParticleType == UCNUCN::UCN0() );
}

#endif /* UCNTimeStopStep_h */
