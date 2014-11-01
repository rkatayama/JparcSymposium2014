/////////////////////////////////////////////////////////////////////////////
//
//   Implemented for the scheme, based on G4ParticleChangeForTransport 10,2014 Ryo Katayama
//
/////////////////////////////////////////////////////////////////////////////

// The class can make time variable updated with long double precision just after AlongStepDoIt.

#ifndef UCNParticleChangeForUCNTime_h
#define UCNParticleChangeForUCNTime_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4TouchableHandle.hh"
#include "G4ParticleChange.hh"

#include "UCNTimeInformation.hh"

class G4MaterialCutsCouple;
class G4VSensitiveDetector;

class UCNParticleChangeForUCNTime: public G4ParticleChange
{ 
  public:
    // default constructor
    UCNParticleChangeForUCNTime();

    // destructor
    virtual ~UCNParticleChangeForUCNTime();

  protected:
    // hide copy constructor and assignment operator as protected
    UCNParticleChangeForUCNTime(const UCNParticleChangeForUCNTime &right);
    UCNParticleChangeForUCNTime & operator=(const UCNParticleChangeForUCNTime &right);

  public: // with description
    // ----------------------------------------------------
    // --- the following methods are for updating G4Step -----   
    // Return the pointer to the G4Step after updating the Step information
    // by using final state information of the track given by a physics
    // process    
    virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
    virtual G4Step* UpdateStepForAtRest(G4Step* Step);
    virtual G4Step* UpdateStepForPostStep(G4Step* Step);
    // A physics process gives the final state of the particle 
    // based on information of G4Track (or equivalently the PreStepPoint)
 
    virtual void Initialize(const G4Track&);
    // Initialize all propoerties by using G4Track information
           
    // ----------------------------------------------------
    //--- methods to keep information of the final state--
    //  IMPORTANT NOTE: Although the name of the class and methods are
    //   "Change", what it stores (and returns in get) are the "FINAL" 
    //   values of the Position, Momentum, etc.

    const G4TouchableHandle& GetTouchableHandle() const;
    void  SetTouchableHandle(const G4TouchableHandle& fTouchable);
    //  Get/Set the touchable of the current particle.
    //  Note: Touchable in PostStepPoint will be updated only after PostStepDoIt

    G4Material* GetMaterialInTouchable() const;
    void SetMaterialInTouchable(G4Material* fMaterial);
    //  Get/Propose the material in the touchable of the current particle.

    const G4MaterialCutsCouple* GetMaterialCutsCoupleInTouchable() const;
    void SetMaterialCutsCoupleInTouchable(const G4MaterialCutsCouple* fMaterialCutsCouple);
    //  Get/Set the materialCutsCouple in the touchable of the current particle.

    G4VSensitiveDetector* GetSensitiveDetectorInTouchable() const;
    void SetSensitiveDetectorInTouchable(G4VSensitiveDetector* fSensitiveDetector);
    //  Get/Set the sensitive detector in the touchable of the current particle.

    G4bool GetMomentumChanged() const;
    void SetMomentumChanged(G4bool b);

  public:
    virtual void DumpInfo() const;

  protected:
    G4TouchableHandle theTouchableHandle;
    //  The changed touchable of a given particle.

  public:

    // Prototype implementation of smooth representation of curved trajectories.
    // Auxiliary points are ThreeVectors for now; change to G4AuxiliaryPoints.

    inline void SetPointerToVectorOfAuxiliaryPoints( std::vector<G4ThreeVector>* theNewVectorPointer );
    inline std::vector<G4ThreeVector>* GetPointerToVectorOfAuxiliaryPoints() const;

  private:
    G4bool     isMomentumChanged;
    //  The flag which is set if momentum is changed in current step
    G4Material* theMaterialChange;
    const G4MaterialCutsCouple* theMaterialCutsCoupleChange;
    G4VSensitiveDetector* theSensitiveDetectorChange;
     // The material (and MaterialCutsCouple) where given track
     // currently locates

  private:
    std::vector<G4ThreeVector>* fpVectorOfAuxiliaryPointsPointer;

		///////////////////////////////////////////
		//add by katayama
	public:
		UCNTimeInformation* info;
		G4bool precise;
		G4bool CheckIt(const G4Track& aTrack);
		G4bool printverbose;
	public:
		UCNTimeInformation* GetUserInfo();
		///////////////////////////////////////////
};

#include "UCNParticleChangeForUCNTime.icc"

#endif
