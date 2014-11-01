//////////////////////////////////////////////////////////////////////
//  the class to store the time with long double precision for UCN
//  7.2014 Masakazu Kurata
////////////////////////////////////////////////////////////////////////
#ifndef UCNTimeInformation_hh
#define UCNTimeInformation_hh 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class UCNTimeInformation : public G4VUserTrackInformation 
{
  
 public:
  UCNTimeInformation();
  UCNTimeInformation(G4Track* aTrack);
  UCNTimeInformation(const UCNTimeInformation* aTrackInfo);
  virtual ~UCNTimeInformation();
  
  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);
  inline int operator ==(const UCNTimeInformation& right) const
    {return (this==&right);}
  
  void UpdateUCNTimeInfo(G4double rtime);
  void UpdateUCNProperTimeInfo(G4double rtime);
  double GetGlobalTime() {return (double)GlobalTime;}
  double GetProperTime() {return (double)ProperTime;}
  double GetGlobalTime0() {return (double)GlobalTime0;}
  double GetProperTime0() {return (double)ProperTime0;}

	void SetPrecise(G4bool);
	inline G4bool GetPrecise();
  
 private:
  G4int                 originalTrackID;
  G4ParticleDefinition* particleDefinition;
  G4ThreeVector         originalPosition;
  G4ThreeVector         originalMomentum;
  G4double              originalEnergy;
  G4double              originalTime;

  long double GlobalTime;  //absolute time to be reserved 
  long double ProperTime;  //absolute proper time to be reserved 
  long double GlobalTime0;  //absolute time to be reserved 
  long double ProperTime0;  //absolute proper time to be reserved 

	G4bool precise;
};


extern G4Allocator<UCNTimeInformation> aTrackInformationAllocator;

inline void* UCNTimeInformation::operator new(size_t)
{ void* aTrackInfo;
  aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
  return aTrackInfo;
}

inline void UCNTimeInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator.FreeSingle((UCNTimeInformation*)aTrackInfo);}

inline G4bool UCNTimeInformation::GetPrecise(){
	return precise;
}

#endif
