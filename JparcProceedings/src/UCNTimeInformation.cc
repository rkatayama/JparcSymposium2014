//////////////////////////////////////////////////////////////////////
//  the class to store the time with long double precision for UCN
//  7.2014 Masakazu Kurata
////////////////////////////////////////////////////////////////////////

#include "UCNTimeInformation.hh"
#include "G4ios.hh"

G4Allocator<UCNTimeInformation> aTrackInformationAllocator;

UCNTimeInformation::UCNTimeInformation()
{
    originalTrackID = 0;
    particleDefinition = 0;
    originalPosition = G4ThreeVector(0.,0.,0.);
    originalMomentum = G4ThreeVector(0.,0.,0.);
    originalEnergy = 0.;
    originalTime = 0.;
    GlobalTime =0.0;
    ProperTime=0.0;
    GlobalTime0 =0.0;
    ProperTime0=0.0;
}

UCNTimeInformation::UCNTimeInformation(G4Track* aTrack)
{
    originalTrackID = aTrack->GetTrackID();
    particleDefinition = aTrack->GetDefinition();
    originalPosition = aTrack->GetPosition();
    originalMomentum = aTrack->GetMomentum();
    originalEnergy = aTrack->GetTotalEnergy();
    originalTime = aTrack->GetGlobalTime();
    GlobalTime =0.0;
    ProperTime=0.0;
    GlobalTime0 =0.0;
    ProperTime0=0.0;
}

UCNTimeInformation::UCNTimeInformation(const UCNTimeInformation* aTrackInfo)
{
    originalTrackID = aTrackInfo->originalTrackID;
    particleDefinition = aTrackInfo->particleDefinition;
    originalPosition = aTrackInfo->originalPosition;
    originalMomentum = aTrackInfo->originalMomentum;
    originalEnergy = aTrackInfo->originalEnergy;
    originalTime = aTrackInfo->originalTime;
    GlobalTime =0.0;
    ProperTime=0.0;
    GlobalTime0 =0.0;
    ProperTime0=0.0;
}

UCNTimeInformation::~UCNTimeInformation(){
  GlobalTime=0.0;
  ProperTime=0.0;
  GlobalTime0 =0.0;
  ProperTime0=0.0;
}

void UCNTimeInformation::UpdateUCNTimeInfo(G4double rtime){
	//printf("UpdateUCNTimeInfo\n");
	GlobalTime0=GlobalTime;
  GlobalTime+=(long double)rtime;

  return;
}

void UCNTimeInformation::UpdateUCNProperTimeInfo(G4double rtime){
	//printf("UpdateUCNProperTimeInfo\n");
	ProperTime0=ProperTime;
  ProperTime+=(long double)rtime;

  return;
}

void UCNTimeInformation::SetPrecise(G4bool Precise){
	precise = Precise;
}
