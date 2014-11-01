//////////////////////////////////////////////////////////////////////
//  the class to make an event ended at the time user wants to stop
//  7.2014 Ryo katayama
//  9.2014 changed as it can get stoptime from its messenger by Ryo Katayama
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "UCNTimeStopStep.hh"
#include "UCNTimeStopStepMessenger.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4StepPoint.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleDefinition.hh"

#include "UCNPrimaryGeneratorAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "UCNDetectorConstruction.hh"
#include "UCNFieldSetup.hh"
#include "UCNPhysicsList.hh"

#include "UCNTimeInformation.hh"

#define neV (1.0e-9*eV)

UCNTimeStopStep* UCNTimeStopStep::theInstance = 0;

	UCNTimeStopStep::UCNTimeStopStep(const G4String& processName)
: G4VContinuousDiscreteProcess(processName),stoptime(1e9)
{
	theMessenger = new UCNTimeStopStepMessenger(this);
	theInstance = this;
	tree = new TTree("katayama","katayama");
	tree->Branch("x",&x,"x/D");
	tree->Branch("y",&y,"y/D");
	tree->Branch("z",&z,"z/D");
	tree->Branch("px",&px,"px/D");
	tree->Branch("py",&py,"py/D");
	tree->Branch("pz",&pz,"pz/D");
	tree->Branch("t",&t,"t/D");
	tree->Branch("spinx",&spinx,"spinx/D");
	tree->Branch("spiny",&spiny,"spiny/D");
	tree->Branch("spinz",&spinz,"spinz/D");
	tree->Branch("phase",&phase,"phase/D");
	tree->Branch("dph",&ddph,"dph/D");
	tree->Branch("dsy",&ddsy,"dsy/D");
	tree->Branch("evt",&evt,"t/I");
}

UCNTimeStopStep::~UCNTimeStopStep(){


	printf("%d\n",strlen(RootFileName.c_str()));

	char str[500];
	printf("textid=%d:edir=%d:rot=%d\n",textid,rot,edir);
	/*
		 if( strlen(RootFileName.c_str())==0 ){
		 strcpy(str,"test.root");
		 }else{
		 strcpy(str,RootFileName.c_str());
		 }
	 */

	//sprintf(str,"data_%d,%d,%0.1f,%0.1f,%0.1f.root",
	//		        (rot==1)?1:-1, (edir==1)?1:-1, b,w,V );
	sprintf(str,"./Repository/root/TextID_%d:%s:%s.root",
			textid,(edir==1)?"Eup":"Edown",
			(rot==1)?"Clockwise":"Counter-Clockwise");

	TFile *file = TFile::Open(str,"RECREATE");
	G4cout<<str<<G4endl;
	tree->Write();
	file->Close();

	delete theMessenger;
	delete tree;
} 

G4VParticleChange* UCNTimeStopStep::PostStepDoIt(const G4Track& aTrack, 
		const G4Step& aStep)
{
	aParticleChange.Initialize(aTrack);

		G4double time = aStep.GetPostStepPoint()->GetGlobalTime();

		if(fabs(stoptime-time)<0.1){

			G4EventManager* evtManager = G4EventManager::GetEventManager();
			G4Event* Evt = (G4Event*)evtManager->GetConstCurrentEvent();
			evt = Evt->GetEventID();

			G4ThreeVector pos = aStep.GetPostStepPoint()->GetPosition();
			G4ThreeVector mom = aStep.GetPostStepPoint()->GetMomentum();

			G4StepPoint * thePostPoint = aStep.GetPostStepPoint();
			G4int stepnr = aTrack.GetCurrentStepNumber();

			double Pol[3];
			double omega = -183.2471999767;
			double omegar = 3000./60.;// 3000 mm/60 mm
			double omegarel = omega*(2.41/1.91)*(3.33e-11/1e-6);
			double domega = 0.5*(omegarel*omegarel/(omega-omegar));
			domega=0; //For Non_rel

			omega+=domega;

			G4ThreeVector C  =  aStep.GetPostStepPoint()->GetPolarization();
			Pol[0]=C.x(); Pol[1]=C.y(); Pol[2]=C.z();

			double dph = TMath::ATan2(TMath::Cos(omega*time/1e9),TMath::Sin(omega*time/1e9)) - TMath::ATan2(C.y(),C.x()) ;
			double dsy = TMath::Cos(omega*time/1e9)-C.y() ;

			//fill the results
			x=pos.x(); y=pos.y(); z=pos.z();
			px=mom.x(); py=mom.y(); pz=mom.z();
			t=time/1e9;
			spinx=Pol[0]; spiny=Pol[1]; spinz=Pol[2];
			phase=TMath::ATan2(C.y(),C.x());

			ddph=dph;
			ddsy=dsy;

			tree->Fill();

			////////////////////////////////////
			G4RunManager* runManager = G4RunManager::GetRunManager();
			UCNPrimaryGeneratorAction* upg = (UCNPrimaryGeneratorAction*)runManager->GetUserPrimaryGeneratorAction();
			UCNDetectorConstruction* udc = (UCNDetectorConstruction*)runManager->GetUserDetectorConstruction();
			UCNPhysicsList* upl = (UCNPhysicsList*)runManager->GetUserPhysicsList();

			if(upg->GetRotation()) printf("Clockwise\n");
			else printf("Counter-clockwise\n");

			UCNFieldSetup *ufs = (UCNFieldSetup*)udc->GetFieldSetup();

			rot = (G4int)upg->GetRotation();
			edir = ufs->GetEFieldDirection();
			b    = upl->GetUCNMaterialBoundary()->GetBParameter();
			w    = upl->GetUCNMaterialBoundary()->GetWParameter();
			V    = upl->GetUCNMaterialBoundary()->GetVParameter();
			textid = upl->GetUCNMaterialBoundary()->GetTextID();

			if(edir) printf("EField is up\n");
			else printf("EField is down\n");

			printf("edir=%d\n",edir);

			printf("b,w,V=%f,%f,%f,\n",b,w,V);
			///////////////////////////////////

			aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
			G4cout<<"Final time = "<<time/1e9<<" sec"<<G4endl;

		}//end if fabs


	return G4VContinuousDiscreteProcess::PostStepDoIt(aTrack, aStep);
}


G4double UCNTimeStopStep::GetMeanFreePath(const G4Track& aTrack,
		G4double ,
		G4ForceCondition* condition)
{
	SetGPILSelection(CandidateForSelection);

	G4double time = aTrack.GetGlobalTime();

//if(time>stoptime*0.99999){//Before correcting
	if(time>stoptime*0.999999){
		*condition = Forced;
		G4double step_length = fabs(time-stoptime)*aTrack.GetVelocity()*0.1;
		return step_length;
	}//end if time>stoptime*0.9999

	return DBL_MAX;

}

G4VParticleChange* UCNTimeStopStep::AlongStepDoIt(
		const G4Track& aTrack,
		const G4Step&)
{
	aParticleChange.Initialize(aTrack);
	return &aParticleChange;
}

G4double UCNTimeStopStep::GetContinuousStepLimit(const G4Track &  aTrack,
		G4double  previousStepSize,
		G4double currentMinimumStep,
		G4double& )
{
	if ((previousStepSize == 0 ) && aTrack.GetCurrentStepNumber()-1) {
		return .1*micrometer;
	}
	return std::min(currentMinimumStep, 10*cm);
}


UCNTimeStopStep* UCNTimeStopStep::GetInstance()
{
	return theInstance;
}

void UCNTimeStopStep::SetStopTime(G4double StopTime){
	stoptime=StopTime;
	G4cout<<"stoptime is "<<stoptime<<G4endl;
}
