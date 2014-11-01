
//////////////////////////////////////////////////////////////////////
//  Material Boundary Process for ultracold neutrons
//  9.9.04 peter fierlinger
////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------
// Change to calculate diffuse scattering using micro-roughness model 6,2013 Ryo Katayama and Dai Sakurai
// --------------------------------------------------------------

#include "G4ios.hh"
#include "UCNMaterialBoundary.hh"
#include "UCNMaterialBoundaryMessenger.hh"
#include "UCNShutterMessenger.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4StepPoint.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleDefinition.hh"

//add by katayama
#include <TRandom.h>/* header file for gRandom */
#include <TMath.h>  /* header file for TMath */
#include <TH1.h>    /* header file for 1-d histogram */
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>/* header file for TCanvas */
#include <TMath.h>  /* header file for TMath */
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TVector3.h"

double microroughness(double *x,double *par);
#define neV (1.0e-9*eV)

/**
 * Add a state (0 - closed, 1 - open) at a given time
 * keeps the array sorted so you don't have to
 */
void UCNShutterStates::AddState(int state, float time) {
	if (n >= UCN_SHUTTERS_STATES_MAX) {
		G4cerr << "Shutter states maximum reached" << G4endl;
	}
	// keep array in order
	int i, j;
	for (i=0; i<n; i++) {
		if (time < times[i]) {
			break;
		}
	}
	// shift tail of array
	for (j=n-1; j>=i; j--) {
		states[j+1] = states[j];
		times[j+1] = times[j];
	}
	// insert
	states[i] =  state;
	times[i] = time;
	n++;
}

/**
 * Get shutter state (0 - closed, 1 - open) at a given time
 * default value if not in definite state
 */
int UCNShutterStates::GetState(float time) {
	int state = UCN_SHUTTERS_DEFAULT_STATE;
	if (loop && n > 0) {
		// in case of loop, offset time (mod)
		float last_t = times[n-1];
		if (time > last_t) {
			int m = floor(time/last_t);
			time -= m * last_t;
			// set to last state in loop
			state = states[n-1];
		}
	}
	for (int i=0; i<n; i++) {
		if (time >= times[i]) {
			state = states[i];
		}
		else {
			break;
		}
	}
	return state;
}

void UCNShutterStates::SetLoop(int loop) {
	this->loop = loop;
}

UCNMaterialBoundary* UCNMaterialBoundary::theInstance = 0;

	UCNMaterialBoundary::UCNMaterialBoundary(const G4String& processName)
: G4VContinuousDiscreteProcess(processName)
	, out("killing.dat",std::ios::out)
{
	theMessenger = new UCNMaterialBoundaryMessenger(this);
	theShutterMessenger = new UCNShutterMessenger(this);
	theInstance = this;
	just_reflected = 0;
	fermipotdiff_shutter = 1e300;
	return_it = 0;
	sh_fermipot = 0;sh_spinflip = 0;sh_loss = 0;sh_diffusion = 0;sh_reflectivity = 0;
	sh_abscs = 0;sh_losscs = 0;sh_scatcs = 0;
	return_it = 0;
	useshutters = 0;

	G4cout<<"textid = "<<textid<<G4endl;
	G4cout<<"diffus_select = "<<diffus_select<<G4endl;


	ntuple = new TNtuple("ntuple","ntuple","RandomTheta:RandomPhi");

	///*
	c = 2.99792458*1.e8*1.e9;//nm/s
	hbar = 6.58211928*1.e-16*1.e9;//neV*s 
	h_plank = hbar*2*TMath::Pi();//neV*s
	m_n = 9.39565378*1.e8*1.e9/c/c;//neV/c^2 
	//*/

}

UCNMaterialBoundary::~UCNMaterialBoundary(){
	delete theMessenger;
	delete theShutterMessenger;
	//add by katayama
	delete ntuple;

} 


void UCNMaterialBoundary::AssignRoughness(G4String strstr){

	sscanf(strstr.c_str(),"%s %d",diffus_select,&textid);
	printf("diffus_select=%s,textid=%d\n",diffus_select,textid);

	//add by katayama
	int wawawa=0;
	char *STR = new char[200];

	switch(textid){

		case 0://DLC
			snprintf(STR,100,"./micro_roughness/b_0.6nm:w_22.0nm:V_220.0neV.txt");//Table has (b,w,V)
			b = 0.6;//nm
			w = 22.;//nm
			V = 250.0;//neV
			break;

		case 1://quartz1
			snprintf(STR,100,"./micro_roughness/b_0.7nm:w_22.0nm:V_94.0neV.txt");//Table has (b,w,V)
			b = 0.7;//nm
			w = 22.0;//nm
			V = 94.0;//neV
			break;

		case 2://quartz2
			snprintf(STR,100,"./micro_roughness/b_0.7nm:w_33.0nm:V_94.0neV.txt");//Table has (b,w,V)
			b = 0.7;//nm
			w = 33.0;//nm
			V = 94.0;//neV
			break;

		case 3://quartz3
			snprintf(STR,100,"./micro_roughness/b_0.7nm:w_45.0nm:V_94.0neV.txt");//Table has (b,w,V)
			b = 0.7;//nm
			w = 45.0;//nm
			V = 94.0;//neV
			break;

		case 4://quartz4
			snprintf(STR,100,"./micro_roughness/b_1.4nm:w_22.0nm:V_94.0neV.txt");//Table has (b,w,V)
			b = 1.4;//nm
			w = 22.0;//nm
			V = 94.0;//neV
			break;

		case 5://quartz5
			snprintf(STR,100,"./micro_roughness/b_1.4nm:w_33.0nm:V_94.0neV.txt");//Table has (b,w,V)
			b = 1.4;//nm
			w = 33.0;//nm
			V = 94.0;//neV
			break;

		case 6://quartz6
			snprintf(STR,100,"./micro_roughness/b_1.4nm:w_45.0nm:V_94.0neV.txt");//Table has (b,w,V)
			b = 1.4;//nm
			w = 45.0;//nm
			V = 94.0;//neV
			break;

		case 7://quartz7
			snprintf(STR,100,"./micro_roughness/b_2.1nm:w_22.0nm:V_94.0neV.txt");//Table has (b,w,V)
			b = 2.1;//nm
			w = 22.0;//nm
			V = 94.0;//neV
			break;

		case 8://quartz8
			snprintf(STR,100,"./micro_roughness/b_2.1nm:w_33.0nm:V_94.0neV.txt");//Table has (b,w,V)
			b = 2.1;//nm
			w = 33.0;//nm
			V = 94.0;//neV
			break;

		case 9://quartz9
			snprintf(STR,100,"./micro_roughness/b_2.1nm:w_45.0nm:V_94.0neV.txt");//Table has (b,w,V)
			b = 2.1;//nm
			w = 45.0;//nm
			V = 94.0;//neV
			break;

		default:
			snprintf(STR,100,"./micro_roughness/b_0.6nm:w_22.0nm:V_220.0neV.txt");
			//Table has (b,w,V)
			b = 0.6;//nm
			w = 22.;//nm
			V = 250.0;//neV
			break;

	}//end switch

	ifstream data(STR);//Table has (b,w,V)

	printf("STR=%s\n",STR);
	delete STR;

	if(data.fail()){
		G4cerr<<"Microrougness File does not exist.\n"<<G4endl;
		exit(0);
	}else{

	while(!data.eof()) {
		data >> Lcondition[wawawa]//lambda
			>> Tcondition[wawawa]//theta
			>> Pcondition[wawawa];//non specular probability
		wawawa++;
		if( wawawa%100000==0 )
		printf("wawawa=%d\n",wawawa);
	}

	data.close();
	}

	////////////////////////////////////////////

}//end AssignMicroRoughness






G4VParticleChange* UCNMaterialBoundary::PostStepDoIt(const G4Track& aTrack, 
		const G4Step& aStep)
{
	aParticleChange.Initialize(aTrack);
	verboselevel=0;//6;
	verboseLevel=0;
	G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
	G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
	G4String volnam1 = pPreStepPoint->GetPhysicalVolume()->GetName();
	G4String volnam2 = pPostStepPoint->GetPhysicalVolume()->GetName();
	if (verboselevel > 5){
		G4cout << "MATERIALBOUNDARYPROCESS: vol1 " << volnam1 << ", vol2 " << volnam2 << " " << aTrack.GetKineticEnergy()/neV<< G4endl;
	}

	// make sure to be at a geometrical boundary
	if (pPostStepPoint->GetStepStatus() != fGeomBoundary) {
		just_reflected = 0;
		return G4VContinuousDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	if (verboselevel > 3){
		G4cout << "MATERIALBOUNDARYPROCESS: vol1 " << volnam1 << ", vol2 " << volnam2 <<  G4endl;
		G4cout << aTrack.GetMomentumDirection() << G4endl;
		G4cout << "pos " << pPreStepPoint->GetPosition() << " " << pPostStepPoint->GetPosition() <<G4endl;
	}

	G4Material* Material1 = pPreStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial();
	G4Material* Material2 = pPostStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial();

	// check, if we are at a shutter volume an if it is open or closed
	// if we use shutters,
	if (useshutters == 1){
		int nr = 0, nr2 = 0;
		sscanf(volnam2, "Shutter%d", &nr);
		sscanf(volnam1, "Shutter%d", &nr2);
		G4double t2 =  aTrack.GetGlobalTime();
		if (nr){
			// retrieve state for shutter nr at time t2
			int state = shutter_states[nr].GetState(t2);

			if (verboseLevel> 3) {
				G4cout << "ucnshutter: set the state for shutter " << nr << " to " << state << G4endl;
			}
			if (state == 1) {
				// the shutter is open.
				if (verboseLevel> 1) {
					G4cout << "enter a shutter, the shutter is open  " << G4endl;
				}
				/////// change the material properties and remember the values until leaving
				G4MaterialPropertiesTable * T;
				T = Material2->GetMaterialPropertiesTable();

				if (T) { 
					G4MaterialPropertyVector* eff1 =T->GetProperty("FERMIPOT");
					if (eff1) {sh_fermipot = eff1->Value(1);}
					G4MaterialPropertyVector* eff2 =T->GetProperty("REFLECTIVITY");
					if (eff2) {sh_reflectivity = eff2->Value(1);}
					G4MaterialPropertyVector* eff3 =T->GetProperty("DIFFUSION");
					if (eff3) {sh_diffusion = eff3->Value(1);}
					G4MaterialPropertyVector* eff4 =T->GetProperty("SPINFLIP");
					if (eff4) {sh_spinflip = eff4->Value(1);}
					G4MaterialPropertyVector* eff5 =T->GetProperty("LOSS");
					if (eff5) {sh_loss = eff5->Value(1);}
					G4MaterialPropertyVector* eff6 =T->GetProperty("LOSSCS");
					if (eff6) {sh_losscs = eff6->Value(1);}
					G4MaterialPropertyVector* eff7 =T->GetProperty("ABSCS");
					if (eff7) {sh_abscs = eff7->Value(1);}
					G4MaterialPropertyVector* eff8 =T->GetProperty("SCATCS");
					if (eff8) {sh_scatcs = eff8->Value(1);}
					T->RemoveProperty("FERMIPOT");
					T->RemoveProperty("REFLECTIVITY");
					T->RemoveProperty("DIFFUSION");
					T->RemoveProperty("SPINFLIP");
					T->RemoveProperty("LOSS");
					T->RemoveProperty("LOSSCS");
					T->RemoveProperty("ABSCS");
					T->RemoveProperty("SCATCS");
					return_it = 1;
				}	
			}	  	  
			if (state == 0){
				// shutter is closed, dont change properties
				if (verboseLevel> 1) {
					G4cout << "the shutter is closed  " << G4endl;
				}
				return_it = 0;
			}

		}  
		else if (nr2){
			if (verboseLevel> 3) {
				G4cout << "UCNSHUTTER: leave a shutter volume " << G4endl;
			}
			if (return_it ==1){	
				const G4int NUM = 2;
				G4double PP[NUM] =  { 1.0, 1.0};
				G4MaterialPropertiesTable * T;
				T = Material1->GetMaterialPropertiesTable();
				G4double FE_POT[NUM] =         {sh_fermipot,sh_fermipot};      // neV
				G4double FE_SPINFLIP[NUM] =      {sh_spinflip,sh_spinflip};   // rel. per wall collision
				G4double FE_ETA[NUM] =          {sh_loss,sh_loss};        // loss coefficient W/V
				G4double FE_DIFFUS[NUM] =       {sh_diffusion,sh_diffusion};    // diffuse scattering probability
				G4double FE_REFLECTIVITY[NUM] =     {sh_reflectivity,sh_reflectivity};   // reflectivity, not used parameter
				G4double FE_ABSCS[NUM] =       {sh_abscs,sh_abscs};    // 1/v loss cross section at room temperature for Be
				G4double FE_LOSSCS[NUM] =       {sh_losscs,sh_losscs};    // loss cross section at room temperature for Be
				G4double FE_SCATCS[NUM] =       {sh_scatcs,sh_scatcs};    // (incoherent) "elastic" scattering cs
				T->AddProperty("REFLECTIVITY", PP, FE_REFLECTIVITY,      NUM);
				T->AddProperty("DIFFUSION",    PP, FE_DIFFUS,            NUM);
				T->AddProperty("FERMIPOT",     PP, FE_POT,               NUM);
				T->AddProperty("SPINFLIP",     PP, FE_SPINFLIP,          NUM);
				T->AddProperty("LOSS",         PP, FE_ETA              , NUM);
				T->AddProperty("LOSSCS",       PP, FE_LOSSCS           , NUM);
				T->AddProperty("ABSCS",        PP, FE_ABSCS            , NUM);
				T->AddProperty("SCATCS",       PP, FE_SCATCS           , NUM);
				Material1->SetMaterialPropertiesTable(T);


				fermipotdiff_shutter = sh_fermipot;

				sh_fermipot = 0;
				sh_spinflip = 0;
				sh_loss = 0;
				sh_diffusion = 0;
				sh_reflectivity = 0;
				sh_abscs = 0;
				sh_losscs = 0;
				sh_scatcs = 0;	
				return_it = 0;
			}
		} 
	}

	if (Material1 == Material2) return G4VContinuousDiscreteProcess::PostStepDoIt(aTrack, aStep);

	// get the normal to the surface
	G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();
	G4ThreeVector thePrePoint = pPreStepPoint->GetPosition();

	G4Navigator* theNavigator =
		G4TransportationManager::GetTransportationManager()->
		GetNavigatorForTracking();

	G4ThreeVector theLocalPoint = theNavigator->
		GetGlobalToLocalTransform().
		TransformPoint(theGlobalPoint);

	G4ThreeVector theLocalNormal;	// Normal points back into volume

	G4bool valid;
	theLocalNormal = -theNavigator->GetLocalExitNormal(&valid);

	if (!valid) {
		G4cout << "local normal: " << theLocalNormal <<G4endl;
	} 

	G4ThreeVector theGlobalNormal = theNavigator->GetLocalToGlobalTransform().
		TransformAxis(theLocalNormal);

	G4double vel = aTrack.GetVelocity();
	G4ThreeVector momdir = aTrack.GetMomentumDirection();
	G4ThreeVector mom = aTrack.GetMomentum();
	G4double tfakt = momdir.dot(theGlobalNormal);
	G4double momnorm = mom.dot(theGlobalNormal);
	G4double velnorm = vel * tfakt;
	G4double enormal = momnorm*momnorm/2./neutron_mass_c2/neV; 
	G4double energy = aTrack.GetKineticEnergy()/neV; 

	
	if( momdir.dot(-theGlobalNormal) <0 ){
	return G4VContinuousDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}


	G4MaterialPropertiesTable* aMaterialPropertiesTable;
	aMaterialPropertiesTable = Material2->GetMaterialPropertiesTable();

	G4double fermipot = 0.;
	G4double pdiffus = 0.;
	G4double pspinflip = 0.;
	G4double pupscatter = 0.;

	// properties in the new volume
	if (aMaterialPropertiesTable) { 
		G4MaterialPropertyVector* eff1 = aMaterialPropertiesTable->GetProperty("FERMIPOT");
		if (eff1) {fermipot = eff1->Value(1);}
	}

	if (aMaterialPropertiesTable) {
		G4MaterialPropertyVector* eff2 = aMaterialPropertiesTable->GetProperty("DIFFUSION");
		if (eff2) {
			//pdiffus = eff2->Value(1);
			pdiffus = 0;
		}
	}

	if (aMaterialPropertiesTable) {
		G4MaterialPropertyVector* eff3 = aMaterialPropertiesTable->GetProperty("SPINFLIP");
		if (eff3) {pspinflip = eff3->Value(1);}
	}

	if (aMaterialPropertiesTable) {
		G4MaterialPropertyVector* eff4 = aMaterialPropertiesTable->GetProperty("LOSS");
		if (eff4) { pupscatter = eff4->Value(1);}
	}

	// properties of the old volume
	G4MaterialPropertiesTable* aMaterialPropertiesTable2;
	aMaterialPropertiesTable2 = Material1->GetMaterialPropertiesTable();
	G4double fermipot_previous = 0.;
	if (aMaterialPropertiesTable2) { 
		G4MaterialPropertyVector* prev1 = aMaterialPropertiesTable2->GetProperty("FERMIPOT");
		if (prev1) {fermipot_previous = prev1->Value(1);}	
	}
	G4double pdiffus_previous = 0.;
	if (aMaterialPropertiesTable2) { 
		G4MaterialPropertyVector* prev2 = aMaterialPropertiesTable2->GetProperty("DIFFUSION");
		if (prev2) {pdiffus_previous = prev2->Value(1);}
	}
	G4double spinflip_previous = 0.;
	if (aMaterialPropertiesTable2) { 
		G4MaterialPropertyVector* prev3 = aMaterialPropertiesTable2->GetProperty("SPINFLIP");
		if (prev3) {spinflip_previous = prev3->Value(1);}
	}
	G4double loss_previous = 0.;
	if (aMaterialPropertiesTable2) { 
		G4MaterialPropertyVector* prev4 = aMaterialPropertiesTable2->GetProperty("LOSS");
		if (prev4) {loss_previous = prev4->Value(1);}
	}
	//// there is the possibility to use the inverted material properties,
	//   eg. for tubes in a vacuum, you dont have to specify inner and outer diameter,
	//   make it easier to build volumes.
	/*int inverted = 0;
		if (inverted == 1){
		G4double a1 = fermipot;
		G4double a2 = pdiffus;
		G4double a3 = pspinflip;
		G4double a4 = pupscatter;
		fermipot = fermipot_previous;
		pdiffus = pdiffus_previous;
		pspinflip = spinflip_previous;
		pupscatter = loss_previous;
		fermipot_previous = a1;
		pdiffus_previous  = a2;
		spinflip_previous = a3;
		loss_previous     = a4;
		}*/

	G4double fermipot_diff = fermipot - fermipot_previous;
	if (verboselevel > 2){ 
		G4cout << "MATERIALBOUNDARY: new fermipot " << fermipot << ", old " << fermipot_previous << G4endl;
		G4cout << "position " << thePrePoint << ", post " << theGlobalPoint << G4endl;
		G4cout << "energy " << energy << G4endl;
	}

	if (fermipotdiff_shutter < 1e300){
		if (verboselevel > 2){ 
			G4cout << "override by shutter " << fermipotdiff_shutter << G4endl;
		}
		fermipot_diff = fermipot - fermipot_previous + fermipotdiff_shutter;
		fermipotdiff_shutter = 1e300;
		//if(fermipot_diff==0) return G4VContinuousDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}




	//add by katayama
	G4double nv = aTrack.GetVelocity();//G4default mm/nsec
	nv*= 1e15;//nm/sec
	G4double thetaMR = TMath::ACos(momdir.dot(-theGlobalNormal));
	G4double NL = h_plank/(m_n*nv);//neutron lambda[nm]



	///// below critical velocity
	if (verboselevel > 5)G4cout << "enrm " << enormal << " " << fermipot_diff << " " << comparepot(enormal, fermipot_diff)<<G4endl;

	/**************************************/
	if (comparepot(enormal, fermipot_diff) == 0){             // reflect from surface
		if (verboselevel >1) G4cout << "MATERIALBOUNDARY: reflect " << G4endl;

		//it add by katayama so that we can obtain the effect of diffuse reflection.
		/*

		/////// losses
		if (loss(pupscatter, velnorm, fermipot_diff) == 1){             // loss on reflection
			if (verboselevel >1) G4cout << "MATERIALBOUNDARY: loss on surface " << G4endl;
			// kill it.
			aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
		}//end if loss(...) == 1

		/////// spinflips
		if (spinflip(pspinflip) == 1) {
			if (verboselevel >1) G4cout << "MATERIALBOUNDARY: spinflip " << G4endl;
			G4ThreeVector spin = aTrack.GetPolarization();
			spin *= -1;
			aParticleChange.ProposePolarization(spin);
		}//end if spinflip(...) == 1

		*/

		////// reflect it

		//G4ThreeVector ref = reflect(0, pdiffus, momdir, theGlobalNormal);//default
		//add by katayama 
		Calc_Pdiffus( &pdiffus, NL, thetaMR );
		G4ThreeVector ref = reflect(0, pdiffus, momdir, theGlobalNormal,NL,b,w,V,thetaMR);

		ref = ref.unit();
		aParticleChange.ProposeMomentumDirection(ref);
		just_reflected = 1;

	}//end if comparepot(...) == 0
	/**************************************/

	/////// above critical velocity
	else {                                      // transmit material

		if (just_reflected == 0){

			// if it is faster than the crticial velocity, there is a probability to be still
			// reflected. this formula is (only) valid for low loss materials
			G4double refl2 = reflectivity(0, fermipot_diff,enormal);

			if (verboselevel > 5)G4cout << "jr " << just_reflected << " " << G4endl; 
			if (verboselevel >1) G4cout << "MATERIALBOUNDARY: reflectivity " << refl2 << G4endl;

			if (G4UniformRand() < refl2){ 

				// Calc_Pdiffus( &pdiffus, vel, theGlobalNormal ); 
				// must calculate vel*=1e15 and thetaMR which is TMath::ACos(momdir.dot(-theGlobalNormal))

				////// reflect it
				//G4ThreeVector ref = reflect(0, pdiffus, momdir, theGlobalNormal);
				//add by katayama 
				Calc_Pdiffus( &pdiffus, NL, thetaMR );
				G4ThreeVector ref = reflect(0, pdiffus, momdir, theGlobalNormal,NL,b,w,V,thetaMR);

				aParticleChange.ProposeMomentumDirection(ref);
				// set the just_reflected variable 1
				just_reflected = 1;

			} else {

				// --- transmission because it is faster than the critical velocity

				// --- kinetic energy in the new media
				G4double enew = transmit(fermipot_diff, energy);

				// --- change of the normal momentum component
				//     p = sqrt(2*m*Ekin)
				G4double m = -sqrt(momnorm*momnorm - neutron_mass_c2*2.*fermipot_diff*neV);

				// --- momentum direction in new media
				G4ThreeVector ref = mom - (momnorm-m)*theGlobalNormal;
				aParticleChange.ProposeMomentumDirection(ref.unit());
				aParticleChange.ProposeEnergy(enew*neV);

				if (verboselevel >1){ 
					G4cout << "energy " << energy << ", enrm " << enormal << ", fpdiff " << fermipot_diff << ", enew " << enew << G4endl;
					G4cout << "MATERIALBOUNDARY: transmitt and set the new energy " << aParticleChange.GetEnergy()/neV << " ekin " << G4endl;
				}
				if (verboselevel > 5)G4cout<<volnam2<<G4endl; 


			}//end else (G4UniformRand() < refl2)

		}//end if justreflected == 0

		////////////////////////////////////
		else if (just_reflected == 1){

			if (verboselevel >1) G4cout << "justreflected " << G4endl;
			just_reflected = 0 ;

		}//end else if just relfected == 1
		////////////////////////////////////

	}//end else (comparepot(enormal, fermipot_diff) == 0)



	return G4VContinuousDiscreteProcess::PostStepDoIt(aTrack, aStep);
}


G4double UCNMaterialBoundary::GetMeanFreePath(const G4Track&,
		G4double ,
		G4ForceCondition* condition)
{
	SetGPILSelection(CandidateForSelection);
	*condition = Forced;
	return DBL_MAX;
}

G4VParticleChange* UCNMaterialBoundary::AlongStepDoIt(
		const G4Track& aTrack,
		const G4Step&)
{
	aParticleChange.Initialize(aTrack);
	return &aParticleChange;
}

G4double UCNMaterialBoundary::GetContinuousStepLimit(const G4Track &  aTrack,
		G4double  previousStepSize,
		G4double currentMinimumStep,
		G4double& )
{
	if ((previousStepSize == 0 || just_reflected==1) && aTrack.GetCurrentStepNumber()-1) {
		return .1*micrometer;
	}
	return std::min(currentMinimumStep, 10*cm);
}

int UCNMaterialBoundary::comparepot(G4double energy, G4double fermipot_diff)
{
	return (energy > fermipot_diff);    
}


int UCNMaterialBoundary::loss(G4double coefficient, G4double velnorm, G4double
		fermipot){

	/// the surface roughness is not taken into account here, 
	// one could use e.g. ultracold neutrons, r.golub, p.35,
	// where mu is increased by roughness parameters sigma and
	// omega, which are related to the height and the distance of
	// "disturbances" on the surface 

	G4double v_bound = sqrt(2.*fermipot*neV/neutron_mass_c2*c_squared);
	G4double v_ratio = velnorm/v_bound;
	G4double loss_von_E = (2*coefficient*v_ratio)/(sqrt(1-(v_ratio*v_ratio)));

	return (G4UniformRand() <= fabs(loss_von_E));
}


int UCNMaterialBoundary::spinflip(G4double coefficient){
	G4double rnd_loss = G4UniformRand();
	if (rnd_loss <= coefficient){
		//G4cout << "MATERIALBOUNDARY: spinflip " << rnd_loss << G4endl;
		return 1;
	}
	else{
		return 0;
	}
}

G4double UCNMaterialBoundary::reflectivity(int, G4double fpot, G4double enormal){

	G4double r = (sqrt(enormal) - sqrt(enormal - fpot)) /
		(sqrt(enormal) + sqrt(enormal - fpot));

	return r*r;
} 


G4ThreeVector UCNMaterialBoundary::reflect(int, G4double coefficient, 
		G4ThreeVector momdir, G4ThreeVector localnormal){
	// reflect specular
	G4double tfakt = momdir.dot(localnormal);
	G4ThreeVector reflected = momdir - 2 * localnormal * tfakt;
	if (verboselevel > 5)G4cout << "specular reflected << " << reflected << G4endl;

	// reflect diffuse 
	if (reflected==momdir || G4UniformRand() < coefficient) {
		G4ThreeVector diffus = cos_diff(localnormal);
		if (verboselevel > 5)G4cout << "diffus localnormal " << localnormal << ", " << diffus << G4endl;
		return diffus;
	}
	return reflected;
}

G4double UCNMaterialBoundary::transmit(G4double fermipot, G4double energy){

	return energy - fermipot;

}


G4ThreeVector UCNMaterialBoundary::cos_diff(G4ThreeVector localnormal)
{
	G4ThreeVector momentum;
	// cosine distribution - Lambert's law
	momentum.setRThetaPhi(1., acos(sqrt(G4UniformRand())), 2.*pi*G4UniformRand());
	momentum.rotateUz(localnormal);
	if(momentum*localnormal<0) {momentum*=-1;G4cout << "!" << G4endl;}
	return momentum;
}

UCNMaterialBoundary* UCNMaterialBoundary::GetInstance()
{
	return theInstance;
}

void UCNMaterialBoundary::SetFermiPotDiff(G4double state)
{
	fermipotdiff_shutter = state;
}

void UCNMaterialBoundary::setVerbose(G4int level)
{
	verboselevel = level ;
}

void UCNMaterialBoundary::SetShutterClose(G4String newval){
	int nr = 0;
	float ti = 0.;
	sscanf(newval, "%d %f", &nr, &ti); // shutternr time
	if (nr >= UCN_SHUTTERS_MAX || nr < 0) {
		G4cerr << "invalid shutter number " << nr << G4endl;
		return;
	}
	G4cout << "ucnshutter: setshutterclose (shutter " << nr
		<< ", t=" << ti << G4endl;
	shutter_states[nr].AddState(0, ti*second);
}

void UCNMaterialBoundary::SetShutterOpen(G4String newval){
	int nr = 0;
	float ti = 0.;
	sscanf(newval, "%d %f", &nr, &ti); // shutternr time
	if (nr >= UCN_SHUTTERS_MAX || nr < 0) {
		G4cerr << "invalid shutter number " << nr << G4endl;
		return;
	}
	G4cout << "ucnshutter: setshutteropen (shutter " << nr
		<< ", t=" << ti << G4endl;
	shutter_states[nr].AddState(1, ti*second);
}

void UCNMaterialBoundary::SetShutterLoop(G4String newval){
	int nr = 0;
	int loop = 0;
	sscanf(newval, "%d %d", &nr, &loop);
	if (nr >= UCN_SHUTTERS_MAX || nr < 0) {
		G4cerr << "invalid shutter number " << nr << G4endl;
		return;
	}
	G4cout << "ucnshutter: setshutter loop " << loop
		<< " (shutter " << nr << ")" << G4endl;
	shutter_states[nr].SetLoop(loop);
}

void UCNMaterialBoundary::SetUseShutters(G4String newval){
	if (newval == "1"){
		useshutters = 1;
	}

}
void UCNMaterialBoundary::SetShutterVerbose(G4String newval){
	if (newval == "1"){
		verboselevel = 1;
	}

}


//add by katayama
G4double UCNMaterialBoundary::Calc_Pdiffus( G4double *coefficient, G4double NL, G4double thetaMR ){

	int i = 0;

	if(strcmp(diffus_select,"lambert")==0){
		//printf("pdiffus=%f\n",*coefficient);
		return 0;
	}


	while(true){
		//NL = h_plank/(m_n*nv);//neutron lambda[nm]
		if ( Lcondition[i] < NL && Tcondition[i] > thetaMR ){
			*coefficient = (Pcondition[i]-Pcondition[i-1])
				/(Tcondition[i]-Tcondition[i-1])
				*(thetaMR-Tcondition[i-1])
				+Pcondition[i-1];
			/*
				 if(II==0){
				 G4cout<<"i am in Calc_Pdiffus function"<<G4endl;
				 G4cout<<"coefficient="<<*coefficient<<G4endl;
				 G4cout<<"coefficient2="<<
				 (Pcondition[i]-Pcondition[i-1])/(Tcondition[i]-Tcondition[i-1])
			 *(thetaMR-Tcondition[i-1]) + Pcondition[i-1]<<G4endl;
			 G4cout<<"(Pcondition[i]-Pcondition[i-1])"<<"/(Tcondition[i]-Tcondition[i-1])"<<
			 "*(thetaMR-Tcondition[i-1])"<<"+Pcondition[i-1]"<<G4endl;
			 printf("Pconditon[i]=%e\n",Pcondition[i]);
			 printf("Pconditon[i-1]=%e\n",Pcondition[i-1]);
			 printf("thetaMR=%e\n",thetaMR);
			 printf("NL=%e\n",NL);
			 printf("Tconditon[i]=%e\n",Tcondition[i]);
			 printf("Tconditon[i-1]=%e\n",Tcondition[i-1]);
			 printf("i=%d\n",i);
			 }
			//*/
			break;
		}
		i++;
	}

}//end Calc_Pdiffus


//add by katayama
double microroughness(double *x,double *par)
{
	double pi = TMath::Pi();

	double theta = x[0]; //final theta angle [rad] <- output1 (rad)
	double phi = x[1];  //final phi angle [rad] <- output2 (rad)
	double lambda = par[0];//corresponding 1/|velocity| (nm)
	double theta_i = par[1];//incident angle(rad)
	double b = par[2];
	double w = par[3];
	double potential = par[4]; //<- fermi potential(neV)
	double k = 2.0*pi/lambda; // wavenumber [nm-1]
	double kc = 6.946916214e-3 * sqrt(potential); //critical wavenumber [nm-1]

	double q2 = k*k * (sin(theta_i)*sin(theta_i) + sin(theta)*sin(theta) - 2.0*sin(theta_i)*sin(theta)*cos(phi)); //squared horizontal momentum transfer
	double sq = exp(-0.5*q2*w*w); //Fourier transform of roughness correlation function
	double factor = (kc*kc*kc*kc * b*b*w*w)/ (8.0*pi * cos(theta_i));

	double tmp;
	double si, si2; //Fresnel transmission factor (incident angle)
	double s, s2; //Fresnel transmission factor (finale angle)

	if (cos(theta_i) >= kc/k){
		tmp = sqrt(cos(theta_i)*cos(theta_i) - (kc*kc)/(k*k));
		si = 2.0*cos(theta_i) / (cos(theta_i) + tmp);
		si2 = si*si;
	}
	else {
		si2 = 4.0*k*k*cos(theta_i)*cos(theta_i)/(kc*kc);
	}

	if (cos(theta) >= kc/k){
		tmp = sqrt(cos(theta)*cos(theta) - (kc*kc)/(k*k));
		s = 2.0*cos(theta) / (cos(theta) + tmp);
		s2 = s*s;
	}
	else {
		s2 = 4.0*k*k*cos(theta)*cos(theta)/(kc*kc);
	}


	double func = factor * si2 * s2 * sq * sin(theta); // sin(theta) is for theta-phi plain plot
	// if this funcition is integrated, this is probability of non specular reflection.

	return func; //<- the angle distribution of non specular reflectivity
}//end microroughness



//add by katayama
G4ThreeVector UCNMaterialBoundary::reflect(int, G4double coefficient, 
		G4ThreeVector momdir, G4ThreeVector localnormal,G4double NL, G4double b, G4double w, G4double fermipot_diff,G4double thetaMR){

	// reflect specular
	G4double tfakt = momdir.dot(localnormal);
	G4ThreeVector reflected = momdir - 2 * localnormal * tfakt;
	if (verboselevel > 5)G4cout << "specular reflected << " << reflected << G4endl;


	// reflect diffuse

	//add by katayama
	if(strcmp(diffus_select,"lambert")==0){
		//printf("lambert lambert.\n");

		if (reflected==momdir || G4UniformRand() < coefficient) {
			G4ThreeVector diffus = cos_diff(localnormal);
			return diffus;
		}

	}
	else if(strcmp(diffus_select,"microroughness")==0){

		//printf("microroughness\n");

		//sakurai
		if ( G4UniformRand() < coefficient ){

			double pi = TMath::Pi();

			TF2 *f = new TF2("f",microroughness,0.,pi/2.,-pi,pi,5);

			f->SetParameters(NL,thetaMR,b,w,fermipot_diff);
			f->GetRandom2(randomtheta,randomphi);

			ntuple->Fill(randomtheta,randomphi);

			//G4cout<<"========================="<<G4endl;
			//G4cout<<"randomtheta="<<randomtheta<<G4endl;
			//G4cout<<"randomphi="<<randomphi<<G4endl;

			delete f;

			G4ThreeVector diffus;

			G4ThreeVector ez = localnormal;	
			G4ThreeVector ex = momdir - localnormal * tfakt;

			if(ex.mag()==0){

				// First, we consider the following situation: using ez = (Ex,Ey,Ez),
				// we make a "ex vector" such that ez.Dot(ex) = 0 and ex.Mag()!=0.
				// If we try to choose the condition as ex = (0,-Ez,Ey),
				// we can seem to easily achieve them.
				// However, there is a problem with this choice.
				// When ez is (1,0,0), ex vector will be zero vector.
				// In order to avoid its problem, we need the below treatment.

				double Z[3],Zmin;
				int Nmin=0;
				int IDmin=1;

				Z[0]=ez.x(); Z[1]=ez.y(); Z[2]=ez.z();

				Zmin = fabs(Z[0]);

				for(int l=0;l<3;l++){

					if( Zmin>fabs(Z[l]) ){
						Zmin = fabs(Z[l]);
					}
					//It is note that theoretically, the least of Zmin is 0.

				}//end for


				for(int l=0;l<3;l++){

					if( Zmin==fabs(Z[l]) ){
						Nmin++;
						IDmin*=(l+1);
					}

				}//end for


				printf("Nmin=%d\n",Nmin);
				printf("IDmin=%d\n",IDmin);

				switch(IDmin){

					case 1://for only x
						ex.set(0,-Z[2],Z[1]);
						//Ofcourse, it will be perpendicular to ez.
						break;
					case 2://for x=y
						ex.set(0,-Z[2],Z[1]);
						//If both ez.x() and ez.y() equals to 0, ex.SetXYZ(0,-1,0) can be perpendicular to ez(0,0,1).
						break;
					case 3://for x=z
						ex.set(0,-Z[2],Z[1]);
						//If both ez.x() and ez.z() equals to 0, ex.SetXYZ(0,0,1) can be perpendicular to ez(0,1,0).
						break;
					case 6://for y=z or x=y=z

						if(Nmin==2){//for y=z
							ex.set(-Z[1],Z[0],0);
							//Such ex will be perpendicular to ez if both y and z is 0.
						}else if(Nmin==3){//for x=y=z.
							//There is no constraint in this situation because "x=y=z=0" never happen.
							ex.set(0,-Z[2],Z[1]);
						}

						break;
					default:
						printf("*********\n");
						break;

				}//end switch


			}//end if ex.mag()==0


			//G4cout<<"ex="<<ex.mag()<<G4endl;
			ex = ex/ex.mag();

			G4ThreeVector ey = localnormal.cross(ex);
			//G4cout<<"ey="<<ey.mag()<<G4endl;
			ey = ey/ey.mag();



			G4double dx ;
			G4double dy ;
			G4double dz ;

			dx = momdir.dot(ex);
			dy = momdir.dot(ey);

			//G4double newdz = TMath::ASin(dz/sin(thetaMR)) + randomphi;
			//G4double newdy = TMath::ACos(dy/sin(thetaMR)) + randomphi;

			//dx = cos(pi - randomtheta);
			//dy = sin(pi - randomtheta)*cos(newdy);
			//dz = sin(pi - randomtheta)*sin(newdz);

			G4double newdx = TMath::ASin(dx/sin(thetaMR)) + randomphi;
			G4double newdy = TMath::ACos(dy/sin(thetaMR)) + randomphi;

			dx = sin(pi - randomtheta)*sin(newdx);
			dy = sin(pi - randomtheta)*cos(newdy);
			dz = cos(pi - randomtheta);

			dz*=-1;

			diffus += dx*ex;
			diffus += dy*ey;
			diffus += dz*ez;


			//G4cout<<"This is diffus"<<G4endl;

			/*
				 G4cout<<"(ex.x,ex.y,ex.x)=("<<ex.getX()<<","<<ex.getY()<<","<<ex.getZ()<<")"<<G4endl;
				 G4cout<<"(ey.x,ey.y,ey.x)=("<<ey.getX()<<","<<ey.getY()<<","<<ey.getZ()<<")"<<G4endl;
				 G4cout<<"(ez.x,ez.y,ez.x)=("<<ez.getX()<<","<<ez.getY()<<","<<ez.getZ()<<")"<<G4endl;
				 G4cout<<"(dif.x,dif.y,dif.z)=("<<diffus.getX()<<","<<diffus.getY()<<","<<diffus.getZ()<<")"<<G4endl;
			 */

			return diffus;

		}//end if G4UniformRand() < coefficient

	}//end if strcmp(diffus_select,"Micro Roughness")

	return reflected;

}//end reflect



