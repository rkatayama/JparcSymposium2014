//
//   UCN detector construction, 9.9.04, peter fierlinger 
// 

#include "UCNDetectorConstruction.hh"
#include "UCNDetectorMessenger.hh"

#include "UCNFieldSetup.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UserLimits.hh"
#include "G4ios.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4Torus.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4SDManager.hh"
/////////////////////////////////////////////////////////////////////////////
//
//

UCNDetectorConstruction::UCNDetectorConstruction()
 : solidWorld(0), logicWorld(0), physiWorld(0),
   solidAbsorber(0),logicAbsorber(0), physiAbsorber(0),
   fEmFieldSetup(0),
   AbsorberMaterial(0),worldchanged(false), WorldMaterial(0)
{
  // default parameter values of the calorimeter

  WorldSizeR = 150.*cm;
  //WorldSizeZ = 3000.*mm;
  WorldSizeZ = 150.*cm;

  AbsorberThickness = 0.0*mm;
  
  AbsorberRadius   = 2.*cm;
  zAbsorber = -20.*cm ;
  
  
  Shutter1x      = 0.*cm;
  Shutter1y      = 3.*cm;
  Shutter1z      = -40.*cm;
  Shutter1Radius = 4.*cm;
  Shutter1Thickness = 0.*cm;
  
  // create commands for interactive definition of the calorimeter  

  detectorMessenger = new UCNDetectorMessenger(this);
  
  DefineMaterials();

  fEmFieldSetup = new UCNFieldSetup() ;
}

//////////////////////////////////////////////////////////////////////////
//
//

UCNDetectorConstruction::~UCNDetectorConstruction()
{ 
  delete detectorMessenger;
  delete fEmFieldSetup ;
}

//////////////////////////////////////////////////////////////////////////
//
//

G4VPhysicalVolume* UCNDetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//////////////////////////////////////////////////////////////////////////////
//
//

void UCNDetectorConstruction::DefineMaterials()
{ 
  #include "UCNDetectorMaterials.icc"
  //This function illustrates the possible ways to define materials
  
}

/////////////////////////////////////////////////////////////////////////
//
//
  
G4VPhysicalVolume* UCNDetectorConstruction::ConstructCalorimeter()
{      
	// complete the Calor parameters definition and Print 

	ComputeCalorParameters();
	//PrintCalorParameters();

	// Cleanup old geometry

	if (physiWorld)
	{
		G4GeometryManager::GetInstance()->OpenGeometry();
		G4PhysicalVolumeStore::GetInstance()->Clean();
		G4LogicalVolumeStore::GetInstance()->Clean();
		G4SolidStore::GetInstance()->Clean();
	}

	G4double X = 200.*mm; // mm 
	G4double Y = 200.*mm; // mm 
	G4double Z = 12000.*mm; // mm


	// World
	solidWorld = new G4Box("World",X,Y,Z);       //its name and size

	logicWorld = new G4LogicalVolume(solidWorld,		//its solid
			WorldMaterial,	//its material
			"World");		//its name

	physiWorld = new G4PVPlacement(0,			//no rotation
			G4ThreeVector(),	//at (0,0,0)
			"World",		//its name
			logicWorld,		//its logical volume
			0,			//its mother  volume
			false,			//no boolean operation
			0);			//copy number                             

	//G4UserLimits * theWorldUserLimits=new G4UserLimits(20.*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(10.0*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(5.0*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(1.0*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(0.5*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(0.3*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(0.1*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(0.075*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(0.06*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(0.05*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(0.04*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(0.03*mm);
	//G4UserLimits * theWorldUserLimits=new G4UserLimits(0.015*mm);
	G4UserLimits * theWorldUserLimits=new G4UserLimits(0.01*mm);
	logicWorld->SetUserLimits(theWorldUserLimits);

	logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
	
	G4VisAttributes * greenAndWire
	  = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
	G4VisAttributes * cyanAndWire
	  = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	G4VisAttributes * cyan2AndWire
	  = new G4VisAttributes(G4Colour(0.0,0.5,1.0));
	G4VisAttributes * cyanAndSolid
	  = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	G4VisAttributes * redAndSolid
	  = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
	G4VisAttributes * blueAndSolid
	  = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	G4VisAttributes * blueAndWire
	  = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	G4VisAttributes * yellowAndSolid
	  = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	G4VisAttributes * yellowAndWire
	  = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	G4VisAttributes * orangeAndWire
	  = new G4VisAttributes(G4Colour(1.0,0.75,0.0));
	G4VisAttributes * redAndWire
	  = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
	G4VisAttributes * grayAndWire
	  = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	G4VisAttributes * grayAndSolid
	  = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	G4VisAttributes * whiteAndSolid
	  = new G4VisAttributes(G4Colour(1.,1.,1.));
	
	redAndSolid->SetVisibility(true);
	blueAndSolid->SetVisibility(true);
	yellowAndSolid->SetVisibility(true);
	
	greenAndWire->SetVisibility(true);
	cyanAndWire->SetVisibility(true);
	orangeAndWire->SetVisibility(true);
	redAndWire->SetVisibility(true);
	grayAndWire->SetVisibility(true);
	
	redAndSolid->SetForceSolid(true);
	blueAndSolid->SetForceSolid(true);
	cyanAndSolid->SetForceSolid(true);
	yellowAndSolid->SetForceSolid(true);
	grayAndSolid->SetForceSolid(true);
	whiteAndSolid->SetForceSolid(true);
	
	yellowAndWire->SetForceWireframe(true);
	greenAndWire->SetForceWireframe(true);
	cyanAndWire->SetForceWireframe(true);
	grayAndWire->SetForceWireframe(true);
	redAndWire->SetForceWireframe(true);
	

	/////////////////////////////////////////////////////////////////////////////////////////////

	bool BOXGUIDE;
	bool Content;
	bool STRAGECELL;

	BOXGUIDE = false;
	Content = false;
	STRAGECELL = true;


	G4RotationMatrix rot;

	/*

	if(BOXGUIDE){

		double L = 5000.*mm; // mm 
		double S = 50.*mm; // mm 
		double T = 10.*mm; // mm 

		G4Box* xtop = new G4Box("xtop",T,S,L);
		G4Box* xbottom = new G4Box("xbottom",T,S,L);
		G4Box* ytop = new G4Box("ytop",S,T,L);
		G4Box* ybottom = new G4Box("ybottom",S,T,L);


		G4LogicalVolume* logicxtop = new G4LogicalVolume(xtop,
				Shutter1Material,"xtop");
		G4LogicalVolume* logicxbottom = new G4LogicalVolume(xbottom,
				Shutter1Material,"xbottom");
		G4LogicalVolume* logicytop = new G4LogicalVolume(ytop,
				Shutter1Material,"ytop");
		G4LogicalVolume* logicybottom = new G4LogicalVolume(ybottom,
				Shutter1Material,"ybottom");


		G4UserLimits * theSolidUserLimits=new G4UserLimits(1.*mm);
		//G4UserLimits * theSolidUserLimits=new G4UserLimits(0.005*mm);
		logicxtop->SetUserLimits(theSolidUserLimits);
		logicxbottom->SetUserLimits(theSolidUserLimits);
		logicytop->SetUserLimits(theSolidUserLimits);
		logicybottom->SetUserLimits(theSolidUserLimits);


		G4double phi, top, back;
		phi = 0;        top = 50.0*mm + T;		back = L;

		rot.rotateX(phi);

		G4VPhysicalVolume* physixtop;
		physixtop = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(top,0,back)),
				"xtop", logicxtop, physiWorld, false, 0);

		G4VPhysicalVolume* physixbottom;
		physixbottom = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(-top,0,back)),
				"xbottom", logicxbottom, physiWorld, false, 0);

		G4VPhysicalVolume* physiytop;
		physiytop = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(0,top,back)),
				"ytop", logicytop, physiWorld, false, 0);

		G4VPhysicalVolume* physiybottom;
		physiybottom = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(0,-top,back)),
				"ybottom", logicybottom, physiWorld, false, 0);



		if(Content){

			G4Box* content = new G4Box("content",S,S,L);//

			G4LogicalVolume* logiccontent = new G4LogicalVolume(content,//
					Vacuum,"content");//

			G4UserLimits * VacuumUserLimits=new G4UserLimits(1.*mm);//
			logiccontent->SetUserLimits(VacuumUserLimits);//

			G4VPhysicalVolume* physiycontent;//
			physiycontent = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(0,0,back)),
					"content", logiccontent, physiWorld, false, 0);//

		}


	}//end if BOXGUIDE
	*/



	if(STRAGECELL){

		// the cylinder of the storage volume  

		G4double cylinderradius = 70.*mm; // mm 
		G4double cylinderinnerradius = 60.*mm; // mm 
		//G4double cylinderthickness = 120.*mm; //mm
		G4double cylinderthickness = 100.*mm; //mm

		// cylinder
		G4Tubs* solidcylinder = new G4Tubs("cylinder",cylinderinnerradius,cylinderradius,
				cylinderthickness/2.,0.,twopi); 

		G4LogicalVolume*   logiccylinder = new G4LogicalVolume(solidcylinder, 
				Shutter1Material,"cylinder");   
		//G4UserLimits * CylinderLimits=new G4UserLimits(1.*mm);//
		//logiccylinder->SetUserLimits(CylinderLimits);//


		G4double phi, x, y, z;
		phi = 0;        x = 0.0*mm;       y = 0.0*mm;   z = 0.0*mm;

		rot.rotateX(phi);
		G4VPhysicalVolume* physicylinder;
		physicylinder = new G4PVPlacement(
				G4Transform3D(rot,G4ThreeVector(x,y,z)),
				"storagecylinder", logiccylinder, 
				physiWorld, false, 0);
	
		logiccylinder->SetVisAttributes(grayAndSolid);

		////////////////////////
		G4double contentradius = 60.*mm; // mm 
		G4double contentinnerradius = 0.*mm; // mm 
		//G4double contentthickness = 120.*mm; //mm
		G4double contentthickness = 100.*mm; //mm

		G4Tubs* contentcylinder = new G4Tubs("content",contentinnerradius,contentradius,
				cylinderthickness/2.,0.,twopi); 

		G4LogicalVolume*   logicalcontent = new G4LogicalVolume(contentcylinder, 
			        WorldMaterial,"content");   
		logicalcontent->SetUserLimits(theWorldUserLimits);//

		G4VPhysicalVolume* physicscontent;
		physicylinder = new G4PVPlacement(
				G4Transform3D(rot,G4ThreeVector(x,y,z)),
				"content", logicalcontent, 
				physiWorld, false, 0);
			
			////////////////////



		G4double diskradius = 70.*mm; // mm 
		G4double diskinnerradius = 0.*mm; // mm 
		G4double diskthickness = 10.*mm; //mm
		// cylinder
		G4Tubs* soliddisk = new G4Tubs("disk",diskinnerradius,diskradius,
				diskthickness/2.,0.,twopi); 
		G4LogicalVolume*   logicdisk = new G4LogicalVolume(soliddisk, 
				Shutter1Material,"disk");   

		G4UserLimits * DiskLimits=new G4UserLimits(1.*mm);//
		logicdisk->SetUserLimits(DiskLimits);//
		
		
		G4VPhysicalVolume* physicydisk;
		G4VPhysicalVolume* physicydisk2;
		G4RotationMatrix rot; rot.rotateX(0.);
		/*
		physicydisk = new G4PVPlacement(
				G4Transform3D(rot,G4ThreeVector(0.,0.,65.*mm)),
				"storagedisk", logicdisk, physiWorld, false, 0);
		physicydisk2 = new G4PVPlacement(
				G4Transform3D(rot,G4ThreeVector(0.,0.,-65.*mm)),
				"storagedisk2", logicdisk, physiWorld, false, 0);		
		 */

		physicydisk = new G4PVPlacement(
				G4Transform3D(rot,G4ThreeVector(0.,0.,55.*mm)),
				"storagedisk", logicdisk, physiWorld, false, 0);
		physicydisk2 = new G4PVPlacement(
				G4Transform3D(rot,G4ThreeVector(0.,0.,-55.*mm)),
				"storagedisk2", logicdisk, physiWorld, false, 0);


	}//end if STRAGECELL


	/////////////////////////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------------------------------------------------
  return physiWorld;
}

///////////////////////////////////////////////////////////////////////////
//
//

void UCNDetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
	// get the pointer to the material table
	const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

	// search the material by its name   
	G4Material* pttoMaterial;
	for (size_t J=0 ; J<theMaterialTable->size() ; J++)
	{ pttoMaterial = (*theMaterialTable)[J];     
		if(pttoMaterial->GetName() == materialChoice)
		{
			AbsorberMaterial = pttoMaterial;
			logicAbsorber->SetMaterial(pttoMaterial); 
		}             
	}
}

void UCNDetectorConstruction::SetAbsorberThickness(G4double val)
{
	// change Absorber thickness and recompute the calorimeter parameters
	AbsorberThickness = val;
	ComputeCalorParameters();
}  

/////////////////////////////////////////////////////////////////////////////
//
//

void UCNDetectorConstruction::SetAbsorberRadius(G4double val)
{
	// change the transverse size and recompute the calorimeter parameters
	AbsorberRadius = val;
	ComputeCalorParameters();
}  

////////////////////////////////////////////////////////////////////////////
//
//

void UCNDetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name   
  G4Material* pttoMaterial;
  for (size_t J=0 ; J<theMaterialTable->size() ; J++)
   { pttoMaterial = (*theMaterialTable)[J];     
     if(pttoMaterial->GetName() == materialChoice)
        {
	  WorldMaterial = pttoMaterial;
          logicWorld->SetMaterial(pttoMaterial); 
        }             
   }
}

////////////////////////////////////////////////////////////////////////////
//
//

void UCNDetectorConstruction::SetWorldSizeZ(G4double val)
{
  worldchanged=true;
  WorldSizeZ = val;
  ComputeCalorParameters();
}  

///////////////////////////////////////////////////////////////////////////
//
//

void UCNDetectorConstruction::SetWorldSizeR(G4double val)
{
  worldchanged=true;
  WorldSizeR = val;
  ComputeCalorParameters();
}  

//////////////////////////////////////////////////////////////////////////////
//
//

void UCNDetectorConstruction::SetAbsorberZpos(G4double val)
{
	zAbsorber  = val;
	ComputeCalorParameters();
}  

////////////////////////////////////////////////////////////////////////////
//
//
  
void UCNDetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//
//
////////////////////////////////////////////////////////////////////////////
