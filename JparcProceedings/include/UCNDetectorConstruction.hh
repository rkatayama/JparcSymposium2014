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
// UCN detector construction, 9.9.04, peter fierlinger

#ifndef UCNDetectorConstruction_h
#define UCNDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class UCNDetectorMessenger;
class UCNFieldSetup;


class UCNDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    UCNDetectorConstruction();
   ~UCNDetectorConstruction();

  public:
     
     void SetAbsorberMaterial (G4String);     
     void SetAbsorberThickness(G4double);     
     void SetAbsorberRadius(G4double);          
      
     void SetAbsorberZpos(G4double);
     
     void SetWorldMaterial(G4String);
     void SetWorldSizeZ(G4double);
     void SetWorldSizeR(G4double);
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     void PrintCalorParameters(); 
                    
     G4Material* GetWorldMaterial()    {return WorldMaterial;};
     G4double GetWorldSizeZ()          {return WorldSizeZ;}; 
     G4double GetWorldSizeR()          {return WorldSizeR;};
     
     
     G4double GetAbsorberZpos()        {return zAbsorber;}; 
     G4double GetzstartAbs()           {return zstartAbs;};
     G4double GetzendAbs()             {return zendAbs;};

     G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};
     G4double    GetAbsorberThickness() {return AbsorberThickness;};      
     G4double GetAbsorberRadius()       {return AbsorberRadius;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetAbsorber()   {return physiAbsorber;};
     G4LogicalVolume* GetLogicalAbsorber()    {return logicAbsorber;};

		 //add by katayama
		 UCNFieldSetup* GetFieldSetup() {return fEmFieldSetup;}

  private:
     
  //G4Tubs*            solidWorld;    // pointer to the solid World 
     G4Box*            solidWorld;    // pointer to the solid World 
     G4LogicalVolume*   logicWorld;    // pointer to the logical World
     G4VPhysicalVolume* physiWorld;    // pointer to the physical World

     G4Tubs*            solidAbsorber;  // pointer to the solid Absorber
     G4LogicalVolume*   logicAbsorber;  // pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber;  // pointer to the physical Absorber



     UCNFieldSetup* fEmFieldSetup;     // pointer to the field helper
     UCNDetectorMessenger* detectorMessenger;  // pointer to the Messenger
     
  //G4Material*        FoilMaterial;
  //G4Material*        GuideMaterial;
     G4Material*        sD2Material;
     G4Material*        AluMaterial;
     G4Material*        IronMaterial;
     G4Material*        Nickel58Material;
     G4Material*        NickelMaterial;
     G4Material*        Shutter_sMaterial;
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;
     G4double           AbsorberRadius;
     G4bool             worldchanged;
     G4Material*        Vacuum;

     G4Material*        WorldMaterial;
     G4double           WorldSizeR;
     G4double           WorldSizeZ;

     G4double           zAbsorber;
     G4double           zstartAbs , zendAbs ;

     G4double           Shutter1x,Shutter1y,Shutter1z,Shutter1Radius ;
     G4double           Shutter1Thickness;
     G4Material*        Shutter1Material;

  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void UCNDetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter

     zstartAbs = zAbsorber-0.5*AbsorberThickness; 
     zendAbs   = zAbsorber+0.5*AbsorberThickness; 

}

#endif
