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
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//

//TETDetectorConstruction.cc
#include "TETDetectorConstruction.hh"

#include "G4VisAttributes.hh"

#include "doiPETAnalysis.hh"

TETDetectorConstruction::TETDetectorConstruction(TETModelImport* _tetData)
:fWorldPhysical(nullptr), fContainer_logic(nullptr), fTetData(_tetData), fTetLogic(nullptr)
{
 // initialisation of the variables for phantom information
 fPhantomSize     = fTetData -> GetPhantomSize();
 fPhantomBoxMin   = fTetData -> GetPhantomBoxMin();
 fPhantomBoxMax   = fTetData -> GetPhantomBoxMax();
 fNOfTetrahedrons = fTetData -> GetNumTetrahedron();
}

TETDetectorConstruction::~TETDetectorConstruction()
{
  delete fTetData;
}

G4VPhysicalVolume* TETDetectorConstruction::Construct()
{
 SetupWorldGeometry();
 ConstructPhantom();
 PrintPhantomInformation();
 //DefineMaterials();
 return fWorldPhysical;
}

void TETDetectorConstruction::SetupWorldGeometry()
{
 // Define the world box (size: 10*10*10 m3)
 //
 G4double worldXYZ = 10. * m;
 G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

 G4VSolid* worldSolid = new G4Box("worldSolid", worldXYZ/2, worldXYZ/2, worldXYZ/2);

 auto* worldLogical = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");

 fWorldPhysical = new G4PVPlacement(nullptr,G4ThreeVector(), worldLogical,"worldPhysical", nullptr, false,0,false);

 // Define the phantom container (10-cm margins from the bounding box of phantom)
 //
 auto* containerSolid = new G4Box("phantomBox", fPhantomSize.x()/2 + 10.*cm,
					           fPhantomSize.y()/2 + 10.*cm,
						   fPhantomSize.z()/2 + 10.*cm);

 fContainer_logic = new G4LogicalVolume(containerSolid, vacuum, "phantomLogical");

 new G4PVPlacement(nullptr, G4ThreeVector(), fContainer_logic, "PhantomPhysical",
	           worldLogical, false, 0);

 fContainer_logic->SetOptimisation(TRUE);
 fContainer_logic->SetSmartless( 0.5 ); // for optimization (default=2)









  /* ----------------------------------------------
     ----------------------------------------------
     ----------------------------------------------
     ----------------------------------------------
     ----------------------------------------------
              DEFINE PET DefineMaterials
     ----------------------------------------------
     ----------------------------------------------
     ----------------------------------------------
     ----------------------------------------------
     ---------------------------------------------- */


  G4NistManager* nist = G4NistManager::Instance();

  //Define air
  air = nist->FindOrBuildMaterial("G4_AIR");

  //Define PMMA
  pmma  = nist->FindOrBuildMaterial("G4_PLEXIGLASS"); //Default 1.19 g/cm3

  //Define water
  water  = nist->FindOrBuildMaterial("G4_WATER");

  //Defining polyethylene from NIST and modifying the density
  polyethylene = nist->BuildMaterialWithNewDensity("polyethylene","G4_POLYETHYLENE",0.959*g/cm3);
  polyethylene->GetIonisation()->SetMeanExcitationEnergy(56*eV);

  //polyethylene_NEMA, defualt density 0.94 g/cm3, and excitation energy = 57.4 eV
  polyethylene_NEMA = nist->FindOrBuildMaterial("G4_POLYETHYLENE");


  //Define expanded polystyrene by modifiying the density to mimic lung phantom used in phantom experiment
  polystyrene = nist->BuildMaterialWithNewDensity( "polystyrene","G4_POLYSTYRENE",0.3*g/cm3);

  isotopes = false;

  //Defile Aluminum material for the detetor cover
  Aluminum = nist->FindOrBuildMaterial("G4_Al", isotopes);

  //Define elements for the GSO  crystal (scintillator)
  O = nist->FindOrBuildElement("O" , isotopes);
  Si = nist->FindOrBuildElement("Si", isotopes);
  Gd = nist->FindOrBuildElement("Gd", isotopes);


  //define GSO crystal for PET detector
  GSO = new G4Material("GSO", 6.7*g/cm3, 3);
  GSO->AddElement(Gd, 2);
  GSO->AddElement(Si, 1);
  GSO->AddElement(O,  5);
  crystalMaterial   = nist->FindOrBuildMaterial("GSO");


  /* ----------------------------------------------
     ----------------------------------------------
     ----------------------------------------------
     ----------------------------------------------
					CONSTRUCT PET
     ----------------------------------------------
     ----------------------------------------------
     ----------------------------------------------
     ---------------------------------------------- */

     //auto* worldLogical = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");

     //WorldPhysical = new G4PVPlacement(nullptr,G4ThreeVector(), worldLogical,"worldPhysical", nullptr, false,0,false);



	//size of the world
	//worldSizeX = 2 * m;//0.5 *mm
	//worldSizeY = 2 * m;//0.5*mm
	//worldSizeZ = 4. * m;

	// Define a solid shape to describe the world volume
	//G4Box* solid_world = new G4Box("world", worldSizeX/2., worldSizeY/2., worldSizeZ/2.);

	// Define a logical volume for the world volume
	//world_logicalV = new G4LogicalVolume(solid_world, air, "world_logicalV", 0,0,0);
	//worldLogical

	// Define the physical world volume

	world_physicalV = new G4PVPlacement(0,G4ThreeVector(),worldLogical, "world_physicalV", 0, false, 0, fCheckOverlaps);
	worldLogical->SetVisAttributes (G4VisAttributes::GetInvisible());

	//NOTE!!!
	//The scanner specification (like size and number of crystals in each detctor) are given in the "doiPETGlobalParameters.hh" header file.

	//Each block detector is identified with its unique number, blockIndex.
	blockIndex = 0;

	//Each crystal is identified with its unique number, crystalIndex
	crystalIndex = 0;




	//Define air volume (box) to fill the detector block. Crystal elements (scintillators) is then placed.
	sizeOfAirBox_DOI = (numberOfCrystal_DOI * sizeOfCrystal_DOI) + (numberOfCrystal_DOI - 1)*crystalGap_DOI;
	sizeOfAirBox_axial = (numberOfCrystal_axial * sizeOfCrystal_axial) + (numberOfCrystal_axial - 1)*crystalGap_axial;
	sizeOfAirBox_tangential = (numberOfCrystal_tangential * sizeOfCrystal_tangential) + (numberOfCrystal_tangential - 1)*crystalGap_tangential;





	//This will pass the size of the detector block so that the position of the PMT will be calculated with respect to the axis of the detector block.
	//see doiPETAnalysis.cc


	//Initialize analysis
	//doiPETAnalysis* ptrAnalysis = doiPETAnalysis::GetInstance();

	doiPETAnalysis* pAnalysis = doiPETAnalysis::GetInstance();
	pAnalysis->GetSizeOfDetector(sizeOfAirBox_DOI,sizeOfAirBox_tangential, sizeOfAirBox_axial);






	G4cout<<"size of crytal element: "<<sizeOfCrystal_tangential<<" "<<sizeOfCrystal_axial<<" "<<sizeOfCrystal_DOI<<G4endl;
	G4cout<<"Size of detector block (without Al cover): "<<sizeOfAirBox_tangential<<" "<<sizeOfAirBox_axial<<" "<<sizeOfAirBox_DOI<<G4endl;



	//Define the size of the detector block.
	sizeOfBlockDetector_DOI = sizeOfAirBox_DOI + AluminumCoverThickness;
	sizeOfBlockDetector_axial = sizeOfAirBox_axial + AluminumCoverThickness;
	sizeOfBlockDetector_tangential = sizeOfAirBox_tangential + AluminumCoverThickness;



	//Define solid shape for the detector block
	G4Box* blockDetector = new G4Box("blockDetector",sizeOfBlockDetector_DOI/2,sizeOfBlockDetector_tangential/2,sizeOfBlockDetector_axial/2);

	//Define the logical volume for the detector block
	blockDetector_logicalV = new G4LogicalVolume(blockDetector,Aluminum,"blockDetector_logicalV", 0,0,0);



	//Define air (box) inside the detector block. Crystal elements will be placed in it.
	G4Box* airBox = new G4Box("airBox", sizeOfAirBox_DOI/2, sizeOfAirBox_tangential/2,sizeOfAirBox_axial/2);

	//Define the logical volume
	airBox_logicalV = new G4LogicalVolume(airBox,air,"airBox_logicalV", 0,0,0);

	//Define its physical volume and place it inside the detector block
	airBox_physicalV = new G4PVPlacement (0,G4ThreeVector(0,0,0),airBox_logicalV,"airBox_physicalV", blockDetector_logicalV,false,0,fCheckOverlaps);




		///////////////////////////////////////// Arrange the PET ring and place the PET detectors in the ring(s) ////////////////////////////////////

		for(G4int Ring = 0; Ring< numberOfRings; Ring++)
		{
			//place the detectors in a ring along the axial direction. Note that the ring gap between two adjcent rings is measured from scintillator to scintillator.  It does not include the Aluminum thickness cover.


// Añade el desplazamiento (offset) al final
G4double Z_OFFSET = 20*cm; // Define tu desplazamiento aquí, por ejemplo 100 mm

                        detectorPositionZ = (Ring-((G4double)numberOfRings)/2 + 9.5)*(sizeOfBlockDetector_axial + ringGap - AluminumCoverThickness) + Z_OFFSET;

			//detectorPositionZ = (Ring-((G4double)numberOfRings)/2 + 9.5)*(sizeOfBlockDetector_axial + ringGap - AluminumCoverThickness);

			for(G4int i = 0; i<numberOfDetector_perRing; i++)
			{
			//The azimuthal angle to arrange the detectors in a ring
			thetaDetector = (double)(i*twopi/numberOfDetector_perRing);

			//The radius of the scanner is measured from opposing crystal (scintillator) faces. It does not include the Aluminum thickness cover.
			detectorPositionX = (scannerRadius + sizeOfBlockDetector_DOI/2 - AluminumCoverThickness/2)*std::cos(thetaDetector);
			detectorPositionY = (scannerRadius + sizeOfBlockDetector_DOI/2 - AluminumCoverThickness/2)*std::sin(thetaDetector);

			//Define the rotation matrix for correct placement of detetors
			G4RotationMatrix rotm_PET = G4RotationMatrix();
			rotm_PET.rotateZ(thetaDetector);
			G4ThreeVector uz_PET = G4ThreeVector(detectorPositionX,detectorPositionY,detectorPositionZ);
			G4Transform3D transform = G4Transform3D(rotm_PET,uz_PET);

			//Define the physical volume of the detectors.
			blockDetector_physicalV = new G4PVPlacement (transform,blockDetector_logicalV,"blockDetector_physicalV", worldLogical,false,blockIndex,fCheckOverlaps);
			blockIndex++;
			//G4cout<<Ring<<" "<<detectorPositionX- ((sizeOfBlockDetector_DOI - AluminumCoverThickness)/2)*cos(thetaDetector)<<" "<<detectorPositionY- ((sizeOfBlockDetector_DOI- AluminumCoverThickness)/2)*sin(thetaDetector)<<" "<<detectorPositionZ<<G4endl;

			}
		}

		//Define the solid crystal
		G4VSolid* CrystalSolid = new G4Box("Crystal", sizeOfCrystal_DOI/2., sizeOfCrystal_tangential/2., sizeOfCrystal_axial/2.);

		//Define the local volume of the crystal
		crystal_logicalV = new G4LogicalVolume(CrystalSolid,crystalMaterial,"Crystal_logicalV", 0,0,0);

		//Place the crystals inside the detectors and give them a unique number with crystalIndex.
		for(G4int i_DOI = 0; i_DOI<numberOfCrystal_DOI; i_DOI++){
			crystalPositionX=(i_DOI-((G4double)numberOfCrystal_DOI)/2 + 0.5)*(sizeOfCrystal_DOI + crystalGap_DOI);
			for(G4int i_axial=0; i_axial< numberOfCrystal_axial;i_axial++){
			crystalPositionZ = (i_axial-((G4double)numberOfCrystal_axial)/2 + 0.5)*(sizeOfCrystal_axial + crystalGap_axial);
			for(G4int i_tan=0; i_tan<numberOfCrystal_tangential;i_tan++){
				crystalPositionY=(i_tan-((G4double)numberOfCrystal_tangential)/2 + 0.5)*(sizeOfCrystal_tangential + crystalGap_tangential);

				//G4cout<<crystalIndex<<" "<<crystalPositionX<<" "<<crystalPositionY<<" "<<crystalPositionZ<<G4endl;
				//place the crystal inside the block detector.
				crystal_physicalV = new G4PVPlacement (0, G4ThreeVector (crystalPositionX,crystalPositionY,crystalPositionZ), crystal_logicalV, "Crystal_physicalV", airBox_logicalV,false,crystalIndex/*,fCheckOverlaps*/);
				crystalIndex++;
			}
			}
		}

		//******************  Visualization *****************************//

		//visualization for the block detector
		G4VisAttributes* blockDetectorVisAtt;
		blockDetectorVisAtt = new G4VisAttributes(G4Colour(1,1.0,1.0));
		blockDetectorVisAtt->SetVisibility (true);
		//blockDetectorVisAtt->SetForceWireframe (true);
		blockDetector_logicalV->SetVisAttributes (blockDetectorVisAtt);
		//blockDetector_logicalV->SetVisAttributes (G4VisAttributes::GetInvisible());

		//visualization for the the box filled with air
		G4VisAttributes* airBoxVisAtt;
		airBoxVisAtt = new G4VisAttributes(G4Colour(1,1.0,1.0));
		airBoxVisAtt->SetVisibility (true);
		airBoxVisAtt->SetForceWireframe (true);
		airBox_logicalV->SetVisAttributes (airBoxVisAtt);
		airBox_logicalV->SetVisAttributes (G4VisAttributes::GetInvisible());

		//visualization for the crystal
		G4VisAttributes* crystalVisAtt;
		crystalVisAtt = new G4VisAttributes(G4Colour(0.88,0.55,1.0));
		//crystalVisAtt->SetVisibility (true);
		crystalVisAtt->SetForceWireframe (true);
		crystal_logicalV->SetVisAttributes (crystalVisAtt);
		crystal_logicalV->SetVisAttributes (G4VisAttributes::GetInvisible());


















}

void TETDetectorConstruction::ConstructPhantom()
{
 // Define the tetrahedral mesh phantom as a parameterised geometry
 //
 // solid and logical volume to be used for parameterised geometry
 G4VSolid* tetraSolid = new G4Tet("TetSolid", G4ThreeVector(),
			           G4ThreeVector(1.*cm,0,0),
			           G4ThreeVector(0,1.*cm,0),
			           G4ThreeVector(0,0,1.*cm));

  G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  fTetLogic = new G4LogicalVolume(tetraSolid, vacuum, "TetLogic");

  // physical volume (phantom) constructed as parameterised geometry
  new G4PVParameterised("wholePhantom",fTetLogic,fContainer_logic,
			  kUndefined, fTetData->GetNumTetrahedron(),
			 new TETParameterisation(fTetData));









}

void TETDetectorConstruction::ConstructSDandField()
{
 // Define detector (Phantom SD) and scorer (eDep)
 //
 G4SDManager* pSDman = G4SDManager::GetSDMpointer();
 G4String phantomSDname = "PhantomSD";

 // MultiFunctional detector
 auto* MFDet = new G4MultiFunctionalDetector(phantomSDname);
 pSDman->AddNewDetector( MFDet );

 // scorer for energy depositon in each organ
 MFDet->RegisterPrimitive(new TETPSEnergyDeposit("eDep", fTetData));

 // attach the detector to logical volume for parameterised geometry (phantom geometry)
 SetSensitiveDetector(fTetLogic, MFDet);
}

void TETDetectorConstruction::PrintPhantomInformation()
{
 // print brief information on the imported phantom
 G4cout<< G4endl;
 G4cout.precision(3);
 G4cout<<"   Phantom name               "<<fTetData->GetPhantomName() << " TET phantom"<<G4endl;
 G4cout<<"   Phantom size               "<<fPhantomSize.x()<<" * "<<fPhantomSize.y()<<" * "<<fPhantomSize.z()<<" mm3"<<G4endl;
 G4cout<<"   Phantom box position (min) "<<fPhantomBoxMin.x()<<" mm, "<<fPhantomBoxMin.y()<<" mm, "<<fPhantomBoxMin.z()<<" mm"<<G4endl;
 G4cout<<"   Phantom box position (max) "<<fPhantomBoxMax.x()<<" mm, "<<fPhantomBoxMax.y()<<" mm, "<<fPhantomBoxMax.z()<<" mm"<<G4endl;
 G4cout<<"   Number of tetrahedrons     "<<fNOfTetrahedrons<<G4endl<<G4endl;
}


//void TETDetectorConstruction::DefineMaterials()
//{
//
//}

////////////////////////////////////////////////////////////////////////////////////////


