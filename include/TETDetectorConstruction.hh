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
// TETDetectorConstruction.hh
// 
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//

//TETDetectorConstruction.hh
#ifndef TETDetectorConstruction_h
#define TETDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include <cmath>

#include "globals.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "TETPSEnergyDeposit.hh"

#include "G4SystemOfUnits.hh"

#include "TETModelImport.hh"
#include "TETParameterisation.hh"
#include "doiPETGlobalParameters.hh"
#include "doiPETAnalysis.hh"


// *********************************************************************
// This is UserDetectorConstruction class that defines the geometry
// -- Construct: construct Geometry by three methods listed below.
//  └-- SetupWorldGeometry: Defines the world box (10*10*10 m3) and,
//                          phantom container which has 10 cm-margins from
//                          the bounding box of phantom
//  └-- ConstructPhantom: Define the phantom geometry by using
//                        G4PVParameterised class
//  └-- PrintPhantomInformation: Print overall phantom information
//
// -- ConstructSDandField: Setup the MultiFunctionalDetector with energy
//                         deposition scorer, and attach it to phantom
//                         geometry
// *********************************************************************

class doiPETAnalysis;

class TETDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	TETDetectorConstruction(TETModelImport* tetData);
	virtual ~TETDetectorConstruction();

	virtual G4VPhysicalVolume* Construct() override;
	virtual void ConstructSDandField() override;

private:
	void SetupWorldGeometry();
	void ConstructPhantom();
	void PrintPhantomInformation();

	void DefineMaterials();

	G4VPhysicalVolume* fWorldPhysical;
	G4LogicalVolume*   fContainer_logic;

	TETModelImport*    fTetData;
	G4ThreeVector      fPhantomSize;
	G4ThreeVector      fPhantomBoxMin, fPhantomBoxMax;
	G4int              fNOfTetrahedrons;

	G4LogicalVolume*   fTetLogic;

	//materials
    G4Material* air;
    G4Material* pmma;
    G4Material* water;
    G4Material* polyethylene;
    G4Material* polyethylene_NEMA;
    //G4Material* inflatedLung;
    G4Material* polystyrene;
    G4Material* Aluminum;

    //elements for GSO
    G4Element*  O;
    G4Element* Si;
    G4Element* Gd;
    G4Material* GSO;

    G4Material* crystalMaterial;

	G4bool  fCheckOverlaps;
    G4bool isotopes;

	//size of world
	G4double worldSizeX;
	G4double worldSizeY;
	G4double worldSizeZ;

	G4int blockIndex;
	//G4int AlCase_Index;
	G4int crystalIndex;

	G4LogicalVolume* world_logicalV;
    G4VPhysicalVolume* world_physicalV;

	G4double sizeOfAirBox_DOI;
    G4double sizeOfAirBox_axial;
    G4double sizeOfAirBox_tangential;

    G4double sizeOfBlockDetector_DOI;
    G4double sizeOfBlockDetector_axial;
    G4double sizeOfBlockDetector_tangential;

	doiPETAnalysis* pAnalysis;

	//detector block
	G4LogicalVolume* blockDetector_logicalV;
	G4VPhysicalVolume* blockDetector_physicalV;

	//air volume to fill the detector block
	G4LogicalVolume* airBox_logicalV;
	G4VPhysicalVolume* airBox_physicalV;

	//crystals
	G4LogicalVolume* crystal_logicalV;
	G4VPhysicalVolume* crystal_physicalV;

	//detector position
	G4double detectorPositionX;
	G4double detectorPositionY;
	G4double detectorPositionZ;

	//crystal position
	G4double crystalPositionX;
	G4double crystalPositionY;
	G4double crystalPositionZ;

    G4double thetaDetector; //The azimuthal angle for arranging the detector in the PET ring

};
#endif
