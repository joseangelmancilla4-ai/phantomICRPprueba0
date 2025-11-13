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


// TETPrimaryGeneratorAction.cc
#include "doiPETPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include <iostream>

doiPETPrimaryGeneratorAction::doiPETPrimaryGeneratorAction(TETModelImport* tetData)
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(nullptr),
      fMessenger(nullptr),
      fTetData(tetData)
{
    fParticleGun = new G4GeneralParticleSource();
    fParticleGun->SetParticleByName("e+");
    fParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(511 * keV);

    DefineCommands();

    G4cout << "[doiPETPrimaryGeneratorAction] Fuente PET inicializada (511 keV e+)." << G4endl;
}

doiPETPrimaryGeneratorAction::~doiPETPrimaryGeneratorAction()
{
    delete fParticleGun;
    delete fMessenger;
}

void doiPETPrimaryGeneratorAction::DefineCommands()
{
    // Prefijo único, evita conflicto con /gps/
    fMessenger = new G4GenericMessenger(this, "/organSource/", "Control de fuentes por órgano");

    fMessenger->DeclareMethod("addOrganSource", &doiPETPrimaryGeneratorAction::AddOrganSource,
                              "Agrega una fuente asociada a un órgano (uso: /organSource/addOrganSource organID intensidad)")
        .SetParameterName("organID", "intensity", true);
}

void doiPETPrimaryGeneratorAction::AddOrganSource(G4int organID, G4double intensity)
{
    fOrganIntensities[organID] = intensity;
    G4cout << "→ Fuente agregada para órgano ID " << organID
           << " con intensidad relativa " << intensity << G4endl;
}

void doiPETPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if (!fTetData || fOrganIntensities.empty()) {
        fParticleGun->GeneratePrimaryVertex(anEvent);
        return;
    }

    // Elegir órgano proporcionalmente a su intensidad
    G4double total = 0.;
    for (auto& kv : fOrganIntensities) total += kv.second;

    G4double r = G4RandFlat::shoot(0., total);
    G4int chosenID = -1, i = 0;
    G4double accum = 0.;
    for (auto& kv : fOrganIntensities) {
        accum += kv.second;
        if (r <= accum) { chosenID = kv.first; break; }
    }

    if (chosenID < 0) chosenID = fOrganIntensities.begin()->first;

    // Buscar un voxel dentro del órgano
    G4ThreeVector pos;
    G4bool found = false;
    for (int n = 0; n < 50000; n++) {
        G4double x = G4RandFlat::shoot(fTetData->GetXMin(), fTetData->GetXMax());
        G4double y = G4RandFlat::shoot(fTetData->GetYMin(), fTetData->GetYMax());
        G4double z = G4RandFlat::shoot(fTetData->GetZMin(), fTetData->GetZMax());
        if (fTetData->GetMaterialIndex(x, y, z) == chosenID) {
            pos = {x, y, z};
            found = true;
            break;
        }
    }

    if (!found) {
        G4cout << "⚠️ No se encontró voxel en órgano ID " << chosenID << G4endl;
        fParticleGun->GeneratePrimaryVertex(anEvent);
        return;
    }

    fParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(pos);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

