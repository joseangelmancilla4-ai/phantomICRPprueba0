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
// TETActionInitialization.cc
// 
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Adaptado por José & GPT para uso conjunto con DOI-PET
//


//TETActionInitialization.cc
#include "TETActionInitialization.hh"

#include "doiPETRunAction.hh"
#include "doiPETEventAction.hh"
#include "doiPETPrimaryGeneratorAction.hh"
#include "doiPETSteppingAction.hh"
#include "TETRunAction.hh"
#include "TETSteppingAction.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

TETActionInitialization::TETActionInitialization(TETModelImport* _tetData, G4String _output, doiPETAnalysis* analysisMan)
 : G4VUserActionInitialization(),
   fTetData(_tetData),
   fOutput(_output),
   analysis(analysisMan)
{ }

void TETActionInitialization::BuildForMaster() const
{
    G4cout << "[DEBUG] BuildForMaster() ejecutándose en hilo maestro." << G4endl;
    SetUserAction(new TETRunAction(fTetData, fOutput));
}

void TETActionInitialization::Build() const
{
    G4cout << "[DEBUG] >>> Entrando a TETActionInitialization::Build()" << G4endl;

    auto* primary = new doiPETPrimaryGeneratorAction(fTetData);
    SetUserAction(primary);
    G4cout << "[DEBUG] >>> doiPETPrimaryGeneratorAction creado e instalado como generador primario." << G4endl;

    SetUserAction(new TETRunAction(fTetData, fOutput));
    G4cout << "[DEBUG] >>> TETRunAction listo." << G4endl;

    SetUserAction(new doiPETRunAction());
    G4cout << "[DEBUG] >>> doiPETRunAction listo." << G4endl;

    SetUserAction(new TETSteppingAction());
    SetUserAction(new doiPETSteppingAction());
    G4cout << "[DEBUG] >>> Stepping actions listas." << G4endl;
}



