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
// TETPrimaryGeneratorAction.hh
// 
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//

//TETPrimaryGeneratorAction.hh
#ifndef doiPETPrimaryGeneratorAction_h
#define doiPETPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "globals.hh"
#include "TETModelImport.hh"
#include <vector>

class G4Event;
class G4GeneralParticleSource;

/// Generador primario adaptado para fuentes volumétricas en órganos específicos
class doiPETPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    // Constructor modificado para recibir el puntero del modelo ICRP
    explicit doiPETPrimaryGeneratorAction(TETModelImport* tetData = nullptr);
    virtual ~doiPETPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    const G4GeneralParticleSource* GetParticleGun() const { return fParticleGun; }
    void ActivityValue();

private:
    G4GeneralParticleSource* fParticleGun;
    TETModelImport* fTetData;  // acceso al modelo voxelizado del fantoma

    // IDs de materiales/órganos relevantes del tórax
    std::vector<G4int> organIDs = {
        8700,  // Heart_wall
        8800,  // Blood_in_heart_chamber
        9700,  // Lung_left
        9900,  // Lung_right
        9500   // Liver
    };
};

#endif


