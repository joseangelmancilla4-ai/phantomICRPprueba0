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

//GEANT4 - Depth-of-Interaction enabled Positron emission tomography (PET) advanced example 

//Authors and contributors

// Author list to be updated, with names of co-authors and contributors from National Institute of Radiological Sciences (NIRS)

// Abdella M. Ahmed (1, 2), Andrew Chacon (1, 2), Harley Rutherford (1, 2),
// Hideaki Tashima (3), Go Akamatsu (3), Akram Mohammadi (3), Eiji Yoshida (3), Taiga Yamaya (3)
// Susanna Guatelli (2), and Mitra Safavi-Naeini (1, 2)

// (1) Australian Nuclear Science and Technology Organisation, Australia
// (2) University of Wollongong, Australia
// (3) National Institute of Radiological Sciences, Japan

//doiPETPrimaryGeneratorAction.cc
#ifndef doiPETPrimaryGeneratorAction_h
#define doiPETPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "TETModelImport.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include <map>
#include <string>

class G4Event;

class doiPETPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction, public G4UImessenger
{
public:
    doiPETPrimaryGeneratorAction(TETModelImport* tetData = nullptr);
    virtual ~doiPETPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);
    virtual void SetNewValue(G4UIcommand*, G4String);

    void AddOrganSource(G4int organID, const G4String& isotope, G4double intensity, const G4String& mode);

private:
    bool GetIsotopeInfo(const std::string& iso, G4int &Z, G4int &A, G4double &Emax_keV) const;

    struct SourceSpec {
        std::string isotope;
        double intensity;
        std::string mode;
    };

    std::map<G4int, SourceSpec> fOrganSources;

    G4GeneralParticleSource* fParticleGun;
    TETModelImport* fTetData;

    // UI
    G4UIdirectory* fDirectory;
    G4UIcommand* fAddSourceCmd;
};

#endif

