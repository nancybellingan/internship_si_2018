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
/// \file runAndEvent/RE02/include/PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.hh 66501 2012-12-19 09:25:23Z gcosmo $
//
 
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include <fstream>
#include <cmath>
#include "PiotPhaseSpace.hh"

class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;
class G4GeneralParticleSource;

//
/// User primary particle generator class
///
/// - void GeneratePrimaries(G4Event*)
///     an incident particle is proton with 150 MeV energy at the position 
///     (x,y,-100 cm) toward the (0,0,1) direction. The x and y positions are
///     uniformly varied from -5 mm to 5 mm, respectively.
//
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();    
   ~PrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

/*	void setGunEnergy(const G4double &gEnergy);
	G4double getGunEnergy() { return fGunEnergy; };

	void updateSigmaGunEnergy()
			{ fSigmaGunEnergy = (fPercentFWHM*fGunEnergy)/(2*sqrt(2*log(2)));};
	// For FWHM = 2.5% percentFWHM = 0.025
	void setFWHMGunEnergy(const G4double& percentFWHM)
					       { fPercentFWHM = percentFWHM;
				             updateSigmaGunEnergy(); };

	std::ifstream fReadPhsFile; */

  private:
/*	G4double fSigmaPosition; // Initial beam spot size in x-y plane.
    G4ParticleGun* fParticleGun; //the particle of the primaries
        G4double fNumOfParticle; //number of particles per event
        G4String fParticleName; //name of the particle
        G4double fGunEnergy;
        G4double fSigmaGunEnergy; //initial gun energy deviation standard
        G4double fPercentFWHM; //Full Width Half Maximum scaled on the gun energy
	G4ParticleTable *fParticleTable;
	G4ParticleDefinition *fParticleDefinition;
        std::vector<double> fProbVec; // is this the vector of probabilities corresponding to the energy (see fEnergyVec)?
        std::vector<double> fEnergyVec;
        std::vector<double> fPDFVec; //need explanation, maybe probability density function
	PrimaryGeneratorMessenger *fMessenger;
	DetectorConstruction  *fGeometry;
	size_t fNotInList;

        G4double gaussianShoot(const G4double &meanValue, //describe the gaussian distribution given a mean and a sigma
                                               const G4double &sigmaValue); */
PrimaryGeneratorMessenger *fMessenger;
           G4GeneralParticleSource* fParticleGun;
};

//

#endif


