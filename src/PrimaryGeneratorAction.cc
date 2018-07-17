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
/// \file runAndEvent/RE02/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.cc 68026 2013-03-13 13:45:22Z gcosmo $
//

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"    
#include <cmath>
#include "config.h"
//============================================================================

PrimaryGeneratorAction::PrimaryGeneratorAction()
 :	G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{

//	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
//	G4ParticleDefinition* particle = particleTable->FindParticle(fParticleName);

        // uncomment for particlegun monoergetic beam
        // fParticleGun = new G4ParticleGun(fNumOfParticle);
       // uncomment for GPS general particle source for am-be
        fParticleGun = new G4GeneralParticleSource();
       // fParticleGun->SetParticleDefinition("neutron");
      // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.));
        // fParticleGun->SetParticleEnergy(fGunEnergy);
        // updateSigmaGunEnergy();

	fMessenger = new PrimaryGeneratorMessenger(this);

//----------------------------------------------------------------------------
	// Checking Gaussian Distribution

	/*
	std::ofstream outw;
	outw.open("distri.out");
	outw << "Sigma Gun: " << fSigmaGunEnergy << std::endl;
	for(size_t i = 0; i<1000000;i++)
	{
		G4double distri = gaussianShoot(fGunEnergy,fSigmaGunEnergy);
		outw << distri << std::endl;
	}
	outw.close();
	*/

//----------------------------------------------------------------------------

	//G4double position = -157.5*cm+75.*cm;

// Initial beam spot size in sigma.; This is not a part of ParticleGun.
//	fSigmaPosition = 10.* mm;
        //fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm, 0.*cm, position));

}

//============================================================================

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fMessenger;
}

//======generate primary particles and relative momentum======================================================================
#include <mutex>
static std::mutex mutexFileWrite;
static int counterForFlileFlush = 0;
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

	counterForFlileFlush++;
	// G4ThreeVector position = fParticleGun->GetParticlePosition();
	fParticleGun->GeneratePrimaryVertex(anEvent);
	G4double energy = fParticleGun->GetParticleEnergy();

	std::lock_guard<std::mutex> lock(mutexFileWrite);
	*(conf()->edistr) << energy << '\n';
	//G4cout << "energy: " << energy << G4endl;

	if(counterForFlileFlush > 128){
		conf()->edistr->flush();
		counterForFlileFlush = 0;
	}
}

//============================================================================

/* void PrimaryGeneratorAction::setGunEnergy(const G4double &gEnergy)
{
	fGunEnergy=gEnergy;
	updateSigmaGunEnergy();
	fParticleGun->SetParticleEnergy(fGunEnergy);
}

//============================================================================

G4double PrimaryGeneratorAction::gaussianShoot(
	const G4double &meanValue, const G4double &sigmaValue)
{
	if(sigmaValue > 0)
	{
		return G4RandGauss::shoot(meanValue,sigmaValue);
	}
	else 
	{
		return meanValue;
	}
} */

