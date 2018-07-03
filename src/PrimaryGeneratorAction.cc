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
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "config.h"
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
#include <time.h>
//============================================================================

PrimaryGeneratorAction::PrimaryGeneratorAction()
    :	G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{

	ambedistribution.resize(conf()->ebin.size());
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

	out.open("EDistribution.dat");


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
	/*
#ifdef FLUKA_PSF
	std::string ifile="./../../../RPTC/ginputs/11x11/140MeV";
	fReadPhsFile.open(ifile.c_str());
	if(!fReadPhsFile)
	{
		std::cout << "Cannot find or open file: "
				  << ifile << std::endl;
		exit(-1);
	}
#endif

#ifdef SPECSCHUSS
	fGeometry = (DetectorConstruction*)
			(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	fEnergyVec = fGeometry->GetEnergyVec();
	fPDFVec 	 = fGeometry->GetPDFVec();
#endif
*/
}

//============================================================================

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fMessenger;
#ifdef FLUKA_PSF
	fReadPhsFile.close();
	G4cout << "WARNING: Not in Particle List: " << fNotInList << G4endl;
#endif
	out.close();
}

//======generate primary particles and relative momentum======================================================================



void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
	static int number = 0;
	// G4ThreeVector position = fParticleGun->GetParticlePosition();
	number++;
	timespec tv;
	clock_gettime (CLOCK_MONOTONIC, &tv);
	G4double energy = fParticleGun->GetParticleEnergy(); //TODO: Check Function call
	out << energy << "," << number << "," << tv.tv_sec << "." << tv.tv_nsec << '\n';
	out.flush();//this will force to write on disk immediately

	for (uint i=0;i<conf()->ebin.size();i++){
		if (i==0){
			if(energy>0 && energy <conf()->ebin[i]){
				ambedistribution[i]++;
			}
		}else {
			if(energy>=conf()->ebin[i-1] && energy<conf()->ebin[i]){
				ambedistribution[i]++;
			}
		}

		fParticleGun->GeneratePrimaryVertex(anEvent);

	}
}

std::vector<G4double> PrimaryGeneratorAction::getambedistribution() {
	return ambedistribution;
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

