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
// $Id: FluenceEnergyDistributionHit.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file FluenceEnergyDistributionHit.hh
/// \brief Definition of the FluenceEnergyDistributionHit class

#ifndef FluenceEnergyDistributionHit_h
#define FluenceEnergyDistributionHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"
#include <map>

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class FluenceEnergyDistributionHit : public G4VHit
{
  public:
    FluenceEnergyDistributionHit();
    FluenceEnergyDistributionHit(const FluenceEnergyDistributionHit&);
    virtual ~FluenceEnergyDistributionHit();

    // operators
    const FluenceEnergyDistributionHit& operator=(const FluenceEnergyDistributionHit&);
    G4int operator==(const FluenceEnergyDistributionHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID(const G4int& trID) { fTrackID = trID; };
    void SetCellFlux(const G4double& CF) { fCellFlux = CF; };
		void SetKineticEnergy(const G4double& kE) { fKineticEnergy = kE; };

    // Get methods
    G4double GetTrackID() const     { return fTrackID; };
    G4double GetCellFlux() const     { return fCellFlux; };
		G4double GetKineticEnergy() const { return fKineticEnergy; };

  private:

			G4int fTrackID;
      G4double fCellFlux;
			G4double fKineticEnergy;	
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<FluenceEnergyDistributionHit> FluenceEnergyDistributionHitsCollection;

extern G4ThreadLocal G4Allocator<FluenceEnergyDistributionHit>* FluenceEnergyDistributionHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* FluenceEnergyDistributionHit::operator new(size_t)
{
  if(!FluenceEnergyDistributionHitAllocator)
      FluenceEnergyDistributionHitAllocator = new G4Allocator<FluenceEnergyDistributionHit>;
  return (void *) FluenceEnergyDistributionHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void FluenceEnergyDistributionHit::operator delete(void *hit)
{
  FluenceEnergyDistributionHitAllocator->FreeSingle((FluenceEnergyDistributionHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
