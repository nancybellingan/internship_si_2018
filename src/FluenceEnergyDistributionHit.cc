// // NOT USED BY  MY CODE, TO FULLY REMOVE
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
// $Id: FluenceEnergyDistributionHit.cc 69706 2013-05-13 09:12:40Z gcosmo $
//
// \file FluenceEnergyDistributionHit.cc
// \brief Implementation of the FluenceEnergyDistributionHit class

/*#include "FluenceEnergyDistributionHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<FluenceEnergyDistributionHit>* FluenceEnergyDistributionHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Distribution of the events constructor
FluenceEnergyDistributionHit::FluenceEnergyDistributionHit()
 : G4VHit(),
   fTrackID(-1),
   fCellFlux(0.),
   fKineticEnergy(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Destructor
FluenceEnergyDistributionHit::~FluenceEnergyDistributionHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// description of G4Hit. the name of the fluence hit distribution is "right"
FluenceEnergyDistributionHit::FluenceEnergyDistributionHit(const FluenceEnergyDistributionHit& right)
  : G4VHit()
{
  fTrackID   = right.fTrackID;
  fCellFlux = right.fCellFlux;
  fKineticEnergy = right.fKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// operator that copy and return the memory address of FluenceEnergyDistributionHit
const FluenceEnergyDistributionHit& FluenceEnergyDistributionHit::operator=(const FluenceEnergyDistributionHit& right)
{
  fTrackID   = right.fTrackID;
  fCellFlux = right.fCellFlux;
  fKineticEnergy = right.fKineticEnergy;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// it compares right and FluenceEnergyDistributionHit and returns "right" memory address
G4int FluenceEnergyDistributionHit::operator==(const FluenceEnergyDistributionHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FluenceEnergyDistributionHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {*/
 /*   G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
	pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FluenceEnergyDistributionHit::Print()
{
  G4cout
     << "  trackID: " << fTrackID << " cellFlux: " << fCellFlux*(1/(cm*cm))
     << " kinEnergy: "
     << std::setw(7) << G4BestUnit(fKineticEnergy,"Energy")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
*/
