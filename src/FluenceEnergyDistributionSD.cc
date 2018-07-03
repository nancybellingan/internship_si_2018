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
// $Id: FluenceEnergyDistributionSD.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file FluenceEnergyDistributionSD.cc
/// \brief Implementation of the FluenceEnergyDistributionSD class

#include "FluenceEnergyDistributionSD.hh"
#include "FluenceEnergyDistributionHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VPVParameterisation.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "fstream"
#include "iomanip"
#include <iostream>
#include "config.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Sensitive Detector
FluenceEnergyDistributionSD::FluenceEnergyDistributionSD(
		const G4String& name, const G4int& ni, const G4int& nj, 
		const G4int& nk, const G4int& depth, const G4int& depi, 
		const G4int& depj,const G4int& depk) 
 : G4VSensitiveDetector(name), fDepth(depth), fDepthi(depi), fDepthj(depj), fDepthk(depk), fWeighted(true), fIntEdep(0.)
{
	G4cout << "Creating FluenceEnergyDistributionSD for: " << name << G4endl;
    fNi=ni;
	fNj=nj;
	fNk=nk;

	fEventManager = (G4EventManager*) 
					G4EventManager::GetEventManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FluenceEnergyDistributionSD::~FluenceEnergyDistributionSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FluenceEnergyDistributionSD::Initialize(G4HCofThisEvent*) //Hits Collection of this event
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool FluenceEnergyDistributionSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  


	G4double edep = aStep->GetTotalEnergyDeposit();
	if(edep == 0.) return false; 
	edep*=aStep->GetPreStepPoint()->GetWeight();
	fIntEdep+=edep;
	G4int  index = GetIndex(aStep);
	fLayerEDep[index]+=edep;


	return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void FluenceEnergyDistributionSD::clear()
{
		fLayerEDep.clear();
}

//=============================================================================

void FluenceEnergyDistributionSD::EndOfEvent(G4HCofThisEvent*)
{
	if(fIntEdep > 0.) {
	G4int evtID = fEventManager->GetConstCurrentEvent()->GetEventID();
	fCountsMap[evtID]=fLayerEDep;
	}
	clear();
	fIntEdep=0.;
}

//=============================================================================

G4double FluenceEnergyDistributionSD::ComputeVolume(G4Step* step, G4int idx)
{
 G4VPhysicalVolume* physVol = step->GetPreStepPoint()->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid = 0;
  if(physParam)
  { // for parameterized volume
    if(idx<0)
    {
      G4ExceptionDescription ED;
      ED << "Incorrect replica number --- GetReplicaNumber : " << idx << G4endl;
      G4Exception("G4PSCellFlux::ComputeVolume","DetPS0001",JustWarning,ED);
    }
    solid = physParam->ComputeSolid(idx, physVol);
    solid->ComputeDimensions(physParam,idx,physVol);
  }
  else
  { // for ordinary volume
    solid = physVol->GetLogicalVolume()->GetSolid();
  }
  return solid->GetCubicVolume();
}

//=============================================================================

inline G4int FluenceEnergyDistributionSD::GetIndex(G4Step* aStep)
{
	const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
	G4int i = touchable->GetReplicaNumber(fDepthi);
	G4int j = touchable->GetReplicaNumber(fDepthj);
	G4int k = touchable->GetReplicaNumber(fDepthk);

	if(i<0||j<0||k<0)
	{
		G4ExceptionDescription ED;
		ED << "GetReplicaNumber is negative" << G4endl
		   << "touchable->GetReplicaNumber(fDepthi) returns i,j,k = "
		   << i << "," << j << "," << k << " for volume "
		   << touchable->GetVolume(fDepthi)->GetName() << ","
		   << touchable->GetVolume(fDepthj)->GetName() << ","
		   << touchable->GetVolume(fDepthk)->GetName() << G4endl;
		  G4Exception("G4PSEnergyDeposit3D::GetIndex","DetPS0006",
													JustWarning,ED);
	}

	return i*fNj*fNk+j*fNk+k;
}

//=============================================================================

void FluenceEnergyDistributionSD::DumpAllDetectorCollects(

                                                                                        std::fstream& out)
{
    //out.open();
    out.precision(6);


    //------- Energy Deposition will be printed in output folder, for both cases of 100 layers and 8 (1+7)
if(conf()->SiLayersDep==1)
{
	CountMap::iterator it    = fCountsMap.begin();
	CountMap::iterator itend = fCountsMap.end();
	out << "Tot Events:" << G4endl;
        out << fCountsMap.size() << G4endl;
	for(;it!=itend;++it) {
out << "Event ID with hit" << G4endl;
		out << std::scientific << it->first << "\t";
		G4int i=0, j=0, k=0;
		G4double tot=0;
		for(;i<fNi;i++){
	        for(;j<fNj;j++)
		    for(;k<10;k++) {
	                G4int ind = CalcIndex(i,j,k);
					out << (it->second)[ind] << "\t";
					tot=+(it->second)[ind];
	            }
		}
		out <<"tot E"<< tot << G4endl;
	}
}




}

