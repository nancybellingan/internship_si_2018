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
/// \file runAndEvent/RE02/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.hh 75682 2013-11-05 09:11:19Z gcosmo $
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MultiFunctionalDetector.hh"
#include <vector>
#include <G4PSPassageCellCurrent.hh>
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4VPrimitiveScorer.hh"
class G4Box;
class G4Orb;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVPlacement;
class G4Material;
class G4SubtractionSolid;
class G4Sphere;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

	DetectorConstruction();
	virtual ~DetectorConstruction();

public:
		  // virtual method from G4VUserDetectorConstruction.
	virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();

public:
  // Get/Set Access methods for data members
	void GetNumberOfSegmentsInPhantom(G4int& nx, G4int& ny, G4int& nz) const
	                                                    { nx=fNx; ny=fNy; nz=fNz;}

	std::vector<G4double> GetEnergyBinning() const { return fEnergyBinning; };
	std::vector<G4String> GetUSDNames()      const { return fUSDNames;      };
	std::vector<G4String> GetMFDNames()      const { return fMFDNames;      };
	std::vector<G4String> GetUSDParticles()  const { return fUSDParticles;  };
	std::vector<G4String> GetMFDParticles()  const { return fMFDParticles;  };
	std::vector<G4LogicalVolume*> GetUSDVolumes() const { return fUSDVolumes; };
	std::vector<G4LogicalVolume*> GetMFDVolumes() const { return fMFDVolumes; };
//	G4VPhysicalVolume GetScorer() const {return physicscorer; };
private:
  // Data members
	G4bool fCheckOverlaps;
	G4int  fNx,fNy,fNz;  // Number of segmentation of water phantom.
	G4double numberOfLayers, numberOfLayers2;
	G4LogicalVolume* LogicFastSi;
	G4LogicalVolume* LogicAlbedoSi= nullptr;
	G4LogicalVolume* fast_housing_log;
	G4LogicalVolume* albedo_housing_log;

	// Data for all Threads
	//==========================================================================

	std::vector<G4double> fEnergyBinning;
	std::vector<G4String> fUSDNames;
	std::vector<G4String> fMFDNames;
	std::vector<G4String> fUSDParticles;
	std::vector<G4String> fMFDParticles;
	std::vector<G4LogicalVolume*> fUSDVolumes;
	std::vector<G4LogicalVolume*> fMFDVolumes;
//	std::vector<G4PSPassageCellCurrent*> totalSphereFlux;
//	std::vector<G4SDParticleWithEnergyFilter*> binfilter;
//	std::vector<G4PSPassageCellCurrent*> fastflux;
//	std::vector<G4SDParticleWithEnergyFilter*> fastbinfilter;
	//std::vector<G4PSPassageCellCurrent*> albedoflux;
//	std::vector<G4SDParticleWithEnergyFilter*> albedobinfilter;

	G4Material* Vacuum;
	G4Material* fAir;
	G4Material* Air;
	G4Material* fWater;
	G4Material* fPMMA;
	G4Element* fO;
	G4Element* fthsc_H;
	G4Material* Concrete;

	G4PVPlacement* fast_leadFront;
	G4PVPlacement* fast_wax;
	G4PVPlacement* fast_ceramic;
	G4PVPlacement* fast_leadBack;
	G4PVPlacement* physiFastSens;
	G4ThreeVector Fast_housing_pos;

	G4PVPlacement* albedo_hullFront;
	G4PVPlacement* albedo_gap;
	G4PVPlacement* albedo_converter;
	G4PVPlacement* albedo_ceramic;
	G4PVPlacement* albedo_hullBack;
	G4PVPlacement* albedo_hole;
	G4PVPlacement* physiAlbedoSens;
	G4ThreeVector albedo_housing_pos;

	//==========================================================================
	G4LogicalVolume* fast_leadFront_log;
	G4LogicalVolume* albedo_hullFront_log;
	G4Box *fSolidPMMAPhantom;
	G4LogicalVolume* fLogicPMMAPhantom;
	G4VPhysicalVolume* fPhysiPMMAPhantom;

	G4Box* outerBox;
	G4Box* innerBox;
	G4SubtractionSolid* solidroom;
	G4LogicalVolume* logicroom;
	G4PVPlacement* physicroom;

	G4Box *fSolidPMMAPhantomScorer;
	G4LogicalVolume* fLogicPMMAPhantomScorer;
	G4VPhysicalVolume* fPhysiPMMAPhantomScorer;

	G4Box *solidcolumn;
	G4LogicalVolume* logiccolumn;
	G4PVPlacement* physiccolumn;

	G4Box *solidfloor;
	G4LogicalVolume* logicfloor;
	G4PVPlacement* physicfloor;

	G4Sphere *solidscorer;
	G4LogicalVolume* logicscorer;
	G4VPhysicalVolume* physicscorer;

	//==========================================================================
};

#endif
