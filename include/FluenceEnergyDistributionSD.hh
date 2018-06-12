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
// $Id: FluenceEnergyDistributionSD.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file FluenceEnergyDistributionSD.hh
/// \brief Definition of the FluenceEnergyDistributionSD class

#ifndef FluenceEnergyDistributionSD_h
#define FluenceEnergyDistributionSD_h 1

#include "G4VSensitiveDetector.hh"

#include "FluenceEnergyDistributionHit.hh"
#include "G4EventManager.hh"
#include "fstream"
#include "iomanip"
#include "iostream"
#include <vector>
#include <map>

class G4Step;
class G4HCofThisEvent;

typedef std::map<G4int,G4double> LayerCollect;
typedef std::map<G4int,LayerCollect> CountMap;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// B2Tracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero 
/// energy deposit.

class FluenceEnergyDistributionSD : public G4VSensitiveDetector
{
  public:
	  FluenceEnergyDistributionSD(const G4String& name,
									const G4int& ni=1,
									const G4int& nj=1,
									const G4int& nk=1,
									const G4int& depth=0,
									const G4int& depi=2,
									const G4int& depj=1,
									const G4int& depk=0);
    virtual ~FluenceEnergyDistributionSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);
        void DumpAllDetectorCollects(std::fstream& out);
 //       void DumpAllDetectorCollects(std::ofstream& out);
	virtual G4double ComputeVolume(G4Step* step,G4int idx);
	virtual void SetDepth(const G4int& depth) { fDepth = depth; };
	virtual void SetWeighted(const G4bool& weighted) { fWeighted = weighted; };

  protected:
	inline virtual G4int GetIndex(G4Step*);

  private:
	inline G4int CalcIndex(G4int& i, G4int& j, G4int& k)
		                        { return i*fNj*fNk+j*fNk+k; };
	inline void clear();

  private:
	G4int  fNi, fNj, fNk;
	G4int  fDepth;
	G4int  fDepthi, fDepthj, fDepthk;
	G4bool fWeighted;
	LayerCollect fLayerEDep;
	CountMap fCountsMap;
	G4EventManager* fEventManager;
	G4double fIntEdep;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
