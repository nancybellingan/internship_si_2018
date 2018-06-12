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
/// \file runAndEvent/RE02/src/Run.cc
/// \brief Implementation of the Run class
//
//
// $Id: Run.cc 76292 2013-11-08 13:10:20Z gcosmo $
//

//=====================================================================
//
//  (Description)
//    Run Class is for accumulating scored quantities which is
//  scored using G4MutiFunctionalDetector and G4VPrimitiveScorer.
//  Accumulation is done using G4THitsMap object.
//
//    The constructor Run(const std::vector<G4String> mfdName)
//  needs a vector filled with MultiFunctionalDetector names which
//  was assigned at instantiation of MultiFunctionalDetector(MFD).
//  Then Run constructor automatically scans primitive scorers
//  in the MFD, and obtains collectionIDs of all collections associated
//  to those primitive scorers. Futhermore, the G4THitsMap objects
//  for accumulating during a RUN are automatically created too.
//  (*) Collection Name is same as primitive scorer name.
//
//    The resultant information is kept inside Run objects as
//  data members.
//  std::vector<G4String> fCollName;            // Collection Name,
//  std::vector<G4int> fCollID;                 // Collection ID,
//  std::vector<G4THitsMap<G4double>*> fRunMap; // HitsMap for RUN.
//
//  The resualtant HitsMap objects are obtain using access method,
//  GetHitsMap(..).
//
//=====================================================================

#include "Run.hh"
#include "FluenceEnergyDistributionSD.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include <G4VHitsCollection.hh>
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "config.h"
#include <string>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//#define MFD 1

//  Constructor.
//   (The vector of MultiFunctionalDetector name has to given.)
Run::Run() : G4Run()
{

  G4SDManager* pSDman = G4SDManager::GetSDMpointer();
//voglio dire nell header che la dimensione sarebbe conf ebin.size ma non me lo fa fare mi dava errore)
  SphereFluxID.resize(conf()->ebin.size());
  FastFluxID.resize(conf()->ebin.size());
  AlbedoFluxID.resize(conf()->ebin.size());
  for (uint i=0;i<conf()->ebin.size();i++){
  G4String evtID =std::to_string(i);

  SphereFluxID[i] = pSDman->GetCollectionID(detectorName.pathSphere+evtID);
  FastFluxID[i] = pSDman->GetCollectionID(detectorName.pathFast+evtID);
  AlbedoFluxID[i] = pSDman->GetCollectionID(detectorName.pathAlbedo+evtID);

}
  //=================================================
  //  Initalize RunMaps for accumulation.
  //  Get CollectionIDs for HitCollections.


// G4cout << "Run Constructor" << G4endl;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Destructor
//    clear all data members.
}
Run::~Run()
{

  //--- Clear HitsMap for RUN
  G4int nMap = fRunMap.size();
  for ( G4int i = 0; i < nMap; i++){
    if(fRunMap[i] ) fRunMap[i]->clear();
  }
  fCollName.clear();
  fCollID.clear();
  fRunMap.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  RecordEvent is called at end of event.
//  For scoring purpose, the resultant quantity in a event,
//  is accumulated during a Run.
void Run::RecordEvent(const G4Event* aEvent)
{
	numberOfEvent++;  // This is an original line.
	//	G4cout << "RECORD" << G4endl;
	//=============================
	// HitsCollection of This Event
	//============================
	aEvent->GetEventID();


	G4HCofThisEvent* pHCE = aEvent->GetHCofThisEvent();
	if (!pHCE) return;


	eventSphereFlux.resize(conf()->ebin.size());
	eventFastFlux.resize(conf()->ebin.size());
	eventAlbedoFlux.resize(conf()->ebin.size());
	totSphereFlux.resize(conf()->ebin.size());
	totFastFlux.resize(conf()->ebin.size());
	totAlbedoFlux.resize(conf()->ebin.size());


  for (uint i=0;i<conf()->ebin.size();i++){
	  if(totSphereFlux[i] == nullptr){
		  totSphereFlux[i] = new G4THitsMap<G4double>();
	  }
	  if(totFastFlux[i] == nullptr){
		  totFastFlux[i] = new G4THitsMap<G4double>();
	  }
	  if(totAlbedoFlux[i] == nullptr){
		  totAlbedoFlux[i] = new G4THitsMap<G4double>();
	  }
	eventSphereFlux[i] = (G4THitsMap<G4double>*)(pHCE->GetHC(SphereFluxID[i]));
	eventFastFlux[i] = (G4THitsMap<G4double>*)(pHCE->GetHC(FastFluxID[i]));
	auto cap = pHCE->GetCapacity();
	volatile auto x2 = pHCE->GetNumberOfCollections();
	eventAlbedoFlux[i] = (G4THitsMap<G4double>*)(pHCE->GetHC(AlbedoFluxID[i]));
	*totSphereFlux[i] += *eventSphereFlux[i];
	*totFastFlux[i] += *eventFastFlux[i];
	*totAlbedoFlux[i] += *eventAlbedoFlux[i];
}



  //auto ref = eventSphereFlux->GetMap();
//	for(int i = 0; i < ref->size(); i++){
//		auto ref2 = ref->at(i);
//		G4cout << "pos " << i << " value:" << ref2 << G4endl;
//	}
	/*
	for (auto pair: *(eventSphereFlux->GetMap())){
		G4double flux =*(pair.second);
		G4int copyNb = *(pair.first);
	}
	*/




	G4Run::RecordEvent(aEvent);

}


template <typename T>
bool sameSize(const T& lhs, const T& rhs){
	if(lhs.size() != rhs.size()){
	      G4cout << "The size of the vector is different" << __FILE__ << __LINE__;
		  return false;
	}
	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Run::Merge(const G4Run * aRun)
{

  G4cout << "MERGE" << G4endl;


  const Run* localRun = static_cast<const Run*>(aRun);
  //=======================================================
  // Merge HitsMap of working threads
  //=======================================================
  bool a = sameSize(totSphereFlux,localRun->totSphereFlux);
  bool b = sameSize(totFastFlux,localRun->totFastFlux);
  bool c = sameSize(totAlbedoFlux,localRun->totAlbedoFlux);
  if(!(a && b && c)){
	  abort();
  }
      for(uint num = 0; num < totSphereFlux.size(); num++){
		 *totSphereFlux.at(num)  += *localRun->totSphereFlux[num];
	  }
	  for(uint num = 0; num < totFastFlux.size(); num++){
		 *totFastFlux.at(num)  += *localRun->totFastFlux[num];
	  }
	  for(uint num = 0; num < totAlbedoFlux.size(); num++){
		 *totAlbedoFlux.at(num)  += *localRun->totAlbedoFlux[num];
	  }

  //SphereFlux += *(localRun->SphereFlux);
  //FastFlux += *(localRun->FastFlux);
  //AlbedoFlux += *(localRun->AlbedoFlux);



  G4Run::Merge(aRun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  Access method for HitsMap of the RUN
//
//-----
// Access HitsMap.
//  By  MultiFunctionalDetector name and Collection Name.
//G4THitsMap<G4double>* Run::GetHitsMap(const G4String& detName,
//					 const G4String& colName){
//    G4String fullName = detName+"/"+colName;
//    return GetHitsMap(fullName);
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//-----
// Access HitsMap.
//  By full description of collection name, that is
//    <MultiFunctional Detector Name>/<Primitive Scorer Name>
//G4THitsMap<G4double>* Run::GetHitsMap(){
//    G4int nCol = fCollName.size();
//    for ( G4int i = 0; i < nCol; i++){
//	if ( fCollName[i] == fullName ){
//	    return fRunMap[i];
//	}
//    }
//    return NULL;
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//-----
// - Dump All HitsMap of this RUN. (for debuging and monitoring of quantity).
//   This method calls G4THisMap::PrintAll() for individual HitsMap.
//void Run::DumpAllScorer(  std::fstream& out){


//    out << "prova dentro run" <<G4endl;
//  // - Number of HitsMap in this RUN.
//  G4int n = GetNumberOfHitsMap();
//  // - GetHitsMap and dump values.
//  for ( G4int i = 0; i < n ; i++ ){
//    G4THitsMap<G4double>* runMap =GetHitsMap(i);
//    if ( runMap ) {
//      out << " PrimitiveScorer RUN "
//	     << runMap->GetSDname() <<","<< runMap->GetName() << G4endl;
//      out << " Number of entries " << runMap->entries() << G4endl;
//      std::map<G4int,G4double*>::iterator itr = runMap->GetMap()->begin();
//      for(; itr != runMap->GetMap()->end(); itr++) {
//	out << "  copy no.: " << itr->first
//	       << "  Run Value : " << *(itr->second)
//	       << G4endl;
//      }
//    }
//  }

//}
/**/


//std::vector<G4THitsMap<G4double>*> Run::GetSphereFlux()
//{
//	return totSphereFlux;
//}
