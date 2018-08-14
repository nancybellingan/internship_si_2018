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
	// get the pointer pSDman for the sensitive scorers
	G4SDManager* pSDman = G4SDManager::GetSDMpointer();


	//if the spherescorer is enabled, resize the vector containing all the IDs of the scorers
	//and get the ID for each of them via pSDman
	if(conf()->SphereScorer==1){
		SphereFluxID.resize(conf()->ebin.size());
		for (uint i=0;i<conf()->ebin.size();i++){
			G4String evtID =std::to_string(i);

			SphereFluxID[i] = pSDman->GetCollectionID(detectorName.pathSphere+evtID);
		}
	}
	//if the DummyScorer is enabled, resize the vector containing all the IDs of the scorers
	//and get the ID for each of them via pSDman, for both fast and albedo
	if(conf()->DummyScorer==1){
		FastFluxID.resize(conf()->ebin.size());
		AlbedoFluxID.resize(conf()->ebin.size());
		for (uint i=0;i<conf()->ebin.size();i++){
			G4String evtID =std::to_string(i);
			FastFluxID[i] = pSDman->GetCollectionID(detectorName.pathFast+evtID);
			AlbedoFluxID[i] = pSDman->GetCollectionID(detectorName.pathAlbedo+evtID);

		}
	}
	// if the energy deposition is wanted, get the ID for both fast and albedo
	if(conf()->SiLayersDep==1){
		EDepFastID = pSDman->GetCollectionID("fastDetDep/EDepFast");
		EDepAlbedoID = pSDman->GetCollectionID("albedoDetDep/EDepAlbedo");
	}

	//=================================================
	//  Initalize RunMaps for accumulation.
	//  Get CollectionIDs for HitCollections.

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
	//
	// Destructor
	//    clear all data members.
}
Run::~Run()
{
}

// elements needed for the multithreading (mutable variable and the counterforfileflush)
#include <mutex>
static std::mutex mutexFileWrite;
static std::mutex mutexFileWrite2;
static std::mutex mutexFileWrite3;
static std::mutex mutexFileWrite4;
static int counterForFlileFlush = 0;
static int counterForFlileFlush2 = 0;
static int counterForFlileFlush3 = 0;
static int counterForFlileFlush4 = 0;

// print function for the number of neutrons per binning for the sphere and dummy scorers
// arguments to pass are the hitsmap, the file where to save, the binning slot aka number
// of scorer and the event
void printOnHit(const G4VHitsCollection* eventMap, std::ofstream* file, const uint binSlot, const G4Event* aEvent ){
	const G4THitsMap<G4double>* castedMap = (const G4THitsMap<G4double>*) eventMap;
	auto map = castedMap->GetMap();
	if(! map->empty()){
		//	*file << binSlot << "prova" << G4endl;
		//We only care about the FIRST (key = 0) layer
		//in case we have > 1 event detected, we print multiple time
		auto iter = map->find(0);
		if( iter != map->end()){
			counterForFlileFlush++;
			auto eventRegistered = *iter->second;
			std::lock_guard<std::mutex> lock(mutexFileWrite);
			// for each time there is a scoring, it prints the number of the binning slot
			// for future histogram
			for(int i=0; i < eventRegistered; i++){
				*file << "EventID= \t" << aEvent->GetEventID() << "\t Bin number \t" << binSlot << G4endl;
			}
		}
	}
}
// here i tried to do the same done for PrintonHit, but adapted to this kind of map which should have 40 different
//first values one for each division of the silicon part, and then as second the energy deposition. i need it over
// the first 10 elements actually, not all.
//probably implement a parallel variable with the sum with the correction factor
void PrintOnDep(const G4VHitsCollection* eventMap, std::ofstream* file, const G4Event* aEvent){
	const G4THitsMap<G4double>* castedMap = (const G4THitsMap<G4double>*) eventMap;
	std::map<G4int,G4double*>::iterator it = castedMap->GetMap()->begin();
	counterForFlileFlush2++;
	std::lock_guard<std::mutex> lock(mutexFileWrite2);
	int index = 0;
	for(; it != castedMap->GetMap()->end(); it++){
		auto edep = *it->second;

		if(edep>0){
			index++;

			*file << "fNz number for event ID \t" <<  aEvent->GetEventID() << "\t = \t"<< it->first << "\t E Dep in MeV = \t" << edep << "\t \t";

		}
	}
	if(index>0){
		*file << "\n"  << G4endl;
	}

}
// this function is to print the total energy deposited by an event over the first 10 segments
// weighted with the correction factor
void PrintDeptot(const G4VHitsCollection* eventMap, std::ofstream* file, const G4Event* aEvent){
	const G4THitsMap<G4double>* castedMap = (const G4THitsMap<G4double>*) eventMap;
	std::map<G4int,G4double*>::iterator it = castedMap->GetMap()->begin();
	counterForFlileFlush3++;
	std::lock_guard<std::mutex> lock(mutexFileWrite3);
	int index = 0;
	double etot = 0;
	for(; it != castedMap->GetMap()->end(); it++){
		if (it->first<10){
			auto edep = *it->second;
			index++;
			//cumulative etot
			etot = etot + edep*correctionfactor.at(it->first);
		}
	}
	if (etot > 0){
		*file << "event ID = \t" << aEvent->GetEventID() << "\t Edeptot (MeV)= \t" << etot << G4endl ;
	}

}

void PrintDeptotfiltered(const G4VHitsCollection* eventMap, std::ofstream* file, const G4Event* aEvent){
	const G4THitsMap<G4double>* castedMap = (const G4THitsMap<G4double>*) eventMap;
	std::map<G4int,G4double*>::iterator it = castedMap->GetMap()->begin();
	counterForFlileFlush4++;
	std::lock_guard<std::mutex> lock(mutexFileWrite4);
	int index = 0;
	double etot = 0;
	for(; it != castedMap->GetMap()->end(); it++){
		if (it->first<10){
			auto edep = *it->second;
			index++;
			//cumulative etot
			etot = etot + edep*correctionfactor.at(it->first);
		}
	}
	if (etot > 0.15){
		*file << "event ID = \t" << aEvent->GetEventID() << "\t Edeptot (MeV)= \t" << etot << G4endl ;
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  RecordEvent is called at end of event.
//  For scoring purpose, the resultant quantity in a event,
//  is accumulated during a Run.
void Run::RecordEvent(const G4Event* aEvent) {
	numberOfEvent++;  // This is an original line.
	//=============================
	// HitsCollection of This Event
	//============================
	aEvent->GetEventID();

	G4HCofThisEvent* pHCE = aEvent->GetHCofThisEvent();
	if (!pHCE) {
		return;
	}

	//when you switch to the other sensor style, change plug in the other function and gg...
	if(conf()->SphereScorer==1){

		eventSphereFlux.resize(conf()->ebin.size());
		totSphereFlux.resize(conf()->ebin.size());

		for (uint binSlot=0;binSlot<conf()->ebin.size();binSlot++){
			if(totSphereFlux[binSlot] == nullptr){

				totSphereFlux[binSlot] = new G4THitsMap<G4double>();
			}
			eventSphereFlux[binSlot] = (G4THitsMap<G4double>*)(pHCE->GetHC(SphereFluxID[binSlot]));
			//			printOnHit(pHCE->GetHC(SphereFluxID[binSlot]),conf()->SphereFlux,binSlot, aEvent);
			printOnHit(eventSphereFlux[binSlot],conf()->SphereFlux,binSlot, aEvent);

			if (conf()->totdata==1){
				if(totSphereFlux[binSlot] == nullptr){

					totSphereFlux[binSlot] = new G4THitsMap<G4double>();
				}
			*totSphereFlux[binSlot] += *eventSphereFlux[binSlot];
			}

		}
	}
	if(conf()->DummyScorer==1){
		//Here the sensor has 8 LAYER, but we do not care
		eventAlbedoFlux.resize(conf()->ebin.size());
		totAlbedoFlux.resize(conf()->ebin.size());
		eventFastFlux.resize(conf()->ebin.size());
		totFastFlux.resize(conf()->ebin.size());

		for (uint binSlot=0;binSlot<conf()->ebin.size();binSlot++){
			if(totAlbedoFlux[binSlot] == nullptr){
				totAlbedoFlux[binSlot] = new G4THitsMap<G4double>();
			}
			if(totFastFlux[binSlot] == nullptr){
				totFastFlux[binSlot] = new G4THitsMap<G4double>();
			}
			eventAlbedoFlux[binSlot] = (G4THitsMap<G4double>*)(pHCE->GetHC(AlbedoFluxID[binSlot]));
			eventFastFlux[binSlot] = (G4THitsMap<G4double>*)(pHCE->GetHC(FastFluxID[binSlot]));

			printOnHit(eventFastFlux[binSlot], conf()->fastFlux,binSlot, aEvent);
			printOnHit(eventAlbedoFlux[binSlot], conf()->albedoFlux,binSlot, aEvent);

			if (conf()->totdata==1){
				if(totAlbedoFlux[binSlot] == nullptr){
					totAlbedoFlux[binSlot] = new G4THitsMap<G4double>();
				}
				if(totFastFlux[binSlot] == nullptr){
					totFastFlux[binSlot] = new G4THitsMap<G4double>();
				}
				*totAlbedoFlux[binSlot] += *eventAlbedoFlux[binSlot];
				*totFastFlux[binSlot] += *eventFastFlux[binSlot];
			}
		}
	}
	if(conf()->SiLayersDep==1){
		auto depMapFast = (G4THitsMap<G4double>*)(pHCE->GetHC(EDepFastID));
		auto depMapAlbedo = (G4THitsMap<G4double>*)(pHCE->GetHC(EDepAlbedoID));

		PrintOnDep(depMapFast,conf()->fastDep, aEvent);
		PrintOnDep(depMapAlbedo,conf()->albedoDep, aEvent);
		PrintDeptot(depMapFast,conf()->fastTotDep, aEvent);
		PrintDeptot(depMapAlbedo,conf()->albedoTotDep, aEvent);
		PrintDeptotfiltered(depMapAlbedo,conf()->albedoTotDepfilt, aEvent);
		PrintDeptotfiltered(depMapFast,conf()->fastTotDepfilt, aEvent);
	}
	if(counterForFlileFlush > 128){
		std::lock_guard<std::mutex> lock(mutexFileWrite);
		conf()->SphereFlux->flush();
		conf()->fastFlux->flush();
		conf()->albedoFlux->flush();

		counterForFlileFlush = 0;
	}
	if(counterForFlileFlush2 > 128){
		std::lock_guard<std::mutex> lock(mutexFileWrite2);
		conf()->albedoDep->flush();
		conf()->fastDep->flush();
		counterForFlileFlush2 = 0;
	}
	if(counterForFlileFlush3 > 128){
		std::lock_guard<std::mutex> lock(mutexFileWrite3);
		conf()->fastTotDep->flush();
		conf()->albedoTotDep->flush();
		counterForFlileFlush3 = 0;
	}
	if(counterForFlileFlush4 > 128){
		std::lock_guard<std::mutex> lock(mutexFileWrite4);
		conf()->albedoTotDepfilt->flush();
		conf()->fastTotDepfilt->flush();
		counterForFlileFlush4 = 0;
	}



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
// Roy this is the part doesnt work in multi threading if i have the tot...flux on.
	const Run* localRun = static_cast<const Run*>(aRun);
	//=======================================================
	// Merge HitsMap of working threads
	//=======================================================
	//bool a = sameSize(totSphereFlux,localRun->totSphereFlux);
	//bool b = sameSize(totFastFlux,localRun->totFastFlux);
	//bool c = sameSize(totAlbedoFlux,localRun->totAlbedoFlux);
	//if(!(a && b && c)){
//		abort();
	//}
	if(conf()->SphereScorer==1 && conf()->totdata==1){
		static std::mutex mutexFileWrite5;
		for(uint num = 0; num < totSphereFlux.size(); num++){
	std::lock_guard<std::mutex> lock(mutexFileWrite5);
			*totSphereFlux.at(num)  += *localRun->totSphereFlux[num];
		}
	}
	if(conf()->DummyScorer==1 && conf()->totdata==1){
		static std::mutex mutexFileWrite5;
		static std::mutex mutexFileWrite6;
		for(uint num = 0; num < totFastFlux.size(); num++){
			std::lock_guard<std::mutex> lock(mutexFileWrite5);

			*totFastFlux.at(num)  += *localRun->totFastFlux[num];
		}
		for(uint num = 0; num < totAlbedoFlux.size(); num++){
			std::lock_guard<std::mutex> lock(mutexFileWrite6);

			*totAlbedoFlux.at(num)  += *localRun->totAlbedoFlux[num];
		}
	}
	G4Run::Merge(aRun);
}

