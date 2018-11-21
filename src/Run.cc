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
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include <G4VHitsCollection.hh>
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "config.h"
#include <string>
#include "G4MTRunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//#define MFD 1

//  Constructor.
//   (The vector of MultiFunctionalDetector name has to given.)
Run::Run() : G4Run()
{
	// get the pointer pSDman for the sensitive scorers



	G4SDManager* pSDman = G4SDManager::GetSDMpointer();
	totSphereFluxv2.resize(conf()->ebin.size());
	totAlbedoFluxv2.resize(conf()->ebin.size());
	totFastFluxv2.resize(conf()->ebin.size());
	totPhantomFluxv2.resize(conf()->ebin.size());
	totbackalbedov2.resize(conf()->ebin.size());
	totFrontalPhantomv2.resize(conf()->ebin.size());

	partialsumAlbedoFluxv2.resize(conf()->ebin.size());
	partialsumFastFluxv2.resize(conf()->ebin.size());
	partialsumPhantomFluxv2.resize(conf()->ebin.size());
	partialsumSphereFluxv2.resize(conf()->ebin.size());
	partialsumFrontalPhantomv2.resize(conf()->ebin.size());
	partialsumbackalbedov2.resize(conf()->ebin.size());

	    //if the spherescorer is enabled, resize the vector containing all the IDs of the scorers
	//and get the ID for each of them via pSDman
	if(conf()->SphereScorer==1){
		SphereFluxID.resize(conf()->ebin.size());
		for (uint i=0;i<conf()->ebin.size();i++){
			G4String evtID =std::to_string(i);
			
			SphereFluxID[i] = pSDman->GetCollectionID(detectorName.pathSphere+evtID);
		}
	}

	if (conf()->backflux==1){
	backalbedoID.resize(conf()->ebin.size());
	PhantomFrontalID.resize(conf()->ebin.size());
	for (uint i=0;i<conf()->ebin.size();i++){
		G4String evtID =std::to_string(i);
		PhantomFrontalID[i] = pSDman->GetCollectionID("frontalphantom/phantomfrontalflux"+evtID);
		backalbedoID[i] = pSDman->GetCollectionID("backalbedo/backalbedoflux"+evtID);
	    }
	}
	if(conf()->phantomscorer ==1){
		PhantomFluxID.resize(conf()->ebin.size());
		for (uint i=0;i<conf()->ebin.size();i++){
			G4String evtID =std::to_string(i);
			
			PhantomFluxID[i] = pSDman->GetCollectionID("phantomscorer/totalphantomFlux"+evtID);
		}
	}

	if (conf()->faston==1){

		//if the DummyScorer is enabled, resize the vector containing all the IDs of the scorers
		//and get the ID for each of them via pSDman, for both fast and albedo
		if(conf()->DummyScorer==1){
			FastFluxID.resize(conf()->ebin.size());
			for (uint i=0;i<conf()->ebin.size();i++){
				G4String evtID =std::to_string(i);
				FastFluxID[i] = pSDman->GetCollectionID(detectorName.pathFast+evtID);

			}
		}
		// if the energy deposition is wanted, get the ID for both fast and albedo
		if(conf()->SiLayersDep==1){
			EDepFastID = pSDman->GetCollectionID("fastDetDep/EDepFast");
		}

	}
	if (conf()->albedoon ==1){


		//if the DummyScorer is enabled, resize the vector containing all the IDs of the scorers
		//and get the ID for each of them via pSDman, for both fast and albedo
		if(conf()->DummyScorer==1){
			AlbedoFluxID.resize(conf()->ebin.size());
			for (uint i=0;i<conf()->ebin.size();i++){
				G4String evtID =std::to_string(i);
				AlbedoFluxID[i] = pSDman->GetCollectionID(detectorName.pathAlbedo+evtID);

			}
		}
		// if the energy deposition is wanted, get the ID for both fast and albedo
		if(conf()->SiLayersDep==1){
			EDepAlbedoID = pSDman->GetCollectionID("albedoDetDep/EDepAlbedo");
		}

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
void Run::printOnHitsphere(const G4VHitsCollection* eventMap, std::ofstream* file, const uint binSlot, const G4Event* aEvent ){
	const G4THitsMap<G4double>* castedMap = (const G4THitsMap<G4double>*) eventMap;
	auto map = castedMap->GetMap();
	if(! map->empty()){
		auto iter = map->find(0);
		partialsumSphereFluxv2[binSlot] = 0;
		counterForFlileFlush++;
		std::lock_guard<std::mutex> lock(mutexFileWrite);
		//		auto itr = castedMap->GetMap()->begin();
		//	for (; itr != castedMap->GetMap()->end();itr++ ){
		partialsumSphereFluxv2[binSlot] = partialsumSphereFluxv2[binSlot] + ((int) lround(*(iter->second)));
		//	}
		*file << "EventID= \t" << aEvent->GetEventID() << "\t Bin number \t" << binSlot << "\t number of events: \t" << partialsumSphereFluxv2[binSlot] << G4endl;
	}
}
void Run::printonHitbackalbedo(const G4VHitsCollection* eventMap, std::ofstream* file, const uint binSlot, const G4Event* aEvent ){
	const G4THitsMap<G4double>* castedMap = (const G4THitsMap<G4double>*) eventMap;
	auto map = castedMap->GetMap();
	if(!map->empty()){
		auto iter = map->find(0);
		counterForFlileFlush++;
		partialsumbackalbedov2[binSlot] =0;
		    std::lock_guard<std::mutex> lock(mutexFileWrite2);
			partialsumbackalbedov2[binSlot] = partialsumbackalbedov2[binSlot] + ((int) lround(*(iter->second)));
			    *file << "EventID= \t" << aEvent->GetEventID() << "\t Bin number \t" << binSlot << "\t number of events: \t"
			          << partialsumbackalbedov2[binSlot] << G4endl;
	//	partialsumFrontalFastv2[binSlot] = 0;
	//	counterForFlileFlush++;
	//	std::lock_guard<std::mutex> lock(mutexFileWrite);
	//	auto itr = castedMap->GetMap()->begin();
	//	for (; itr != castedMap->GetMap()->end();itr++ ){
	//		partialsumFrontalFastv2[binSlot] = partialsumFrontalFastv2[binSlot] + ((int) lround((*itr->second)));
	//	}
	//	*file << "EventID= \t" << aEvent->GetEventID() << "\t Bin number \t" << binSlot << "\t number of events: \t" << partialsumFrontalFastv2[binSlot] << G4endl;

	}

}
void Run::printonHitfrontalphantom(const G4VHitsCollection* eventMap, std::ofstream* file, const uint binSlot, const G4Event* aEvent ){
	const G4THitsMap<G4double>* castedMap = (const G4THitsMap<G4double>*) eventMap;
	auto map = castedMap->GetMap();
	if(! map->empty()){
		auto iter = map->find(0);
		    std::lock_guard<std::mutex> lock(mutexFileWrite3);
		    counterForFlileFlush++;
			partialsumFrontalPhantomv2[binSlot] =0;
			partialsumFrontalPhantomv2[binSlot] = partialsumFrontalPhantomv2[binSlot] + ((int) lround(*(iter->second)));

		*file << "EventID= \t" << aEvent->GetEventID() << "\t Bin number \t" << binSlot << "\t number of events: \t"
			          << partialsumFrontalPhantomv2[binSlot] << G4endl;

	}
}

void Run::printOnHitphantom (const G4VHitsCollection* eventMap, std::ofstream* file, const uint binSlot, const G4Event* aEvent ){
	const G4THitsMap<G4double>* castedMap = (const G4THitsMap<G4double>*) eventMap;
	auto map = castedMap->GetMap();
	if(! map->empty()){

		partialsumPhantomFluxv2[binSlot] = 0;
		counterForFlileFlush++;
		std::lock_guard<std::mutex> lock(mutexFileWrite);
		auto itr = castedMap->GetMap()->begin();
		for (; itr != castedMap->GetMap()->end();itr++ ){
			partialsumPhantomFluxv2[binSlot] = partialsumPhantomFluxv2[binSlot] + ((int) lround((*itr->second)));
		}
		*file << "EventID= \t" << aEvent->GetEventID() << "\t Bin number \t" << binSlot << "\t number of events: \t" << partialsumPhantomFluxv2[binSlot] << G4endl;
	}

}
void Run::printOnHitfast(const G4VHitsCollection* eventMap, std::ofstream* file, const uint binSlot, const G4Event* aEvent ){
	const G4THitsMap<G4double>* castedMap = (const G4THitsMap<G4double>*) eventMap;
	auto map = castedMap->GetMap();
	if(! map->empty()){
		auto iter = map->find(0);
		if( iter != map->end()){
			counterForFlileFlush++;
			partialsumFastFluxv2[binSlot] =0;
			//	auto eventRegistered = castedMap->second;
			//	auto key = *castedMap->first;
			std::lock_guard<std::mutex> lock(mutexFileWrite4);
			//	auto itr = castedMap->GetMap()->begin();
			// for (; itr != castedMap->GetMap()->end();itr++ ){
			partialsumFastFluxv2[binSlot] = partialsumFastFluxv2[binSlot] + ((int) lround(*(iter->second)));
		}
		*file << "EventID= \t" << aEvent->GetEventID() << "\t Bin number \t" << binSlot << "\t number of events: \t"
		      << partialsumFastFluxv2[binSlot] << G4endl;
		// for each time there is a scoring, it prints the number of the binning slot
		// for future histogram
		//		for(int i=0; i<eventRegistered; i++){
		//			*file << "EventID= \t" << aEvent->GetEventID() << "\t Bin number \t" << binSlot << G4endl;
		//		}

	}


}

void Run::printOnHitalbedo(const G4VHitsCollection* eventMap, std::ofstream* file, const uint binSlot, const G4Event* aEvent ){
	const G4THitsMap<G4double>* castedMap = (const G4THitsMap<G4double>*) eventMap;
	auto map = castedMap->GetMap();
	if(! map->empty()){
		auto iter = map->find(0);
		if( iter != map->end()){
			counterForFlileFlush++;
			partialsumAlbedoFluxv2[binSlot] =0;
			//	auto eventRegistered = castedMap->second;
			//	auto key = *castedMap->first;

			//	auto itr = castedMap->GetMap()->begin();
			// for (; itr != castedMap->GetMap()->end();itr++ ){
			partialsumAlbedoFluxv2[binSlot] = partialsumAlbedoFluxv2[binSlot] + ((int) lround(*(iter->second)));
		}
		*file << "EventID= \t" << aEvent->GetEventID() << "\t Bin number \t" << binSlot << "\t number of events: \t" << partialsumAlbedoFluxv2[binSlot] << G4endl;
		// for each time there is a scoring, it prints the number of the binning slot
		// for future histogram
		//		for(int i=0; i<eventRegistered; i++){
		//			*file << "EventID= \t" << aEvent->GetEventID() << "\t Bin number \t" << binSlot << G4endl;
		//		}
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
// print function for the number of neutrons per binning for the sphere and dummy scorers
// arguments to pass are the hitsmap, the file where to save, the binning slot aka number
// of scorer and the event

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
		
		for (uint binSlot=0;binSlot<conf()->ebin.size();binSlot++){
			mutexFileWrite2.lock();

			eventSphereFlux[binSlot] = (G4THitsMap<G4double>*)(pHCE->GetHC(SphereFluxID[binSlot]));
			printOnHitsphere(eventSphereFlux[binSlot],conf()->SphereFlux,binSlot, aEvent);
			mutexFileWrite2.unlock();

			if (conf()->totdata==1){
				if((!totSphereFluxv2[binSlot])){
					totSphereFluxv2[binSlot] = 0;
				}
				if(partialsumSphereFluxv2[binSlot] && partialsumSphereFluxv2[binSlot]!=0 ){
					totSphereFluxv2[binSlot] = totSphereFluxv2[binSlot] + partialsumSphereFluxv2[binSlot];
			}
			}
		}
		std::fill(partialsumSphereFluxv2.begin(), partialsumSphereFluxv2.end(), 0);
	}
	if (conf()->phantomscorer==1){
		eventPhantomFlux.resize(conf()->ebin.size());
		for (uint binSlot=0;binSlot<conf()->ebin.size();binSlot++){
			eventPhantomFlux[binSlot] = (G4THitsMap<G4double>*)(pHCE->GetHC(PhantomFluxID[binSlot]));
			mutexFileWrite4.lock();
			printOnHitphantom(eventPhantomFlux[binSlot],conf()->phantomFlux,binSlot, aEvent);
   mutexFileWrite4.unlock();
            if (conf()->totdata==1){
				if((!totPhantomFluxv2[binSlot])){
					totPhantomFluxv2[binSlot] = 0;
				}
				if(partialsumPhantomFluxv2[binSlot] && partialsumPhantomFluxv2[binSlot]!=0 ){
					totPhantomFluxv2[binSlot] = totPhantomFluxv2[binSlot] + partialsumPhantomFluxv2[binSlot];
			}
			}

		}
		std::fill(partialsumPhantomFluxv2.begin(), partialsumPhantomFluxv2.end(), 0);


	    }


if (conf()->faston==1){
	if(conf()->DummyScorer==1){

		eventFastFlux.resize(conf()->ebin.size());
		for (uint binSlot=0;binSlot<conf()->ebin.size();binSlot++){
			mutexFileWrite.lock();
			eventFastFlux[binSlot] = (G4THitsMap<G4double>*)(pHCE->GetHC(FastFluxID[binSlot]));
			printOnHitfast(eventFastFlux[binSlot], conf()->fastFlux,binSlot, aEvent);
mutexFileWrite.unlock();

			if (conf()->totdata==1){
				if(!totFastFluxv2[binSlot]){
					totFastFluxv2[binSlot] = 0;
				}
				if(partialsumFastFluxv2[binSlot] && partialsumFastFluxv2[binSlot]!=0 ){
					totFastFluxv2[binSlot] = totFastFluxv2[binSlot] + partialsumFastFluxv2[binSlot];
			}

			}

		}
		    std::fill(partialsumFastFluxv2.begin(), partialsumFastFluxv2.end(), 0);
	}



	if(conf()->SiLayersDep==1){
		auto depMapFast = (G4THitsMap<G4double>*)(pHCE->GetHC(EDepFastID));
mutexFileWrite2.lock();
		PrintOnDep(depMapFast,conf()->fastDep, aEvent);
		PrintDeptot(depMapFast,conf()->fastTotDep, aEvent);
		PrintDeptotfiltered(depMapFast,conf()->fastTotDepfilt, aEvent);
		mutexFileWrite2.unlock();

	}
}


if (conf()->albedoon==1){
	if(conf()->DummyScorer==1){
		eventAlbedoFlux.resize(conf()->ebin.size());
		for (uint binSlot=0;binSlot<conf()->ebin.size();binSlot++){
			eventAlbedoFlux[binSlot] = (G4THitsMap<G4double>*)(pHCE->GetHC(AlbedoFluxID[binSlot]));
mutexFileWrite.lock();
			printOnHitalbedo(eventAlbedoFlux[binSlot], conf()->albedoFlux,binSlot, aEvent);
mutexFileWrite.unlock();

			if (conf()->totdata==1){
				if(!totAlbedoFluxv2[binSlot]){
					totAlbedoFluxv2[binSlot] = 0;
				}
				if(partialsumAlbedoFluxv2[binSlot] && partialsumAlbedoFluxv2[binSlot]!=0 ){
					totAlbedoFluxv2[binSlot] = totAlbedoFluxv2[binSlot] + partialsumAlbedoFluxv2[binSlot];
			}
			}

		}
		    std::fill(partialsumAlbedoFluxv2.begin(), partialsumAlbedoFluxv2.end(), 0);
	}



	if(conf()->SiLayersDep==1){
		auto depMapAlbedo = (G4THitsMap<G4double>*)(pHCE->GetHC(EDepAlbedoID));
	mutexFileWrite4.lock();
		PrintOnDep(depMapAlbedo,conf()->albedoDep, aEvent);
		PrintDeptot(depMapAlbedo,conf()->albedoTotDep, aEvent);
		PrintDeptotfiltered(depMapAlbedo,conf()->albedoTotDepfilt, aEvent);
		    mutexFileWrite4.unlock();
	}

}



        if(conf()->backflux==1){
			eventFrontalPhantom.resize(conf()->ebin.size());
			eventbackalbedo.resize(conf()->ebin.size());

			for (uint binSlot=0;binSlot<conf()->ebin.size();binSlot++){

				eventFrontalPhantom[binSlot] = (G4THitsMap<G4double>*)(pHCE->GetHC(PhantomFrontalID[binSlot]));
				printonHitfrontalphantom(eventFrontalPhantom[binSlot],conf()->phantomfrontal,binSlot, aEvent);
				eventbackalbedo [binSlot] = (G4THitsMap<G4double>*)(pHCE->GetHC(backalbedoID[binSlot]));
				printonHitbackalbedo(eventbackalbedo[binSlot],conf()->backalbedo,binSlot, aEvent);

				if (conf()->totdata==1){
					if((!totbackalbedov2[binSlot])){
						totbackalbedov2[binSlot] = 0;
					}
					if((!totFrontalPhantomv2[binSlot])){
						totFrontalPhantomv2[binSlot] = 0;
					}
					if(partialsumbackalbedov2[binSlot] && partialsumbackalbedov2[binSlot]!=0 ){
						totbackalbedov2[binSlot] = totbackalbedov2[binSlot] + partialsumbackalbedov2[binSlot];
				}
					if(partialsumFrontalPhantomv2[binSlot] && partialsumFrontalPhantomv2[binSlot]!=0 ){
						totFrontalPhantomv2[binSlot] = totFrontalPhantomv2[binSlot] + partialsumFrontalPhantomv2[binSlot];
				}
				}
			}
			std::fill(partialsumbackalbedov2.begin(), partialsumbackalbedov2.end(), 0);
			            std::fill(partialsumFrontalPhantomv2.begin(), partialsumFrontalPhantomv2.end(), 0);
		}




	if(counterForFlileFlush > 128){
		std::lock_guard<std::mutex> lock(mutexFileWrite);
		conf()->SphereFlux->flush();
		conf()->fastFlux->flush();
		conf()->albedoFlux->flush();
		conf()->phantomFlux->flush();

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


/*template <typename T>
bool sameSize(const T& lhs, const T& rhs){
	if(lhs.size() != rhs.size()){
		G4cout << "The size of the vector is different" << __FILE__ << __LINE__;
		return false;
	}
	return true;
}
*/


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run * aRun)
{	    G4cout << "MERGE" << G4endl;
	    const Run* localRun = static_cast<const Run*>(aRun);
		    //=======================================================
		    // Merge HitsMap of working threads
		    //=======================================================



		    if(conf()->SphereScorer==1 && conf()->totdata==1){
				for(uint num = 0; num < conf()->ebin.size(); num++){

					if (localRun->totSphereFluxv2[num]){
	//					G4cout << "merge1" << G4endl;
						if (totSphereFluxv2[num]){
							totSphereFluxv2[num]  += localRun->totSphereFluxv2[num];
						}else{
							totSphereFluxv2[num]  = localRun->totSphereFluxv2[num];
						}
					}
				}
				// fSumDose    += localRun->fSumDose;
			}
			if(conf()->phantomscorer==1 && conf()->totdata==1){
	//			G4cout << "merge2" << G4endl;
				for(uint num = 0; num < conf()->ebin.size(); num++){
					if (localRun->totPhantomFluxv2[num]){
						totPhantomFluxv2[num]  += localRun->totPhantomFluxv2[num];
					}
				}
			}
			if(conf()->DummyScorer==1 && conf()->totdata==1){
				if(conf()->faston==1){
	//			G4cout << "merge3" << G4endl;
					static std::mutex mutexFileWrite5;
				for(uint num = 0; num < conf()->ebin.size(); num++){
					if (localRun->totFastFluxv2[num]){
						totFastFluxv2[num]  += localRun->totFastFluxv2[num];
					}
				}
				}
				if (conf()->albedoon==1){
				for(uint num = 0; num < conf()->ebin.size(); num++){
		//			G4cout << "merge4" << G4endl;
					if (localRun->totAlbedoFluxv2[num]){
						totAlbedoFluxv2[num]  += localRun->totAlbedoFluxv2[num];
					}
				}
			}
			}

			if(conf()->backflux==1 && conf()->totdata==1){
				static std::mutex mutexFileWrite6;
				for(uint num = 0; num < conf()->ebin.size(); num++){
					if (localRun->totbackalbedov2[num]){
						totbackalbedov2[num]  += localRun->totbackalbedov2[num];
					}
				}

				for(uint num = 0; num < conf()->ebin.size(); num++){
		//			G4cout << "merge4" << G4endl;
					if (localRun->totFrontalPhantomv2[num]){
						totFrontalPhantomv2[num]  += localRun->totFrontalPhantomv2[num];
					}
				}

			}


			G4Run::Merge(aRun);
}




