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
/// \file runAndEvent/RE02/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
// $Id: RunAction.cc 75682 2013-11-05 09:11:19Z gcosmo $
// 
#include "RunAction.hh"
#include "Run.hh"

//-- In order to obtain detector information.
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4SDManager.hh"
#include "DetectorConstruction.hh"
#include "FluenceEnergyDistributionSD.hh"
#include "FluenceEnergyDistributionHit.hh"
#include "G4THitsMap.hh"
#include "DetectorCollect.hh"
#include "SDCollect.hh"
#include <sys/stat.h>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"    
#include <chrono>
#include <ctime>
#include <time.h>
#include <fstream>
#include <iostream>
#include "config.h"
#include <string>
//=======================================================================
// RunAction
//  
//
//
//=======================================================================

//#define MFD 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Constructor
RunAction::RunAction()
  : G4UserRunAction(),
    fNx(0), fNy(0), fNz(0)
{
  // - Prepare data member for Run.
  //   vector represents a list of MultiFunctionalDetector names.
    // if the Multi functional detector is activated in the main file, then this macro happens. if not it uses the user detector
/*
#ifdef MFD
	fSDName.push_back(G4String("PhantomSD"));
#else
	fDetector = (DetectorConstruction*)
		          (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
#endif*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Destructor.
RunAction::~RunAction()
{
    // if the debug is activated from main, it is given the information that the destructor is happening

	G4cout << "Destructor RunAction" << G4endl;

  fSDName.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//== 
G4Run* RunAction::GenerateRun()
{
  // Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  //  Detail description can be found in Run.hh/cc.

	return (new Run());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//==
void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//== 
void RunAction::EndOfRunAction(const G4Run* aRun)
{
 Run* re02Run=(Run*)aRun;
  if(!IsMaster()) return;
 // std::time_t timeint = std::time(0);  // t is an integer type
 // G4String timenow = std::to_string(timeint);
  // G4String path = "./";
   auto flux1 = re02Run->GetSphereFlux();
 // pathtime = path + conf()->timenow;
 // mkdir(pathtime, 0777);
  //- Run object.
  //--- Dump all scored quantities involved in Run if is not being used the Multi functional detector
 /* if (conf()->SphereScorer==1){
	  std::fstream outp1;
	  G4int sum = 0;
	  G4String outfiler1=conf()->timenow +"/outputsphere1.dat";
	  //  outfiler.append(std::chrono::system_clock::now());
	  outp1.open(outfiler1.c_str(),std::fstream::out | std::fstream::in  |   std::fstream::trunc);
	  outp1 << "prova" << G4endl;
	  //   G4int n = re02Run->GetNumberOfHitsMap();
	  //  outpr << n << G4endl;
	  //  std::vector<G4THitsMap<G4double>*> ref = re02Run->GetSphereFlux();
	  auto flux1 = re02Run->GetSphereFlux();
	  for (uint i=0;i<flux1.size();i++){
		  auto energy = conf()->ebin[i];
		  for(auto line : *(flux1[i]) ){
			  outp1 << "bin n° " << i <<" energy: " << energy << " id evento ? " << line.first << " boh " << *line.second << G4endl;
			  sum += *line.second;
		  }

		  //gg all
	  }
	  outp1 << "tot events" << sum << G4endl;
	  outp1.close();
  } */
  /*
  if (conf()->DummyScorer==1){
	  std::fstream outp1;
	  G4String outfiler2=pathtime +"/outputfastdummy.dat";
	  G4String outfiler3=pathtime +"/outputalbedodummy.dat";
	  outp1.open(outfiler2.c_str(),std::fstream::out | std::fstream::in  |   std::fstream::trunc);
	  outp1 << "prova" << G4endl;
	  auto flux2 = re02Run->GetFastFlux();
	  for (uint i=0;i<flux2.size();i++){
		  auto energy = conf()->ebin[i];
		  for(auto line : *(flux2[i]) ){
			  outp1 << "bin n° " << i <<" energy: " << energy << " id evento ? " << line.first << " boh " << *line.second << G4endl;
		  }
		  //gg all
	  }
	  outp1.close();
	  outp1.open(outfiler3.c_str(),std::fstream::out | std::fstream::in  |   std::fstream::trunc);
	  outp1 << "prova" << G4endl;
	  auto flux3 = re02Run->GetAlbedoFlux();
	  for (uint i=0;i<flux3.size();i++){
		  auto energy = conf()->ebin[i];
		  for(auto line : *(flux3[i]) ){
			  outp1 << "bin n° " << i <<" energy: " << energy << " id evento ? " << line.first << " boh " << *line.second << G4endl;
		  }
		  //gg all
	  }
	  outp1.close();
  }
*/
/*if(conf()->SiLayersDep ==1){


G4SDManager* pSDman = G4SDManager::GetSDMpointer();
	const DetectorConstruction* detector =
		(const DetectorConstruction*)
		(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//	detector->GetNumberOfSegmentsInPhantom(fNx,fNy,fNz); //Fill fNx,y,z.
	

	std::vector<G4String> sds = detector->GetUSDNames();


    //------- Energy Deposition will be printed in output folder, for both cases of 100 layers and 8 (1+7)
if(conf()->SphereScorer==1){
	for(size_t i = 1; i< sds.size(); i++) {
		FluenceEnergyDistributionSD* sensDet
	               = (FluenceEnergyDistributionSD*) (pSDman->
			                       FindSensitiveDetector(sds[i]));
		G4String outfile= pathtime + "/" + sds[i] + "ErgfileVoxels.dat";
//		std::fstream outp(outfile.c_str(),std::fstream::out |
//										  std::fstream::in  |
//										  std::fstream::trunc);
                std::fstream outp;
                outp.open(outfile.c_str(),std::fstream::out | std::fstream::in  |   std::fstream::trunc);
                outp << "prova" << G4endl;
                outp << sds.size() << G4endl;
 //  outp << totflux << G4endl;
// outp << G4THitsMap->GetHitsMap(sds[i]) << G4endl;
	      sensDet->FluenceEnergyDistributionSD::DumpAllDetectorCollects(outp);
		outp.close();
    }
} else{
	for(size_t i = 0; i< sds.size(); i++) {
		FluenceEnergyDistributionSD* sensDet
		           = (FluenceEnergyDistributionSD*) (pSDman->
		                           FindSensitiveDetector(sds[i]));
		G4String outfile=pathtime + "/" + sds[i] + "ErgfileVoxels.dat";
//		std::fstream outp(outfile.c_str(),std::fstream::out |
//										  std::fstream::in  |
//										  std::fstream::trunc);
		        std::fstream outp;
				outp.open(outfile.c_str(),std::fstream::out | std::fstream::in  |   std::fstream::trunc);
				outp << "prova" << G4endl;
				outp << sds.size() << G4endl;
 //  outp << totflux << G4endl;
// outp << G4THitsMap->GetHitsMap(sds[i]) << G4endl;
				sensDet->FluenceEnergyDistributionSD::DumpAllDetectorCollects(outp);
		outp.close();
	}
}
}*/




/*

std::vector<G4String> mfdName = fDetector->GetUSDNames();
G4int nMfd = mfdName.size();
for ( G4int idet = 0; idet < nMfd ; idet++)
    {
    G4String detName = mfdName[idet];
    DetectorCollect* DetCol = (DetectorCollect*)(pSDman->FindSensitiveDetector(detName));
//	SDCollect* SDDetCol = (SDCollect*)(pSDman->FindSensitiveDetector(detName));
    std::fstream outp2;
    G4String outfile2="./outputs/" + mfdName[idet] + "multidet.dat";
    outp2.open(outfile2.c_str(), std::fstream::out | std::fstream::in  |   std::fstream::trunc);
    outp2 << "prova2" << G4endl;
    DetCol->SDCollect::DumpCollects(out2);
    //DetectorCollect::DumpDetectorCollect(outp2);
    outp2.close();

 //with DetCol -> i am communicating which detector should be taken in consideration already, so the whole dumpdetectorcollect happen
// for that one only.


}

*/



}


// --
