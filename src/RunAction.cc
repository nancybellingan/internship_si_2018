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
#include <mutex>
#include <iostream>
#include "config.h"
#include <string>
//=======================================================================
// RunAction
//  
//
//
//=======================================================================
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Constructor
RunAction::RunAction()
  : G4UserRunAction(),
    fNx(0), fNy(0), fNz(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Destructor.
RunAction::~RunAction()
{
    // if the debug is activated from main, it is given the information that the destructor is happening
	G4cout << "Destructor RunAction" << G4endl;
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

  auto flux1 = re02Run->GetSphereFlux();

  if (conf()->SphereScorer==1 && conf()->totdata == 1){
//static std::mutex mutexFileWrite4;
	  std::fstream outp1;
//	  G4int sum = 0;
	  G4String outfiler1=conf()->timenow +"/outputspheretot.dat";
	  outp1.open(outfiler1.c_str(),std::fstream::out | std::fstream::in  |   std::fstream::trunc);
	  std::vector<G4THitsMap<G4double>*> ref = re02Run->GetSphereFlux();
//	  auto flux1 = re02Run->GetSphereFlux();
//	  std::lock_guard<std::mutex> lock(mutexFileWrite4);
	  for (uint i=0;i<ref.size();i++){
		  auto energy = conf()->ebin[i];
		  for(auto line : *(ref[i]) ){

			  outp1 << "bin n° " << i <<" energy: " << energy << " key " << line.first << " boh " << *line.second << G4endl;
//			  sum += *line.second;
		  }

		  //gg all
	  }
//	  outp1 << "tot events" << sum << G4endl;
	  outp1.close();
  }

  if (conf()->DummyScorer==1 && conf()->totdata == 1){
	  std::fstream outp1;
	  G4String outfiler2=conf()->timenow +"/outputfastdummytot.dat";
	  G4String outfiler3=conf()->timenow +"/outputalbedodummytot.dat";
	  outp1.open(outfiler2.c_str(),std::fstream::out | std::fstream::in  |   std::fstream::trunc);
	  outp1 << "prova" << G4endl;
	  auto flux2 = re02Run->GetFastFlux();
	  for (uint i=0;i<flux2.size();i++){
		  auto energy = conf()->ebin[i];
		  for(auto line : *(flux2[i]) ){
			  if (line.first == 0) {
			  outp1 << "bin n° " << i <<" energy: " << energy << "number events in the bin \t" << *line.second << G4endl;
		  }
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
			  if (line.first == 0) {
			  outp1 << "bin n° " << i <<" energy: " << energy << "number events in the bin \t" << *line.second << G4endl;
		  }
		  }
		  //gg all
	  }
	  outp1.close();
  }


}



