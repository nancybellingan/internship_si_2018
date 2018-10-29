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
// $Id: PhysicsList.icc 67969 2013-03-13 09:44:42Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   PhysicsList
//
// Author: 2006 G.Folger
//
// based on QGSP_BIC
// Modified:
// 26.04.2007 G.Folger: Enable quasielastic for QGS string model
// 16.05.2007 V.Ivanchenko: rename EM builders
// 04.06.2010 G.Folger: Use new ctor for builders
// 16.08.2010 H.Kurashige: Remove inclusion of G4ParticleWithCuts 
// 16.10.2012 A.Ribon: Use new default stopping and ion physics
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
//#include "G4EmStandardPhysics.hh"
// G4EmStandardPhysics replaced with
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
//#include "G4HadronElasticPhysicsHP.hh"
// G4HadronElasticPhysiscsHP extended to 
#include "ThermalNeutronScattering.hh"

#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
//#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "config.h"
#include "PhysicsList.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4UserSpecialCuts.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4ProcessManager.hh"

template<class T> TPhysicsList<T>::TPhysicsList(G4int ver):  T()
{
  // default cut value  (1.0mm) 
  // defaultCutValue = 1.0*CLHEP::mm;

  G4DataQuestionaire it(photon, neutron);
  G4cout << "<<< Geant4 Physics List simulation engine: PhysicsList 2.0"<<G4endl;
  G4cout <<G4endl;



  
  
  if (conf()->lightsim == 1){
this->defaultCutValue = 200*CLHEP::mm;
} else {
this->defaultCutValue =  1.0*CLHEP::mm;
}
//    
    
  this->SetVerboseLevel(ver);

  // EM Physics
  this->RegisterPhysics( new G4EmStandardPhysics_option3(ver) );

  // Synchroton Radiation & GN Physics
  this->RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays
  this->RegisterPhysics( new G4DecayPhysics(ver) );

   // Hadron Elastic scattering
  this->RegisterPhysics( new ThermalNeutronScattering(ver) );

   // Hadron Physics
  this->RegisterPhysics(  new G4HadronPhysicsQGSP_BIC_HP(ver));

  // Stopping Physics
  this->RegisterPhysics( new G4StoppingPhysics(ver) );

  // Ion Physics
  this->RegisterPhysics( new G4IonPhysics(ver));
  

}

template<class T> TPhysicsList<T>::~TPhysicsList()
{
}

template<class T> void TPhysicsList<T>::SetCuts()
{





  if (this->verboseLevel >1){
    G4cout << "PhysicsList::SetCuts:";
  }  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 

  this->SetCutsWithDefault();   

	//Set proton cut value to 0 for producing low energy recoil nucleus
	this->SetCutValue(0, "proton");

 // if (this->verboseLevel >0)
 G4VUserPhysicsList::DumpCutValuesTable();  




  if (conf()->lightsim == 1)
	{
	if (conf()->EnableRoomv2 ==1 || conf()->EnableRoom ==1){
 G4Region* detectorconcrete1 = G4RegionStore::GetInstance()->GetRegion("concretewalls");
  G4ProductionCuts* detectorCutsconcrete = new G4ProductionCuts();
 G4double prodcutconcrete = 50*CLHEP::mm;
 detectorCutsconcrete->SetProductionCut(prodcutconcrete);
   detectorCutsconcrete->SetProductionCut(50*CLHEP::mm,"neutron");
  detectorconcrete1->SetProductionCuts(detectorCutsconcrete);

 } 
} 

  
 // G4ProcessManager* pmanager = G4Gamma::Gamma()->GetProcessManager();
// pmanager->AddProcess(new G4UserSpecialCuts(),-1,-1,1);
 
 
 
	if(conf()->phantomon == 1){
	 G4Region* phantomarea = G4RegionStore::GetInstance()->GetRegion("phantomregion");
   G4ProductionCuts* detectorCutsphantom = new G4ProductionCuts();
  G4double prodcutphantom = 0.1*CLHEP::mm;
   detectorCutsphantom->SetProductionCut(prodcutphantom);
 detectorCutsphantom->SetProductionCut(100000000,"gamma");
    detectorCutsphantom->SetProductionCut(0.1*CLHEP::mm,"neutron");
	phantomarea->SetProductionCuts(detectorCutsphantom);
     } 
     
  G4double prodcut = 0.001*CLHEP::mm;   
   G4Region* detector = G4RegionStore::GetInstance()->GetRegion("detectorRegion");
   G4ProductionCuts* detectorCuts = new G4ProductionCuts();
   detectorCuts->SetProductionCut(prodcut);
   detectorCuts->SetProductionCut(0,"proton");
   detectorCuts->SetProductionCut(100000000,"gamma");
   detectorCuts->SetProductionCut(0.001*CLHEP::mm,"neutron");
   detector->SetProductionCuts(detectorCuts);

}



// 2002 by J.P. Wellisch
