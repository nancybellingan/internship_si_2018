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
// $Id: ThermalNeutronScattering.cc 71037 2013-06-10 09:20:54Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   ThermalNeutronScattering
//
// Author: 3 June 2010 V. Ivanchenko
//
// Modified:
// 03.06.2011 V.Ivanchenko change design - now first default constructor 
//            is called, HP model and cross section are added on top
//
//----------------------------------------------------------------------------
//
// HP model for n with E < 20 MeV

#include "ThermalNeutronScattering.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4Neutron.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElastic.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScattering.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include "config.h"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(ThermalNeutronScattering);

G4ThreadLocal G4bool ThermalNeutronScattering::wasActivated = false;
G4ThreadLocal G4HadronElasticPhysics* ThermalNeutronScattering::mainElasticBuilder = 0;


//=========================================================================
ThermalNeutronScattering::ThermalNeutronScattering(G4int ver)
  : G4VPhysicsConstructor("hElasticWEL_CHIPS_HP"), verbose(ver)
{
  if(verbose > 1) { 
    G4cout << "### ThermalNeutronScattering: " << GetPhysicsName() 
	   << G4endl; 
  }
  mainElasticBuilder = new G4HadronElasticPhysics(verbose);
}

//=========================================================================
ThermalNeutronScattering::~ThermalNeutronScattering()
{
  delete mainElasticBuilder;
}

//=========================================================================
void ThermalNeutronScattering::ConstructParticle()
{
  // G4cout << "G4HadronElasticPhysics::ConstructParticle" << G4endl;
  mainElasticBuilder->ConstructParticle();
}

//========================================================================
void ThermalNeutronScattering::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;
  //Needed because this is a TLS object and this method is called by all threads
  if ( ! mainElasticBuilder ) mainElasticBuilder = new G4HadronElasticPhysics(verbose);
  mainElasticBuilder->ConstructProcess();

  mainElasticBuilder->GetNeutronModel()->SetMinEnergy(19.5*MeV);

  G4HadronicProcess* hel = mainElasticBuilder->GetNeutronProcess();
  G4NeutronHPElastic* hp = new G4NeutronHPElastic();
	hp->SetMinEnergy( 4.0*eV);
  hel->RegisterMe(hp);
  hel->AddDataSet(new G4NeutronHPElasticData());

	G4NeutronHPThermalScattering* thermal = new G4NeutronHPThermalScattering();
	thermal->SetMaxEnergy( 4.0*eV);
	hel->RegisterMe(thermal);
	hel->AddDataSet(new G4NeutronHPThermalScatteringData() );

  if(verbose > 1) {
    G4cout << "### G4HadronElasticPhysicsHP + ThermalNeutronScattering is constructed " 
	   << G4endl;
  }
}

void ThermalNeutronScattering::SetCuts(){



	if (conf()->lightsim == 1)
	  {
	  if (conf()->EnableRoomv2 ==1 || conf()->EnableRoom ==1){
   G4Region* detectorconcrete2 = G4RegionStore::GetInstance()->GetRegion("concretewalls");
    G4ProductionCuts* detectorCutsconcrete2 = new G4ProductionCuts();
	 detectorCutsconcrete2->SetProductionCut(50*CLHEP::mm,"neutron");
	detectorconcrete2->SetProductionCuts(detectorCutsconcrete2);

   }
  }
	if(conf()->phantomon == 1){
	 G4Region* phantomarea1 = G4RegionStore::GetInstance()->GetRegion("phantomregion");
   G4ProductionCuts* detectorCutsphantom1 = new G4ProductionCuts();
    detectorCutsphantom1->SetProductionCut(0.1*CLHEP::mm,"neutron");
	phantomarea1->SetProductionCuts(detectorCutsphantom1);
	 }


   G4Region* detector1 = G4RegionStore::GetInstance()->GetRegion("detectorRegion");
   G4ProductionCuts* detectorCuts1 = new G4ProductionCuts();
   detectorCuts1->SetProductionCut(0.001*CLHEP::mm,"neutron");
   detector1->SetProductionCuts(detectorCuts1);

}




