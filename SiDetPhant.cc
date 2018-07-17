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
/// \file runAndEvent/RE02/RE02.cc
/// \brief Main program of the runAndEvent/RE02 example
//
//
// $Id: RE02.cc 76292 2013-11-08 13:10:20Z gcosmo $
//
// 

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"


#include "G4MTRunManager.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"    

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "config.h"
int main(int argc,char** argv) 
{
	//thi is to STOP on start, wait a USER input, so we an attach a debugger.
	printf("Type any char (and enter) to continue \n");
	//int x;
	//scanf("%d",&x);
	// Construct the default run manager


	if (conf()->multithreading ==0){
		G4RunManager * runManager = new G4RunManager;
		DetectorConstruction* detector = new DetectorConstruction;
		// set mandatory initialization classes
printf("using single thread");
        runManager->SetUserInitialization(detector);
		runManager->SetUserInitialization(new PhysicsList());

#ifdef G4VIS_USE
		G4VisManager* visManager = new G4VisExecutive;
		visManager->Initialize();
#endif

		runManager->SetUserInitialization(new ActionInitialization);
		runManager->Initialize();
		// get the pointer to the User Interface manager
		G4UImanager * pUI = G4UImanager::GetUIpointer();

		if(argc==1)
		{
#ifdef G4UI_USE // interactive mode : define UI session
			G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
			pUI->ApplyCommand("/control/execute macros/vis.mac");
#endif
			ui->SessionStart();
			delete ui;
#endif
		}
		else // batch mode
		{
			G4String command = "/control/execute ";
			G4String fileName = argv[1];
			pUI->ApplyCommand(command+fileName);
		}

#ifdef G4VIS_USE
		delete visManager;
#endif
		delete runManager;
	} else if (conf()->multithreading==1) {
		printf("using multithreading with \n");
		printf("%d",conf()->numbercores);

		G4MTRunManager * runManager = new G4MTRunManager;
		runManager->SetNumberOfThreads(conf()->numbercores);
		DetectorConstruction* detector = new DetectorConstruction;
		// set mandatory initialization classes

		runManager->SetUserInitialization(detector);
		runManager->SetUserInitialization(new PhysicsList());

#ifdef G4VIS_USE
		G4VisManager* visManager = new G4VisExecutive;
		visManager->Initialize();
#endif

		runManager->SetUserInitialization(new ActionInitialization);
		runManager->Initialize();
		// get the pointer to the User Interface manager
		G4UImanager * pUI = G4UImanager::GetUIpointer();

		if(argc==1)
		{
#ifdef G4UI_USE // interactive mode : define UI session
			G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
			pUI->ApplyCommand("/control/execute macros/vis.mac");
#endif
			ui->SessionStart();
			delete ui;
#endif
		}
		else // batch mode
		{
			G4String command = "/control/execute ";
			G4String fileName = argv[1];
			pUI->ApplyCommand(command+fileName);
		}

#ifdef G4VIS_USE
		delete visManager;
#endif
		delete runManager;

	}




	ConfigHandler::closeFile();

	return 0;
}
