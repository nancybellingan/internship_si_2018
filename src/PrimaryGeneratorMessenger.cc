






#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
//#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
	PrimaryGeneratorAction *PGA) : primGenAction(PGA)
{


	beamParamDir = new G4UIdirectory("/beam/");
	beamParamDir->SetGuidance("set parameters of beam");

	energyDir = new G4UIdirectory("/beam/energy/");
	energyDir->SetGuidance("set energy of beam");

	particlePosDir = new G4UIdirectory("/beam/position/");
	particlePosDir->SetGuidance("set position of particle");

	momentumDir = new G4UIdirectory("/beam/momentum/");
	momentumDir->SetGuidance("set momentum of particle");

	//=========================================================================
	// Messenger commands for Beam Energy

	meanKinEnergy = new G4UIcmdWithADoubleAndUnit("/beam/energy/meanEnergy",this);
	meanKinEnergy->SetGuidance("set mean kinetic energy");
	meanKinEnergy->SetParameterName("Energy",false);
	meanKinEnergy->SetDefaultUnit("MeV");
	meanKinEnergy->SetUnitCandidates("eV keV MeV GeV TeV");
	meanKinEnergy->AvailableForStates(G4State_PreInit, G4State_Idle);

	sigmaEnergy = new G4UIcmdWithADouble("/beam/energy/EnergyFWHM",this);
	sigmaEnergy->SetGuidance("set energy FWHM");
	sigmaEnergy->SetParameterName("percent",false);
//	sigmaEnergy->SetDefaultUnit("keV");
//	sigmaEnergy->SetUnitCandidates("eV keV MeV GeV TeV");
	sigmaEnergy->AvailableForStates(G4State_PreInit, G4State_Idle);

	//=========================================================================
	// Messenger commands for Particle Position

	posX = new G4UIcmdWithADoubleAndUnit("/beam/position/posX",this);
	posX->SetGuidance("set X coordinate of particle");
	posX->SetParameterName("position",false);
	posX->SetDefaultUnit("mm");
	posX->SetUnitCandidates("um mm cm m");
	posX->AvailableForStates(G4State_PreInit, G4State_Idle);

	posY = new G4UIcmdWithADoubleAndUnit("/beam/position/posY",this);
	posY->SetGuidance("set Y coordinate of particle");
	posY->SetParameterName("position",false);
	posY->SetDefaultUnit("mm");
	posY->SetUnitCandidates("um mm cm m");
	posY->AvailableForStates(G4State_PreInit, G4State_Idle);

	posZ = new G4UIcmdWithADoubleAndUnit("/beam/position/posZ",this);
	posZ->SetGuidance("set Z coordinate of particle");
	posZ->SetParameterName("position",false);
	posZ->SetDefaultUnit("mm");
	posZ->SetUnitCandidates("um mm cm m");
	posZ->AvailableForStates(G4State_PreInit, G4State_Idle);

	sigmaPosX = new G4UIcmdWithADoubleAndUnit("/beam/position/posX/sigmaPosX",this);
	sigmaPosX->SetGuidance("set sigma of X coordinate");
	sigmaPosX->SetParameterName("position",false);
	sigmaPosX->SetDefaultUnit("mm");
	sigmaPosX->SetUnitCandidates("mm cm m");
	sigmaPosX->AvailableForStates(G4State_PreInit, G4State_Idle);

	sigmaPosY = new G4UIcmdWithADoubleAndUnit("/beam/position/posY/sigmaPosY",this);
	sigmaPosY->SetGuidance("set sigma of Y coordinate");
	sigmaPosY->SetParameterName("position",false);
	sigmaPosY->SetDefaultUnit("mm");
	sigmaPosY->SetUnitCandidates("mm cm m");
	sigmaPosY->AvailableForStates(G4State_PreInit, G4State_Idle);

	sigmaPosZ = new G4UIcmdWithADoubleAndUnit("/beam/position/posZ/sigmaPosZ",this);
	sigmaPosZ->SetGuidance("set sigma of Z coordinate");
	sigmaPosZ->SetParameterName("position",false);
	sigmaPosZ->SetDefaultUnit("mm");
	sigmaPosZ->SetUnitCandidates("mm cm m");
	sigmaPosZ->AvailableForStates(G4State_PreInit, G4State_Idle);

	//=========================================================================
	// Messenger commands for particle momentum
/*
	momX = new G4UIcmdWithADouble("/beam/momentum/momX",this);
	momX->SetGuidance("set X momentum of particle");
	momX->SetParameterName("position",false);
	momX->AvailableForStates(G4State_PreInit, G4State_Idle);

	momY = new G4UIcmdWithADouble("/beam/momentum/momY",this);
	momY->SetGuidance("set Y momentum of particle");
	momY->SetParameterName("position",false);
	momY->AvailableForStates(G4State_PreInit, G4State_Idle);

	momZ = new G4UIcmdWithADouble("/beam/momentum/momZ",this);
	momZ->SetGuidance("set Z momentum of particle");
	momZ->SetParameterName("position",false);
	momZ->AvailableForStates(G4State_PreInit, G4State_Idle);

	sigmaMomX = new G4UIcmdWithADouble("/beam/momentum/sigmaMomX",this);
	sigmaMomX->SetGuidance("set momentum in X direction");
	sigmaMomX->SetParameterName("momentum",false);
	sigmaMomX->AvailableForStates(G4State_PreInit,G4State_Idle);

	sigmaMomY = new G4UIcmdWithADouble("/beam/momentum/sigmaMomY",this);
	sigmaMomY->SetGuidance("set momentum in Y direction");
	sigmaMomY->SetParameterName("momentum",false);
	sigmaMomY->AvailableForStates(G4State_PreInit,G4State_Idle);

	sigmaMomZ = new G4UIcmdWithADouble("/beam/momentum/sigmaMomZ",this);
	sigmaMomZ->SetGuidance("set momentum in Z direction");
	sigmaMomZ->SetParameterName("momentum",false);
	sigmaMomZ->AvailableForStates(G4State_PreInit,G4State_Idle);
*/
}


PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
	delete beamParamDir;
	delete energyDir;
	delete meanKinEnergy;
	delete sigmaEnergy;
/*	delete particlePosDir;
	delete posX;
	delete posY;
	delete posZ;
	delete sigmaPosX;
	delete sigmaPosY;
	delete sigmaPosZ;
	delete momX;
	delete momY;
	delete momZ; 
	delete sigmaMomX;
	delete sigmaMomY;
	delete sigmaMomZ;
*/
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand *command,
					G4String newValue)
{
	newValue.clear();
	if(command == meanKinEnergy)
	{	
        //	primGenAction->setGunEnergy(meanKinEnergy->GetNewDoubleValue(newValue));
        //	G4cout << "Mean Gun Energy set to " <<  newValue << " MeV" << G4endl;
	}
	else if(command == sigmaEnergy)
	{
        //	primGenAction->setFWHMGunEnergy(sigmaEnergy->GetNewDoubleValue(newValue));
        //	G4cout << "Beam Energy Distribution Width set to " << newValue
        //			<< G4endl;
	}
/*	else if(command == posX)
	{
		primGenAction->setXPos(posX->GetNewDoubleValue(newValue));
	}
	else if(command == posY)
	{
		primGenAction->setYPos(posY->GetNewDoubleValue(newValue));
	}
	else if(command == posZ)
	{
		primGenAction->setZPos(posZ->GetNewDoubleValue(newValue));
	}
	else if(command == sigmaPosX)
	{
		primGenAction->setSigmaXPos(sigmaPosX->GetNewDoubleValue(newValue));
		G4cout << "Beam X Position Distribution Width set to " 
				<< newValue << " keV" << G4endl;
	}
	else if(command == sigmaPosY)
	{
		primGenAction->setSigmaYPos(sigmaPosX->GetNewDoubleValue(newValue));
		G4cout << "Beam Y Position Distribution Width set to " 
				<< newValue << " keV" << G4endl;
	}
	else if(command == sigmaPosZ)
	{
		primGenAction->setSigmaZPos(sigmaPosX->GetNewDoubleValue(newValue));
		G4cout << "Beam Z Position Distribution Width set to " 
				<< newValue << " keV" << G4endl;
	}
	else if(command == sigmaMomX)
	{
		primGenAction->setSigmaXMom(sigmaMomX->GetNewDoubleValue(newValue));
		G4cout << "Beam X Momentum Distribution Width set to " 
				<< newValue << G4endl;
	}
	else if(command == sigmaMomY)
	{
		primGenAction->setSigmaYMom(sigmaMomY->GetNewDoubleValue(newValue));
		G4cout << "Beam X Momentum Distribution Width set to " 
				<< newValue << G4endl;
	} */
}
