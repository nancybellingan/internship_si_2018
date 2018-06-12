



#ifndef PRIMARYGENERATORACTIONMESSENGER_HH
#define PRIMARYGENERATORACTIONMESSENGER_HH 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
//#include " G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"


class PrimaryGeneratorMessenger : public G4UImessenger
{
  public:
	PrimaryGeneratorMessenger(PrimaryGeneratorAction *PGA);
	~PrimaryGeneratorMessenger();

	void SetNewValue(G4UIcommand* command, G4String newValue);

  private:

	PrimaryGeneratorAction *primGenAction;
	G4UIdirectory *beamParamDir;
	G4UIdirectory *energyDir;
	G4UIdirectory *particlePosDir;
	G4UIdirectory *momentumDir;

	G4UIcmdWithADoubleAndUnit *meanKinEnergy;
	G4UIcmdWithADouble *sigmaEnergy;

	G4UIcmdWithADoubleAndUnit *posX;
	G4UIcmdWithADoubleAndUnit *posY;
	G4UIcmdWithADoubleAndUnit *posZ;
	G4UIcmdWithADoubleAndUnit *sigmaPosX;
	G4UIcmdWithADoubleAndUnit *sigmaPosY;
	G4UIcmdWithADoubleAndUnit *sigmaPosZ;


	G4UIcmdWithADouble *momX;
	G4UIcmdWithADouble *momY;
	G4UIcmdWithADouble *momZ;

	G4UIcmdWithADouble *sigmaMomX;
	G4UIcmdWithADouble *sigmaMomY;
	G4UIcmdWithADouble *sigmaMomZ;

};


#endif // end class PrimaryGeneratorMessenger
