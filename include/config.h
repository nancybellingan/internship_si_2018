#ifndef CONFIG_H
#define CONFIG_H
#include <vector>
#include <chrono>
#include <ctime>
#include <time.h>
#include <vector>
typedef unsigned long long u64;
typedef long long i64;

#include <G4String.hh>
struct DetectorName{
	G4String pathSphere ="mySphereScorer/totSphereFlux";
	G4String pathFast = "fastDet/totfastflux";
	G4String pathAlbedo = "albedoDet/totalbedoflux";
	G4String pathFastDep = "fastDet/totaldep";
	G4String pathAlbedoDep = "albedoDet/totaldep";
	G4String pathAlbedoDepScorer ="EDepFast";
};
static const std::vector<double> correctionfactor{ // for the silicon layers correction factor
	1, // from 0 to 10 um
	1, // from 10 to 20 um
	1, // from 20 to 30 um
	1, // from 30 to 40 um
	0.94, // from 40 to 50 um
	0.7059, // from 50 to 60 um
	0.4052, // from 60 to 70 um
	0.1831, // from 70 to 80 um
	0.0737, // from 80 to 90 um
	0.02458 // from 90 to 100 um
};

static const DetectorName detectorName;


struct Conf{
	bool DefMaterials = 1;
	bool SiLayersDep = 0;
	bool EnableRoom = 1;
	bool SphereScorer = 0;
	double Sourcexcm = 0;
	double Sourceycm = 0;
	double Sourcezcm = 0;
	double distancephantsurf = 100;
	bool DummyScorer = 0;
	std::vector<double> ebin;
	G4String folder;
	bool print_stored_trajectories = false;
	G4String timenow = std::to_string(std::time(0));
	bool Iondummy = false;
	mutable std::ofstream* edistr = nullptr;
	mutable std::ofstream* SphereFlux = nullptr;
	mutable std::ofstream* fastFlux = nullptr;
	mutable std::ofstream* albedoFlux = nullptr;
	mutable std::ofstream* albedoDep = nullptr;
	mutable std::ofstream* fastDep = nullptr;
	mutable std::ofstream* albedoTotDep = nullptr;
	mutable std::ofstream* fastTotDep = nullptr;
	mutable std::ofstream* albedoTotDepfilt = nullptr;
	mutable std::ofstream* fastTotDepfilt = nullptr;
	mutable std::ofstream* phantomFlux = nullptr;
	mutable std::ofstream* phantomfrontal = nullptr;
	mutable std::ofstream* backalbedo = nullptr;
	bool multithreading = false;
	int numbercores = 0;
	bool albedocentre = true;
	bool totdata = false;
	bool phantomon = true;
	bool EnableRoomv2 = false;
	bool phantomscorer = true;
	bool lightsim = true;
	int sensorposz = 0;
	bool backflux = 0;
	bool albedoon = 1;
	bool faston = 1;
};
const Conf* conf();

class ConfigHandler {
public:
	static ConfigHandler& getInstance();
	static Conf conf;
	static void closeFile();
private:
	ConfigHandler();
};



#endif // CONFIG_H
