#ifndef CONFIG_H
#define CONFIG_H
#include <vector>
#include <chrono>
#include <ctime>
#include <time.h>
typedef unsigned long long u64;
typedef long long i64;

#include <G4String.hh>
struct DetectorName{
	G4String pathSphere ="mySphereScorer/totSphereFlux";
	G4String pathFast = "fastDet/totfastflux";
	G4String pathAlbedo = "albedoDet/totalbedoflux";
	G4String pathFastDep = "fastdep/totaldep";
	G4String pathAlbedoDep = "albedodep/totaldep";
};

static const DetectorName detectorName;


struct Conf{
	bool enableRoom = false;
	i64 eventNumber = 5000;
	bool isotropic = true;
	double accuracy = 1;
	bool SiLayersDep = 0;
	bool  DefMaterials = 1;
	bool EnableRoom = 1;
	bool SphereScorer = 0;
	double Sourcexcm = 0;
	double Sourceycm = 0;
	double Sourcezcm = 0;
	bool DummyScorer = 0;
	std::vector<double> ebin;
	G4String folder;
	bool print_stored_trajectories = false;
	G4String timenow = std::to_string(std::time(0));
	bool Protondummy = false;
};
const Conf* conf();

class ConfigHandler {
public:
	static ConfigHandler& getInstance();
	static Conf conf;
private:
	ConfigHandler();
};



#endif // CONFIG_H
