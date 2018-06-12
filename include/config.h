#ifndef CONFIG_H
#define CONFIG_H
#include <vector>
typedef unsigned long long u64;
typedef long long i64;

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
	std::vector<double> ebin;
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
