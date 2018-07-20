#include "config.h"
#include <fstream>      // std::ifstream
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>

using namespace std;
Conf ConfigHandler::conf;

// this file is created in order to pass certain informations in the code and be easy to access and modify

std::map<std::string,std::string> settings;
ConfigHandler &ConfigHandler::getInstance() {
	static ConfigHandler instance = ConfigHandler();
	return instance;
}

void ConfigHandler::closeFile() {
	conf.edistr->close();
	conf.SphereFlux->close();
	conf.fastFlux->close();
	conf.albedoFlux->close();
}

void swap(const std::string& key, std::string& val) {
	auto it = settings.find(key);
	if (it != settings.end()){
		val = it->second;
	}
}

void swap(const std::string& key, bool& val){
	auto it = settings.find(key);
	if (it != settings.end()){
		val = stoll(it->second);
	}
}

void swap(const std::string& key, i64& val){
	auto it = settings.find(key);
	if (it != settings.end()){
		val = stoll(it->second);
	}

}
void swap(const std::string& key, int& val){
	auto it = settings.find(key);
	if (it != settings.end()){
		val = stoi(it->second);
	}

}
void swap(const std::string& key, double& val){
	auto it = settings.find(key);
	if (it != settings.end()){
		val = std::stod(it->second);
	}
}

void loadBinning(std::vector<double>& binning, const std::string& fileName ){

	std::ifstream file (fileName, std::ifstream::in);
	std::string   line;
	while(std::getline(file, line)) 	{
		std::stringstream   linestream(line);
		std::string         val;
		// Read the integers using the operator >>
		linestream >> val;
		binning.push_back( std::stod(val));
	}

}

std::ofstream* openFile(const G4String& name){
	auto file = new std::ofstream();
	file->open(name);
	if(!file->is_open()) {
		std::cout << "impossible to open " << name << endl;
		exit(0);
	}
	return file;
}

ConfigHandler::ConfigHandler() {
// it reads the parameters from an input file that override the default variables
	std::ifstream file ("config.ini", std::ifstream::in);
	std::string   line;
// the first item of the row is the key to read (parameter) and the second is the value
	while(std::getline(file, line)) 	{
		std::stringstream   linestream(line);
		std::string         data;
		std::string         key;
		std::string         val;

		// Read the integers using the operator >>
		linestream >> key >> val;

		settings.insert(std::pair<std::string,std::string>(key,val));
	}

	swap("enableRoom",conf.enableRoom); //true in order to consider the room geometry
	swap("isotropic",conf.isotropic);
	swap("eventNumber",conf.eventNumber);
	swap("accuracy",conf.accuracy);
	swap("SiLayersDep",conf.SiLayersDep); //true for the energy deposition on the Si divided in 40 layers
	swap("DefMaterials",conf.DefMaterials);
	swap("EnableRoom",conf.EnableRoom);
	swap("SphereScorer", conf.SphereScorer); //true for the dummy sphere around the Am-Be source for E Kin
	swap("Sourcexcm",conf.Sourcexcm);
	swap("Sourceycm",conf.Sourceycm);
	swap("Sourcezcm",conf.Sourcezcm);
	swap("DummyScorer",conf.DummyScorer); //true for the dummy Si scorer only for E kin
	swap("print_stored_trajectories",conf.print_stored_trajectories);
	swap("Iondummy",conf.Iondummy); //true if only the E kin of ions in the Si dummy is wanted
	swap("multithreading",conf.multithreading); //if use MT Manager
	swap("numbercores",conf.numbercores); //  number of cores to use for the MT
	//Now load the binning Info
	loadBinning(conf.ebin,"configbinning.ini"); //load the binning energies


	G4String pathtimenow = "./" + conf.timenow; //create a new folder for the ouputs with unique name
	mkdir(pathtimenow, 0777);

	G4String outputEDistr ="./" + conf.timenow + "/Edistr.dat"; //output file for initial E Kin of neutrons
	conf.edistr = openFile(outputEDistr);

	G4String outfiler1="./" + conf.timenow + "/outputSphere.dat"; //output for dummy sphere
	conf.SphereFlux = openFile(outfiler1);

	G4String fastFlux="./" + conf.timenow + "/outputFlux.dat"; //output for dummy fast
	conf.fastFlux = openFile(fastFlux);

	G4String albedoFlux="./" + conf.timenow + "/outputAlbedo.dat"; // output for dummy albedo
	conf.albedoFlux = openFile(albedoFlux);

	G4String albedoDep="./" + conf.timenow + "/outputAlbedoDept.dat"; //output for E deposition in albedo
	conf.albedoDep =openFile(albedoDep);

	G4String fastDep="./" + conf.timenow + "/outputFastDept.dat"; // output for E deposition in fast
	conf.fastDep =openFile(fastDep);
}

const Conf *conf() {
	return &ConfigHandler::getInstance().conf;
}
