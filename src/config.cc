#include "config.h"
#include <fstream>      // std::ifstream
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace std;
Conf ConfigHandler::conf;



std::map<std::string,std::string> settings;
ConfigHandler &ConfigHandler::getInstance() {
	static ConfigHandler instance = ConfigHandler();
	return instance;
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

ConfigHandler::ConfigHandler() {

	std::ifstream file ("config.ini", std::ifstream::in);
	std::string   line;

	while(std::getline(file, line)) 	{
		std::stringstream   linestream(line);
		std::string         data;
		std::string         key;
		std::string         val;


		// Read the integers using the operator >>
		linestream >> key >> val;

		settings.insert(std::pair<std::string,std::string>(key,val));
	}

	swap("enableRoom",conf.enableRoom);
	swap("isotropic",conf.isotropic);
	swap("eventNumber",conf.eventNumber);
	swap("accuracy",conf.accuracy);
	swap("SiLayersDep",conf.SiLayersDep);
	swap("DefMaterials",conf.DefMaterials);
	swap("EnableRoom",conf.EnableRoom);
	swap("SphereScorer", conf.SphereScorer);
	swap("Sourcexcm",conf.Sourcexcm);
	swap("Sourceycm",conf.Sourceycm);
	swap("Sourcezcm",conf.Sourcezcm);
	swap("DummyScorer",conf.DummyScorer);
	//forgive me nancy
	swap("print_stored_trajectories",conf.print_stored_trajectories);
	swap("Protondummy",conf.Protondummy);
	//Now load the binning Info
	loadBinning(conf.ebin,"configbinning.ini");

}

const Conf *conf() {
	return &ConfigHandler::getInstance().conf;
}
