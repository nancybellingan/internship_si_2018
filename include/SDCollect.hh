
#ifndef SDCOLLECT_HH_
#define SDCOLLECT_HH_ 1

#include "globals.hh"
#include <map>
#include <vector>
#include <fstream>

class SDCollect
{
	public:
		
		SDCollect(const std::vector<G4double>& thisBinning);
		~SDCollect();
		
		void SetParticleName(const G4String& particleName) 
							{ fParticleName = particleName; };
		G4String GetParticleName() { return fParticleName; };	
		void SetFlux(const G4double& flux) 
							{ fFlux = flux; };
		void AddToFlux(const G4double& flux) 
							{ fFlux += flux; };
		G4double GetFlux() { return fFlux; };

		void SetkEnergy (const G4double& kenergy)
							{ fkenergy = kenergy; };

		G4double GetkEnergy() { return fkenergy; };
		

		void SetSpectralFluxMap(const std::map<G4double,G4double>& thisMap)
							{ fSpectralFluxMap = thisMap; };
		std::map<G4double,G4double> GetSpectralFluxMap() const
							{ return fSpectralFluxMap; };
		void AddToSpectralFluxMap(const G4double& energy, const G4double& flux);
		void SumSpectralFluxMapVal(const G4double& key, const G4double& flux)
							{ fSpectralFluxMap[key] += flux; };
		
		void PrintSpectralFluxMap() const;
		

		void SetBinningVector(const std::vector<G4double>& thisVector)
							{ fBinningVector = thisVector; };
		std::vector<G4double> GetBinningVector() const 
							{ return fBinningVector; };

		SDCollect operator + ( SDCollect b );
		SDCollect operator += ( SDCollect b );

	void DumpCollects(std::fstream& out);
	private:
		void InitSpectralFluxMap();


	private:
		G4String fParticleName;
		G4double fFlux;
		G4double fDoseDeposit;
		G4double fkenergy;
		std::map<G4double,G4double> fSpectralFluxMap;
		std::vector<G4double> fBinningVector;
	
	
};

#endif
// End of class SDCollect
