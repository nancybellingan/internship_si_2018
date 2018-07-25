/*

#ifndef DETECTORCOLLECT_HH_
#define DETECTORCOLLECT_HH_ 1

#include "SDCollect.hh"
#include <vector>
#include <fstream>

class map;

class DetectorCollect
{
	public:
		DetectorCollect(const G4String& detName,
										const std::vector<G4double>& eBinning,
										const std::vector<G4String>& scorParticles);
		~DetectorCollect();

		void SetDetectorName(const G4String& detName) { fDetectorName = detName; };
		G4String GetDetectorName() { return fDetectorName; };

		void SetBinningVector(const std::vector<G4double>& binningVector)
													{ fBinningVector = binningVector; };
		std::vector<G4double> GetBinningVector() { return fBinningVector; };

		void SetScorParticleVector(const std::vector<G4String>& particleVector)
													{ fScorParticleVector = particleVector; };
		std::vector<G4String> GetScorParticleVector() { return fScorParticleVector; };

		void SetSDCollectVector(const std::vector<SDCollect*> sdCollectVector)
													{ fSDCollectVector = sdCollectVector; };
		std::vector<SDCollect*> GetSDCollectVector() { return fSDCollectVector; };
		SDCollect* FindParticleSDCollect(const G4String& particleName);

		        void SetTotalSDCollect( SDCollect & sdCollect){
					fTotalSDCollect = &sdCollect;
				};
                SDCollect* GetTotalSDCollect() { return fTotalSDCollect; };

		void DumpDetectorCollect(std::fstream& out);

	private:

		G4String fDetectorName;
		std::vector<G4double> fBinningVector;
		std::vector<G4String> fScorParticleVector;
		std::vector<SDCollect*> fSDCollectVector;
                SDCollect* fTotalSDCollect;

};

#endif // end of class DetectorCollect
*/
