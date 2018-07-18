// NOT USED BY  MY CODE, TO FULLY REMOVE
#include "G4SystemOfUnits.hh"

#include "DetectorCollect.hh"
#include "SDCollect.hh"
#include <fstream>
#include <iomanip>
#include <map>

//constructing the DetectorCollect class
/*DetectorCollect::DetectorCollect(const G4String& detName,
						 const std::vector<G4double>& eBinning,
						 const std::vector<G4String>& scorParticles)
	: fDetectorName(detName), fBinningVector(eBinning), 
	  fScorParticleVector(scorParticles)
{
	fTotalSDCollect = new SDCollect(eBinning);
	fTotalSDCollect->SetParticleName("all");
// for all the particles, do: the Binning, give a unique name/ID, ..
	for(size_t i = 0; i < scorParticles.size(); i++)
	{
		SDCollect* sdCollect = new SDCollect(eBinning);
		sdCollect->SetParticleName(scorParticles[i]);
		fSDCollectVector.push_back(sdCollect);
	}
}

//=============================================================================

DetectorCollect::~DetectorCollect()
{}

//=============================================================================

//for each particle, check if it matches the wanted one
SDCollect* DetectorCollect::FindParticleSDCollect(
		 	 const G4String& particleName)
{

	for(size_t i = 0 ; i < fSDCollectVector.size() ; i++)
        {
		if(fSDCollectVector[i]->GetParticleName() == particleName)
		{
			return fSDCollectVector[i];

                }

	}
	G4cout << "Cannot find SDCollect for Particle: " << particleName
				 << "\nreturn empty SDCollect" << G4endl;
	SDCollect* emptySDCollect = new SDCollect(fBinningVector);
	return emptySDCollect;
}

*/

