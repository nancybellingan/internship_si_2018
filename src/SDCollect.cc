

#include <vector>
#include <map>
#include "SDCollect.hh"
#include "G4SystemOfUnits.hh"

SDCollect::SDCollect(const std::vector<G4double>& thisBinning)
	: fParticleName("any"), fFlux(0.), fkenergy(0.),
		fBinningVector(thisBinning)
{
	InitSpectralFluxMap();
}

/*SDCollect::SDCollect()
	: fParticleName("any"), fFlux(0.), fDoseDeposit(0.)
{	
}
*/

SDCollect::~SDCollect()
{}
// goes through each binning / energy step and add the flux in the right energy range for the whole mapped Energies
void SDCollect::AddToSpectralFluxMap(const G4double& energy, const G4double& flux)
{
	G4double lowerBinEdge = 0;
	std::map<G4double,G4double>::iterator ithBin = fSpectralFluxMap.begin();
	if(energy < ithBin->first) return;
	
	for(; ithBin != fSpectralFluxMap.end(); ++ithBin)
	{
                if(energy >= lowerBinEdge && energy < ithBin->first)
		{
			fSpectralFluxMap[lowerBinEdge] += flux;
		}
		lowerBinEdge = ithBin->first;
	}	
}

void SDCollect::InitSpectralFluxMap()
{
	std::vector<G4double>::iterator ithBin = fBinningVector.begin();
	for(; ithBin != fBinningVector.end();++ithBin)
	{
		fSpectralFluxMap[*ithBin];
	}
}
//it gives back B, which includes particle name, flux map, dose and energy deposition, etc
SDCollect SDCollect::operator + (SDCollect b)
{
        (*this).fParticleName = b.fParticleName; //the particle name is the same, but the flux, dose deposition, energy deposition are added to previous
	(*this).fFlux = (*this).fFlux + b.fFlux;
	(*this).fkenergy = (*this).fkenergy;

//	if(this.fSpectralFluxMap.size() == b.fSpectralFluxMap.size())
	if((*this).fBinningVector == b.fBinningVector)
	{
		for(size_t i=0; i < b.fBinningVector.size(); ++i)
		{
                        (*this).fSpectralFluxMap[i] = (*this).fSpectralFluxMap[i]+ b.fSpectralFluxMap[i];
		}
	}
	else {
		G4cerr << "ERROR: Cannot add to SpectralMaps!" << G4endl;
		G4cerr << "Assigning Map2 to Map1" << G4endl;
	}

	(*this).fBinningVector = b.fBinningVector;
	return *this;
}


SDCollect SDCollect::operator += (SDCollect b)
{
	*this = *this + b;
	return *this;
}

void SDCollect::DumpCollects(
		std::fstream& out)
{
	auto it = fSpectralFluxMap.begin();
	auto itEnd = fSpectralFluxMap.end();

	for(; it != itEnd; ++it)
	{
out << "flux, energy " << it->first << it->second << G4endl;
// G4cout << "flux, energy " << SDCollect::getFlux() << SD:Collect::getkEnergy() << g4endl;
	}
}

