// src/PiotPhaseSpace.cc

#include "PiotPhaseSpace.hh"

PiotPhaseSpace::PiotPhaseSpace()
: particleName("spaceFile.spec"), atomNum(0), nuclNum(0),
  xPos(0.0), yPos(0.0), zPos(0.0),
  xMom(0.0), yMom(0.0), zMom(0.0),
  kinEnergy(0.0)
{;}

//----------------------------------------------------------------------------
// With particle name
PiotPhaseSpace::PiotPhaseSpace(
                    const std::string &partName,   	// PDG ID
                    const unsigned int &atomicNum,   	// atomic number
                    const unsigned int &nucleonNum,  	// nucleon number
                    const double &xPosition,         	// X position 
                    const double &yPosition,         	// Y position
                    const double &zPosition,         	// Z position
                    const double &xMomentum,         	// x-momentum
                    const double &yMomentum,         	// y-momentum
                    const double &zMomentum,         	// z-momentum
                    const double &kineticEnergy)      // kinetic energy
: particleName(partName), 
  atomNum(atomicNum), 
  nuclNum(nucleonNum),
  xPos(xPosition), 
  yPos(yPosition), 
  zPos(zPosition),
  xMom(xMomentum), 
  yMom(yMomentum), 
  zMom(zMomentum),
  kinEnergy(kineticEnergy)
{;}

//----------------------------------------------------------------------------
// With particle name
PiotPhaseSpace::PiotPhaseSpace(const PiotPhaseSpace &pps)
: particleName(pps.particleName), atomNum(pps.atomNum), 
  nuclNum(pps.nuclNum), xPos(pps.xPos), yPos(pps.yPos), zPos(pps.zPos),
  xMom(pps.xMom), yMom(pps.yMom), zMom(pps.zMom),
  kinEnergy(pps.kinEnergy)
{;}

//----------------------------------------------------------------------------

PiotPhaseSpace::~PiotPhaseSpace()
{;}

