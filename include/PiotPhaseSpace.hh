// include/PiotPhaseSpace.hh



#ifndef PIOTPHASESPACE_HH
#define PIOTPHASESPACE_HH 1

#include <iostream>
#include <string>

class PiotPhaseSpace
{
  public:
	PiotPhaseSpace();	// default constructor initializing everthing with 0

	PiotPhaseSpace(const PiotPhaseSpace &pps); // copy constructor

	// With particle name
	PiotPhaseSpace( const std::string &partName,			// particle name
					const unsigned int &atomicNum,	// atomic number
					const unsigned int &nucleonNum,	// nucleon number
					const double &xPosition,		// X position 
					const double &yPosition,		// Y position
					const double &zPosition,		// Z position
					const double &xMomentum,		// x-momentum
					const double &yMomentum,		// y-momentum
					const double &zMomentum,		// z-momentum
					const double &kineticEnergy);		// kinetic energy

	virtual ~PiotPhaseSpace();

	//------------------------------------------------------------------------
	// methods to get private variables

	std::string getParticleName() const { return particleName; };
	unsigned int getAtomNum() const { return atomNum; };
	unsigned int getNuclNum() const { return nuclNum; };

	double getXPos() const { return xPos; };
	double getYPos() const { return yPos; };
	double getZPos() const { return zPos; };
	double getXMom() const { return xMom; };
	double getYMom() const { return yMom; };
	double getZMom() const { return zMom; };
	double getKinEnergy() const { return kinEnergy; };

	//------------------------------------------------------------------------
	// methods to set private variable 

	void setParticleName(const std::string &partName) 
								{ particleName = partName; };
	
	void setAtomNum(const unsigned long &atomicNum) 
								{ atomNum = atomicNum; };
		
	void setNuclNum(const unsigned long &nucleonNum) 
								{ nuclNum = nucleonNum; };


	void setXPos(const double &xPosition) 
							{ xPos = xPosition; };

	void setYPos(const double &yPosition) 
							{ yPos = yPosition; };

	void setZPos(const double &zPosition) 
							{ zPos = zPosition; };

	void setXMom(const double &xMomentum) 
							{ xMom = xMomentum; };

	void setYMom(const double &yMomentum) 
							{ yMom = yMomentum; };

	void setZMom(const double &zMomentum) 
							{ zMom = zMomentum; };

	void setKinEnergy(const double &kineticEnergy) 
							{ kinEnergy = kineticEnergy; };
	
	void setAll( const std::string &particleName,		// particle name
				 const unsigned int &atomicNum,	// atomic number
				 const unsigned int &nucleonNum,	// nucleon number
				 const double &xPosition,		// X position 
				 const double &yPosition,		// Y position
				 const double &zPosition,		// Z position
				 const double &xMomentum,		// x-momentum
				 const double &yMomentum,		// y-momentum
				 const double &zMomentum,		// z-momentum
				 const double &kineticEnergy);	// kinetic energy

	//------------------------------------------------------------------------

	friend std::ostream &operator <<(std::ostream &os, const PiotPhaseSpace &pps)
	{
		os << pps.particleName << "\t" << pps.atomNum << "\t" 
		   << pps.nuclNum << "\t" << pps.xPos  << "\t" << pps.yPos    << "\t"
		   << pps.zPos    << "\t" << pps.xMom  << "\t" << pps.yMom    << "\t" 
		   << pps.zMom    << "\t" << pps.kinEnergy;

		return os;
	}

	friend std::istream &operator >>(std::istream &is, PiotPhaseSpace &pps)
	{
		is >> pps.particleName >> pps.atomNum >> pps.nuclNum 
	  	   >> pps.xPos  >> pps.yPos  >> pps.zPos    >> pps.xMom
		   >> pps.yMom  >> pps.zMom  >> pps.kinEnergy;

		return is;
	}


	//------------------------------------------------------------------------

  private:
	std::string particleName;
	unsigned int atomNum;
	unsigned int nuclNum;

	double  xPos;
	double  yPos;
	double  zPos;
	double  xMom;
	double  yMom;
	double  zMom;
	double  kinEnergy;

};

//============================================================================
// inline Functiions

inline void PiotPhaseSpace::setAll( 
             const std::string &partName,       // particle Name
             const unsigned int &atomicNum,      // atomic number
             const unsigned int &nucleonNum,     // nucleon number
             const double &xPosition,            // X position 
             const double &yPosition,            // Y position
             const double &zPosition,            // Z position
             const double &xMomentum,            // x-momentum
             const double &yMomentum,            // y-momentum
             const double &zMomentum,            // z-momentum
             const double &kineticEnergy)            // kinetic energy
{
    setParticleName(partName);
    setAtomNum(atomicNum);
    setNuclNum(nucleonNum);
    setXPos(xPosition);
    setYPos(yPosition);
    setZPos(zPosition);
    setXMom(xMomentum);
    setYMom(yMomentum);
    setZMom(zMomentum);
    setKinEnergy(kineticEnergy);
}

#endif // PiotPhaseSpace

