#! /bin/bash


runs=$1
fwhm=$2
energy1=$3 
NofParticles1=$4
energy2=$5 
NofParticles2=$6
energy3=$7 
NofParticles3=$8
ident=${9}

datafolder=data-`date +%Y%m%d`
mkdir ./${datafolder}

for (( i=0; i<=${runs}; i++))
do

simID=$[ ( $RANDOM % 10000 )  + 1 ]
ergfile1=${energy1}MeV${simID}sID${i}${ident}
macfile1=${ergfile1}.mac
logfile1=${ergfile1}.log
echo "/gun/particle proton " >> ${macfile1}
echo "/beam/energy/meanEnergy ${energy1} MeV" >> ${macfile1}
echo "/beam/energy/EnergyFWHM ${fwhm}" >> ${macfile1}
echo "/random/setSeeds ${simID} 2" >> ${macfile1}
echo "/run/beamOn ${NofParticles1}" >> ${macfile1}

./SiDetPhant ${macfile1} > ./${datafolder}/${logfile1}
mv Ergfile.dat ./${datafolder}/${ergfile1}.dat
touch Ergfile.dat
mv ${macfile1} ./${datafolder}/${macfile1}

simID=$[ ( $RANDOM % 10000 )  + 1 ]
ergfile2=${energy2}MeV${simID}sID${i}${ident}
macfile2=${ergfile2}.mac
logfile2=${ergfile2}.log
echo "/gun/particle proton " >> ${macfile2}
echo "/beam/energy/meanEnergy ${energy2} MeV" >> ${macfile2}
echo "/beam/energy/EnergyFWHM ${fwhm}" >> ${macfile2}
echo "/random/setSeeds ${simID} 2" >> ${macfile2}
echo "/run/beamOn ${NofParticles2}" >> ${macfile2}

./SiDetPhant ${macfile2} > ./${datafolder}/${logfile2}
mv Ergfile.dat ./${datafolder}/${ergfile2}.dat
touch Ergfile.dat
mv ${macfile2} ./${datafolder}/${macfile2}

simID=$[ ( $RANDOM % 10000 )  + 1 ]
ergfile3=${energy3}MeV${simID}sID${i}${ident}
macfile3=${ergfile3}.mac
logfile3=${ergfile3}.log
echo "/gun/particle proton " >> ${macfile3}
echo "/beam/energy/meanEnergy ${energy3} MeV" >> ${macfile3}
echo "/beam/energy/EnergyFWHM ${fwhm}" >> ${macfile3}
echo "/random/setSeeds ${simID} 2" >> ${macfile3}
echo "/run/beamOn ${NofParticles3}" >> ${macfile3}

./SiDetPhant ${macfile3} > ./${datafolder}/${logfile3}
mv Ergfile.dat ./${datafolder}/${ergfile3}.dat
touch Ergfile.dat
mv ${macfile3} ./${datafolder}/${macfile3}
done
