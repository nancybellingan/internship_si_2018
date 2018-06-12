#! /bin/bash


runs=${1}
fwhm=${2}
energy1=${3} 
NofParticles1=${4}
ident=${5}

datafolder=data-`date +%Y%m%d`
mkdir -p ./${datafolder}
touch Ergfile.dat

for (( i=0; i<=${runs}; i++))
do

simID=$[ ( $RANDOM % 100000 )  + 1 ]
simID2=$[ ( $RANDOM % 100000 )  + 1 ]
ergfile1=${energy1}MeV${simID}sID${i}${ident}
macfile1=${ergfile1}.mac
logfile1=${ergfile1}.log
echo "/gun/particle proton " >> ${macfile1}
echo "/gun/energy ${energy1} MeV" >> ${macfile1}
echo "/beam/energy/EnergyFWHM ${fwhm}" >> ${macfile1}
echo "/random/setSeeds ${simID} ${simID2}" >> ${macfile1}
echo "/run/beamOn ${NofParticles1}" >> ${macfile1}

./SiDetPhant ${macfile1} > ./${datafolder}/${logfile1}
mv Ergfile.dat ./${datafolder}/${ergfile1}.dat
touch Ergfile.dat
mv totED.txt ./${datafolder}/${ergfile1}.dep
mv ${macfile1} ./${datafolder}/${macfile1}

done
