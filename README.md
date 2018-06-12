# SiDetPhant
Geant4 simulation of silicon PIN diodes in PMMA phantom



setup:

2 sensors close to each other, both 1.5 cm away from the centre of the phantom, on the surface facing the source.



Room geometry:

irregular, simplified to a normal box with weighting the shortest and the longest width and depth: 1/3 weight given by shortest, 2/3 given by longest due to the geometry.

if we as width the shorter horizonzal axis, 3.80 mt - 4 mt ---> weighted avg is 3.93 mt
for the depth 12.88-13.07 ---> 13,01 mt
height is 4 mt



for the column blocking the door: 
avg distance from wall opposite to neutron source: 360 cm
avg distance from wall on side of neutron source: 883 cm
column size: 58 cm along depth, 80 cm along width

distance from centre of the room of the centre of the colum: 6.501 mt (half depth) - 3.6 mt (space wall-column) - 0.58 m /2 (half of column size along that axis)
distance along widhth: widht - column size all divided by 2


4mm thick PVC floor above the concrete.

to add: steel bars










sphere for scorer was added:
need to implement how to track the particle crossing the section




idea behind the sphere scorer:
1)get binning values to be read and put into a vector
2) link the volume to a multifunctional detector
3) link i-primitive with energy interval filter, one for each binning level
4) record the neutron crossing the area in each level of energy
5) basically get the number of neutrons vs the energy level just like this
6) divide by the surface area to get the fluence 
