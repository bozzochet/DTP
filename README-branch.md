
###LEGACY branch

PosSim class was created merging together code parts in digitization and analysis src files.
Originally the measures were done at the end of every event while at the moment, Digitization.cpp stores measures at each step (i.e. at each event's single hit).
After this change PosSim does not returns the same results and because of this the previous asset is maintained in this branch, to be able the check the differences between this and the present state of the project.

_Segm_ branch could return the code to the previous behaviour for PosSim but it needs substantial modification to the saving procedure of the output root file.
