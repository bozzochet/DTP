
###SEGM branch

Segm branch implements TimeSegm class in Digitization.cpp.
This branch has a different behaviour in noise management that improves the parallel working of PosSim and Time* classes.
PosSim is not working properly in master branch while the code in here for Digitization.cpp returns to the orginal behaviour (is the one in _legacy_ branch) where PosSim was working well.

The problem of this branch is execution time.
The new code is very slow and the cause could be the deletion and reallocation of TGraph private members of TimeSegm::Group (it is a nested class) in every Digitization step or the use of TGraph objects themselves.
The modifications done would not have broken the code, so TimeSegm should work fine; anyway segm C was added afterwards and not tested properly yet.
Becuase of this, TimeSegm class has to be used keeping in mind it could still have some unexpected behaviours.
