
### SEGM branch

This branch has a different behaviour in noise management that improves the parallel working of PosSim and Time* classes.

The problem of this branch is execution time.
The new code is very slow and the cause could be the deletion and reallocation of TGraph private members of TimeSegm::Group (it is a nested class) in every Digitization step or the use of TGraph objects themselves.
Also the use of Group::Sort() could slow the program.

The modifications done would not have broken the code, so TimeSegm should work fine; anyway segm C was added afterwards and not tested properly yet.
Becuase of this, TimeSegm class has to be used keeping in mind it could still have some unexpected behaviours.

At this point TimeSegm has not implemented completely yet because to see the real effects of segmentation, measures has to be computed at the end of each event and not every step, as set now.
This different behaviour was the original one of PosSim measures (look at legacy branch) and it will make PosSim measures work again theorically.
The main problem is that the measure saving process in the output root file has to be done at the end of every event instead of each hit as for the measurement.
