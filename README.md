This Repository contains the scripts needed to run FASER tracker alignment.


How to compile Calypso:
```
Checkout calypso
git clone –recursive ssh://git@gitlab.cern.ch:7999/faser/calypso.git or from your fork  
Add ssh://git@gitlab.cern.ch:7999/keli/calypso.git
Checkout the branch: “globalalign”
Make “build” and “run” folder
Compile calypso:
In build folder, do:
cmake -DINSTALL_CONDB=ON -DCMAKE_INSTALL_PREFIX=../run ../calypso
Then make & make install
```
