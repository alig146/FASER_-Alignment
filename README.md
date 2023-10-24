This Repository contains the scripts needed to run FASER tracker alignment.


Get Calypso:
```
git clone –recursive https://gitlab.cern.ch/faser/calypso.git (or from your fork)  
git remote add [name1] https://gitlab.cern.ch/keli/calypso.git
git fetch [name1]
git checkout -b name2 [name1]/globalalign
```

Setup Environment:
```
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'

setupATLAS
asetup --input=calypso/asetup.faser Athena,22.0.49 
```

Compile Calypso:
```
mkdir build run
cd build
cmake -DINSTALL_CONDB=ON -DCMAKE_INSTALL_PREFIX=../run ../calypso
make -j8
make install
```


For Millepede-II installation see [here](https://www.desy.de/~kleinwrt/MP2/doc/html/index.html). 
