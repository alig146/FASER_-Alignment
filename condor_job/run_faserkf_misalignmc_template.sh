#!/bin/sh
workdir=`pwd`
echo 'start iteration ' $1 in $workdir
cd /home/agarabag/ift_alignment/ #use work dir
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
export ALRB_localConfigDir=$HOME/localConfig
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
asetup --input=calypso/asetup.faser Athena,22.0.49
cd run/
source ./setup.sh
export ATLAS_POOLCOND_PATH=${workdir}/data
echo 'pool path ' ${ATLAS_POOLCOND_PATH}
export EOS_MGM_URL=root://eosuser.cern.ch
export outputpath=/data/agarabag/alignment_output/tracking_mc_clus_algn_S1L0_iter4/
cd /home/agarabag/condor_jobs/actual_align_mc_test/
cp faserkf_misalignmc_666NAME666.py ${workdir}/
cp aligndb_template_*sh ${workdir}/
cp WriteAlignmentConfig_Faser0*py ${workdir}/
cp inputforalign.txt ${workdir}/
cp material-maps.json ${workdir}/
cd ${workdir}
mkdir ./data
mkdir ./data/sqlite200
mkdir ./result
mkdir ./result/666NAME666
cp /cvmfs/faser.cern.ch/repo/sw/database/DBRelease/current/sqlite200/ALLP200_221206.db ./data/sqlite200/ALLP200.db
cat aligndb_template_head.sh > ./aligndb_copy.sh
echo "python WriteAlignmentConfig_Faser02.py 'AlignDbTool.AlignmentConstants={ $(cat inputforalign.txt) }' >& writeAlignment_Faser02.log" >>./aligndb_copy.sh
cat aligndb_template_tail.sh >>./aligndb_copy.sh
chmod 755 ./aligndb_copy.sh
./aligndb_copy.sh >log
python faserkf_misalignmc_666NAME666.py
cp kfalignment_mc_666NAME666.root ${outputpath}/kfalignment_mc_666NAME666.root
tail -n 600 *out >> log
cp log ${outputpath}/log_666NAME666
