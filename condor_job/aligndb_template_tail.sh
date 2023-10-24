
AtlCoolCopy "sqlite://;schema=./FASER-01_ALLP200.db;dbname=OFLP200" "sqlite://;schema=data/sqlite200/ALLP200.db;dbname=OFLP200" -fnp -ot "TRACKER-ALIGN-01"
AtlCoolCopy "sqlite://;schema=./FASER-02_ALLP200.db;dbname=OFLP200" "sqlite://;schema=data/sqlite200/ALLP200.db;dbname=OFLP200" -fnp -ot "TRACKER-ALIGN-02"
mkdir data/poolcond
cp -n PoolFileCatalog.xml data/poolcond/
ln -s ${PWD}/data/poolcond/PoolFileCatalog.xml ${PWD}/data/poolcond/PoolCat_oflcond.xml
cp -n *_Align.pool.root data/poolcond/

