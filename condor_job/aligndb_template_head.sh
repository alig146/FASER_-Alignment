#!/bin/sh

rm -f *_Align.pool.root
rm -f *_ALLP200.db
rm -f PoolFileCatalog.xml
rm -f writeAlignment_*.log

python WriteAlignmentConfig_Faser01.py 
