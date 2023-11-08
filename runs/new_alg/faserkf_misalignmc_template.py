#!/usr/bin/env python

import sys
import os
import time
from AthenaCommon.Logging import log, logging
from AthenaCommon.Constants import INFO
from AthenaCommon.Configurable import Configurable
from CalypsoConfiguration.AllConfigFlags import ConfigFlags
from AthenaConfiguration.TestDefaults import defaultTestFiles
from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
from AthenaPoolCnvSvc.PoolWriteConfig import PoolWriteCfg
from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
from TrackerPrepRawDataFormation.TrackerPrepRawDataFormationConfig import FaserSCT_ClusterizationCfg
from TrackerSpacePointFormation.TrackerSpacePointFormationConfig import TrackerSpacePointFinderCfg
from TrackerSegmentFit.TrackerSegmentFitConfig import SegmentFitAlgCfg
from FaserActsKalmanFilter.GhostBustersConfig import GhostBustersCfg
from TruthSeededTrackFinder.TruthSeededTrackFinderConfig import TruthSeededTrackFinderCfg
from FaserSCT_GeoModel.FaserSCT_GeoModelConfig import FaserSCT_GeometryCfg
from FaserActsGeometry.ActsGeometryConfig import ActsTrackingGeometryToolCfg
from FaserActsGeometry.ActsGeometryConfig import ActsExtrapolationToolCfg
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from MagFieldServices.MagFieldServicesConfig import MagneticFieldSvcCfg
from AthenaConfiguration.ComponentFactory import CompFactory

THistSvc=CompFactory.getComps("THistSvc")

GlobalChi2Alignment, FaserActsExtrapolationTool, FaserActsTrackingGeometryTool=CompFactory.getComps("GlobalChi2Alignment", "FaserActsExtrapolationTool","FaserActsTrackingGeometryTool")
from FaserActsGeometry.ActsGeometryConfig import ActsTrackingGeometrySvcCfg

def FaserActsAlignmentCondAlgCfg(flags, **kwargs):
   acc = ComponentAccumulator()
   acc.addCondAlgo(CompFactory.FaserActsAlignmentCondAlg(name="FaserActsAlignmentCondAlg", **kwargs))
   return acc

def GlobalChi2Alignment_OutputESDCfg(flags):
    acc = ComponentAccumulator()
    itemList = ["xAOD::EventInfo#*",
                "xAOD::EventAuxInfo#*",
                "TrackCollection#*"
                ]

    acc.merge(OutputStreamCfg(flags, "ESD",itemList))
    ostream = acc.getEventAlgo("OutputStreamESD")
    ostream.TakeItemsFromInput = True
    return acc

def GlobalChi2AlignmentCfg(flags, **kwargs):
    acc = FaserSCT_GeometryCfg(flags)
    acc.merge(MagneticFieldSvcCfg(flags))
    acts_tracking_geometry_svc = ActsTrackingGeometrySvcCfg(flags)
    acc.merge(acts_tracking_geometry_svc )
    track_seed_tool = CompFactory.CircleFitTrackSeedTool()
    track_seed_tool.removeIFT = True
    sigma_loc0 = 1.9e-2
    sigma_loc1 = 9e-1
    sigma_phi = 3.3e-2
    sigma_theta = 2.9e-4
    p = 1000
    sigma_p = 0.1 * p
    sigma_qop = sigma_p / (p * p)
    initial_variance_inflation = [100, 100, 100, 100, 1000]
    track_seed_tool.covLoc0 = initial_variance_inflation[0] * sigma_loc1 * sigma_loc1
    track_seed_tool.covLoc1 = initial_variance_inflation[1] * sigma_loc0 * sigma_loc0
    track_seed_tool.covPhi = initial_variance_inflation[2] * sigma_phi * sigma_phi
    track_seed_tool.covTheta = initial_variance_inflation[3] * sigma_theta * sigma_theta
    track_seed_tool.covQOverP = initial_variance_inflation[4] * sigma_qop * sigma_qop
    track_seed_tool.std_cluster = 0.0231
    track_seed_tool.TrackCollection = "Segments"
    trajectory_states_writer_tool = CompFactory.RootTrajectoryStatesWriterTool()
    trajectory_states_writer_tool.noDiagnostics = kwargs.pop("noDiagnostics", True)
    trajectory_states_writer_tool1 = CompFactory.RootTrajectoryStatesWriterTool()
    trajectory_states_writer_tool1.noDiagnostics = kwargs.pop("noDiagnostics", True)
    trajectory_states_writer_tool1.FilePath = "/pool/track_states_ckf1_10root"
    
    trajectory_summary_writer_tool = CompFactory.RootTrajectorySummaryWriterTool()
    trajectory_summary_writer_tool.noDiagnostics = kwargs.pop("noDiagnostics", True)
    trajectory_summary_writer_tool1 = CompFactory.RootTrajectorySummaryWriterTool()
    trajectory_summary_writer_tool1.FilePath = "/pool/track_summary_ckf1_0.root"
    trajectory_summary_writer_tool1.noDiagnostics = kwargs.pop("noDiagnostics", True)
    actsExtrapolationTool = CompFactory.FaserActsExtrapolationTool("FaserActsExtrapolationTool")
    actsExtrapolationTool.MaxSteps = 1000
    result, actsTrackingGeometryTool = ActsTrackingGeometryToolCfg(flags)
    actsExtrapolationTool.TrackingGeometryTool = actsTrackingGeometryTool
    acc.merge(result)
    #actsExtrapolationTool.TrackingGeometryTool = CompFactory.FaserActsTrackingGeometryTool("TrackingGeometryTool")
    
    trajectory_performance_writer_tool = CompFactory.PerformanceWriterTool("PerformanceWriterTool")
    trajectory_performance_writer_tool.ExtrapolationTool = actsExtrapolationTool
    trajectory_performance_writer_tool.noDiagnostics = kwargs.pop("noDiagnostics", True)
    #track_finder_tool = CompFactory.SPSimpleTrackFinderTool()
    #track_finder_tool = CompFactory.SegmentFitClusterTrackFinderTool()
    globalchi2fit = CompFactory.GlobalChi2Alignment(**kwargs)
    globalchi2fit.BiasedResidual = False
    globalchi2fit.ExtrapolationTool = actsExtrapolationTool
    globalchi2fit.TrackingGeometryTool=actsTrackingGeometryTool
    kalman_fitter1 = CompFactory.KalmanFitterTool(name="fitterTool1")
    kalman_fitter1.noDiagnostics = True
    kalman_fitter1.ActsLogging = "INFO"
    kalman_fitter1.SummaryWriter = False
    kalman_fitter1.StatesWriter = False
    kalman_fitter1.SeedCovarianceScale = 10
    kalman_fitter1.isMC = True
    kalman_fitter1.RootTrajectoryStatesWriterTool = trajectory_states_writer_tool1
    kalman_fitter1.RootTrajectorySummaryWriterTool = trajectory_summary_writer_tool1
    globalchi2fit.KalmanFitterTool1 = kalman_fitter1
    globalchi2fit.TrackSeed = track_seed_tool
    globalchi2fit.ActsLogging = "INFO"
    globalchi2fit.isMC = True
    globalchi2fit.nMax = 10
    globalchi2fit.chi2Max = 100000
    globalchi2fit.maxSteps = 5000
    histSvc= CompFactory.THistSvc()
    histSvc.Output +=  [ "GlobalChi2Alignment DATAFILE='kfalignment_mc_666NAME666.root' OPT='RECREATE'"]
    acc.addService(histSvc)
    acc.addEventAlgo(globalchi2fit)
#    acc.merge(CKF2GlobalAlignment_OutputESDCfg(flags))
    return acc

if __name__ == "__main__":
   a = time.time()

   log.setLevel(INFO)
   Configurable.configurableRun3Behavior = True
   
   # Configure
   #ConfigFlags.Input.Files = [ '/eos/home-k/keli/Faser/skim/8023/skim_8023_666JOB666.DAOD.pool.root']
#   ConfigFlags.Input.Files = [ 
#      '/eos/home-k/keli/Faser/skim//8023/skim_8023_0.DAOD.pool.root'
#      ,'/eos/home-k/keli/Faser/skim//8023/skim_8023_1.DAOD.pool.root'
#   ]
   ConfigFlags.Input.Files = [ '/data/agarabag/mc_data/singlemuon_newfieldmap/666NAME1666/my100GeV_1TeV_1mEvt_666NAME1666.RDO.pool.root','/data/agarabag/mc_data/singlemuon_newfieldmap/666NAME2666/my100GeV_1TeV_1mEvt_666NAME2666.RDO.pool.root','/data/agarabag/mc_data/singlemuon_newfieldmap/666NAME3666/my100GeV_1TeV_1mEvt_666NAME3666.RDO.pool.root','/data/agarabag/mc_data/singlemuon_newfieldmap/666NAME4666/my100GeV_1TeV_1mEvt_666NAME4666.RDO.pool.root' ]


   ConfigFlags.Input.ProjectName = "data21"
   ConfigFlags.GeoModel.Align.Dynamic = False
   ConfigFlags.GeoModel.FaserVersion = "FASERNU-03"
   ConfigFlags.IOVDb.GlobalTag = "OFLCOND-FASER-03"
   ConfigFlags.Detector.GeometryFaserSCT = True

#### FOR MC
   ConfigFlags.Input.isMC = True
   ConfigFlags.IOVDb.DatabaseInstance = "OFLP200"
#### FOR DATA
#   ConfigFlags.Input.isMC = False
#   ConfigFlags.IOVDb.DatabaseInstance = "CONDBR3"

   ConfigFlags.Beam.NumberOfCollisions = 0.
   ConfigFlags.addFlag("Input.InitialTimeStamp", 0)
   ConfigFlags.fillFromArgs(sys.argv[1:])
   #ConfigFlags.TrackingGeometry.MaterialSource = "material-maps.json"
   ConfigFlags.TrackingGeometry.MaterialSource = "/cvmfs/faser.cern.ch/repo/sw/database/DBRelease/current/acts/material-maps.json"
   ConfigFlags.Exec.MaxEvents = -1
   #ConfigFlags.Concurrency.NumThreads = 1
   ConfigFlags.lock()
   
   from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
   # Core components
   acc = MainServicesCfg(ConfigFlags)
   #acc.merge(PoolWriteCfg(ConfigFlags))
   acc.merge(FaserGeometryCfg(ConfigFlags))

   from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
   acc.merge(PoolReadCfg(ConfigFlags))


   acc.merge(FaserSCT_ClusterizationCfg(ConfigFlags))
   acc.merge(TrackerSpacePointFinderCfg(ConfigFlags))
   acc.merge(SegmentFitAlgCfg(ConfigFlags, SharedHitFraction=0.61, MinClustersPerFit=5, TanThetaXZCut=0.083))
   acc.merge(GhostBustersCfg(ConfigFlags))
###########################

#   acc.merge(FaserSCT_ClusterizationCfg(ConfigFlags))
#   acc.merge(FaserSCT_ClusterizationCfg(ConfigFlags, DataObjectName="SCT_xAODs", checkBadChannels=True))


   acc.merge(GlobalChi2AlignmentCfg(ConfigFlags))


#### FOR DATA
   # replicaSvc = acc.getService("DBReplicaSvc")
   # replicaSvc.COOLSQLiteVetoPattern = ""
   # replicaSvc.UseCOOLSQLite = True
   # replicaSvc.UseCOOLFrontier = False
   # replicaSvc.UseGeomSQLite = True


   logging.getLogger('forcomps').setLevel(INFO)
   acc.foreach_component("*").OutputLevel = INFO
   acc.foreach_component("*ClassID*").OutputLevel = INFO   
   #acc.getService("StoreGateSvc").Dump = True
   #acc.getService("ConditionStore").Dump = True
   #acc.printConfig(withDetails=True)
   #ConfigFlags.dump()
   
   # Reduce event loop printout
   # eventLoop = CompFactory.AthenaEventLoopMgr()
   # eventLoop.EventPrintoutInterval = 500
   # acc.addService(eventLoop)

   # print time to run
   b = time.time()
   from AthenaCommon.Logging import log
   log.info(f"Finish execution in {b-a} seconds")

   # Execute and finish
   sc = acc.run()
   if sc.isSuccess():
      log.info("Execution succeeded")
      sys.exit(0)
   else:
      log.info("Execution failed, return 1")
      sys.exit(1)
