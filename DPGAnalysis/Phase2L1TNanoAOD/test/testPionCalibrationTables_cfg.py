"""Check calibration table construction without FastPUPPI configuration imports."""

import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(0))
process.load("DPGAnalysis.Phase2L1TNanoAOD.pionCalibration_cff")
process.tables = cms.Path(process.pionCalibrationTask)
process.out = cms.OutputModule(
    "NanoAODOutputModule",
    fileName = cms.untracked.string("testPionCalibrationTables.root"),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*_*_TEST"),
)
process.end = cms.EndPath(process.out)
