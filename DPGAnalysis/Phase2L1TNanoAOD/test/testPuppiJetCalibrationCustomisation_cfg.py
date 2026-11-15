"""Check the cmsDriver-facing PUPPI jet calibration customisation."""

import FWCore.ParameterSet.Config as cms

from DPGAnalysis.Phase2L1TNanoAOD.customisePuppiJetCalibration import (
    customisePuppiJetCalibrationNano,
)


process = cms.Process("L1TReRun")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(0))
process.nanoAOD_step = cms.Path()
process.NANOAODoutput = cms.OutputModule(
    "NanoAODOutputModule",
    fileName = cms.untracked.string("testPuppiJetCalibrationCustomisation.root"),
    outputCommands = cms.untracked.vstring("drop *"),
)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

process = customisePuppiJetCalibrationNano(process)

assert process.sc4RawJetTable.src == cms.InputTag("l1tSC4PFL1PuppiEmulator")
assert process.sc4RawJetTable.name.value() == "L1puppiJetSC4Raw"
assert "sc4RawJetTable" in process.puppiJetCalibrationTask.moduleNames()
assert (
    "keep nanoaodFlatTable_sc4RawJetTable_*_L1TReRun"
    in process.NANOAODoutput.outputCommands
)
