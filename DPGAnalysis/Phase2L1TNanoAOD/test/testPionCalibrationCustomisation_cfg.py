"""Check the cmsDriver-facing pion calibration customisation."""

import FWCore.ParameterSet.Config as cms

from DPGAnalysis.Phase2L1TNanoAOD.customisePionCalibration import (
    customisePionCalibrationNano,
)

process = cms.Process("L1TReRun")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(0))
process.load("L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_cff")
process.nanoAOD_step = cms.Path()
process.NANOAODoutput = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string("testPionCalibrationCustomisation.root"),
    outputCommands=cms.untracked.vstring("drop *"),
)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

process = customisePionCalibrationNano(process)

assert process.l1tLayer1BarrelRaw.gctHadInputConversionParameters.gctHadCorrector.value() == ""
assert process.l1tLayer1Barrel.gctHadInputConversionParameters.gctHadCorrector.value() != ""
assert process.l1tLayer1HGCalRaw.hgcalInputConversionParameters.corrector.value() == ""
assert process.l1tLayer1HGCalNoTKRaw.hgcalInputConversionParameters.corrector.value() == ""
assert process.l1tLayer1HGCal.hgcalInputConversionParameters.corrector.value() != ""
assert process.l1tLayer1HGCalNoTK.hgcalInputConversionParameters.corrector.value() != ""
for module in (
    "l1tLayer1BarrelRaw",
    "l1tLayer1HGCalRaw",
    "l1tLayer1HGCalNoTKRaw",
    "rawHGCalTable",
    "rawHGCalEMTable",
    "hgcalTable",
):
    assert module in process.pionCalibrationTask.moduleNames()
assert set(process.calibResponseTable.particlePdgIds) == {11, 22, 111, 130, 211}
assert process.rawBarrelCaloTable.src.getProcessName() == ""
assert process.barrelCaloTable.src.getProcessName() == ""
for tags in process.calibResponseTable.objects.parameters_().values():
    for tag in tags:
        assert tag.getProcessName() == ""
assert (
    "keep nanoaodFlatTable_*Table_*_L1TReRun"
    in process.NANOAODoutput.outputCommands
)
