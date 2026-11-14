"""Construct the optional scalar tables and external ValueMap readers."""

import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(0))
process.load("DPGAnalysis.Phase2L1TNanoAOD.puppiScan_cff")
process.tables = cms.Path(process.puppiScanTask)
process.out = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string("testPuppiScanTables.root"),
    outputCommands=cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*_*_TEST"),
)
process.end = cms.EndPath(process.out)

# Reference tables must never consume maps keyed to the pre-PUPPI PF product.
assert not process.puppiScanRefTable.externalVariables.parameterNames_()
assert not process.puppiScanFinalTable.externalVariables.parameterNames_()
assert set(process.puppiScanPFTable.externalVariables.parameterNames_()) == {
    "sumRaw", "region", "slot", "association",
}
for table in (process.puppiScanPFTable, process.puppiScanRefTable, process.puppiScanFinalTable, process.puppiScanVertexTable):
    assert table.cut.value() == ""
    for variable in table.variables.parameters_().values():
        if variable.type.value() == "float":
            assert variable.precision.value() == -1
