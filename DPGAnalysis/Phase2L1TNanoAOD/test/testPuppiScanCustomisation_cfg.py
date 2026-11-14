"""Check cmsDriver attachment, exact vertex source and supported-mode guards."""

import FWCore.ParameterSet.Config as cms
from DPGAnalysis.Phase2L1TNanoAOD.customisePuppiScan import customisePuppiScanNano

process = cms.Process("L1TReRun")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(0))
process.load("L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_cff")
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer
process.l1tLayer2Deregionizer = l1tDeregionizerProducer.clone()
process.nanoAOD_step = cms.Path()
process.NANOAODoutput = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string("testPuppiScanCustomisation.root"),
    outputCommands=cms.untracked.vstring("drop *"),
)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Loading an opt-in fragment must preserve customised standard-menu modules.
from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nanotables_cff import vtxTable
process.vtxTable = vtxTable.clone(src=cms.InputTag("preserveUserVertices"))
original_vertex_table = process.vtxTable.dumpPython()

barrel = process.l1tLayer1Barrel
barrel.vtxCollection = cms.InputTag("alternateVertices", "selected", "INPUT")
original_parameters = barrel.puAlgoParameters.dumpPython()
process = customisePuppiScanNano(process)
assert process.vtxTable.dumpPython() == original_vertex_table
assert barrel.producePuppiDiagnostics.value()
assert barrel.puAlgoParameters.dumpPython() == original_parameters
assert process.puppiScanVertexTable.src == barrel.vtxCollection
assert "keep nanoaodFlatTable_puppiScan*Table_*_L1TReRun" in process.NANOAODoutput.outputCommands
assert set(process.puppiScanTask.moduleNames()) == {
    "puppiScanPFTable", "puppiScanRefTable", "puppiScanFinalTable", "puppiScanVertexTable",
}

# Each unsupported mode must fail before touching existing tables.
def requires_error(target, attribute, replacement, expected):
    existed = hasattr(target, attribute)
    original = getattr(target, attribute, None)
    setattr(target, attribute, replacement)
    try:
        try:
            customisePuppiScanNano(process)
        except ValueError as error:
            assert expected in str(error), str(error)
        else:
            raise AssertionError("Unsupported PUPPI mode was accepted")
    finally:
        if existed:
            setattr(target, attribute, original)
        else:
            delattr(target, attribute)

requires_error(barrel, "vtxCollectionEmulation", cms.bool(False), "VertexWord")
requires_error(barrel, "nVtx", cms.int32(2), "one primary vertex")
requires_error(barrel.puAlgoParameters, "useMLAssociation", cms.bool(True), "cut-based")
requires_error(barrel.puAlgoParameters, "fakePuppi", cms.bool(True), "fakePuppi=False")
requires_error(barrel.puAlgoParameters, "nOut", cms.uint32(1), "nOut == nIn")
try:
    customisePuppiScanNano(process)
except ValueError as error:
    assert "will not replace" in str(error)
else:
    raise AssertionError("Duplicate customisation was accepted")
