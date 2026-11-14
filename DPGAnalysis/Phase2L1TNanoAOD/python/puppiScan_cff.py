"""Optional, unfiltered Layer-1 PF and PUPPI tables for offline scan studies."""

import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import ExtVar, Var
from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nanotables_cff import vtxTable as _vtxTable


# No Nano-level selection or float mantissa reduction is applied. These are
# scalar accessors from the existing EDM products, not a PUPPI reimplementation.
puppiScanPFTable = cms.EDProducer(
    "SimpleTriggerL1PFCandidateFlatTableProducer",
    src = cms.InputTag("l1tLayer1Barrel", "PF"),
    cut = cms.string(""),
    name = cms.string("PuppiScanPF"),
    doc = cms.string("Barrel fiducial PF candidates before PUPPI; retains upstream PF selection"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        pt = Var("pt", float, doc = "Unweighted PF transverse momentum [GeV]", precision = -1),
        eta = Var("eta", float, doc = "Global PF eta; vertex direction for charged candidates", precision = -1),
        phi = Var("phi", float, doc = "Global PF phi [rad]; vertex direction for charged candidates", precision = -1),
        caloEta = Var("caloEta()", float, doc = "Global calorimeter eta stored by the Layer-1 PF producer", precision = -1),
        caloPhi = Var("caloPhi()", float, doc = "Global calorimeter phi [rad] stored by the Layer-1 PF producer", precision = -1),
        pdgId = Var("pdgId", int, doc = "Signed PDG identifier"),
        charge = Var("charge", int, doc = "Electric charge in units of e"),
        z0 = Var("z0()", float, doc = "Charged-candidate longitudinal vertex position [cm]; zero for neutrals", precision = -1),
        hwZ0 = Var("hwZ0()", int, doc = "Signed PF longitudinal vertex code, 0.05 cm/count; zero for neutrals"),
        hwPt = Var("hwPt()", int, doc = "PF transverse-momentum code, 0.25 GeV/count"),
        hwEta = Var("hwEta()", int, doc = "Region-local calorimeter eta code; not global eta"),
        hwPhi = Var("hwPhi()", int, doc = "Region-local calorimeter phi code; not global phi"),
        puppiWeight = Var("puppiWeight()", float, doc = "Stored PUPPI weight; one for the unweighted PF product", precision = -1),
        hwPuppiWeight = Var("hwPuppiWeight()", int, doc = "Stored hardware PUPPI weight; unfilled zero for the PF product"),
    ),
    externalVariables = cms.PSet(
        sumRaw = ExtVar(cms.InputTag("l1tLayer1Barrel", "puppiSumRaw"), int,
            doc = "Actual 12-bit fixed-point neutral alpha sum used by PUPPI; -1 for charged or ineligible candidates"),
        region = ExtVar(cms.InputTag("l1tLayer1Barrel", "puppiRegion"), int,
            doc = "Index of the originating Layer-1 PF region in the configured region order"),
        slot = ExtVar(cms.InputTag("l1tLayer1Barrel", "puppiSlot"), int,
            doc = "Combined charged/neutral PUPPI input slot; neutral slots start at configured nTrack"),
        association = ExtVar(cms.InputTag("l1tLayer1Barrel", "puppiAssociation"), int,
            doc = "Charged primary-vertex association decision (0 or 1); -1 for neutral or ineligible candidates"),
    ),
)

# The PUPPI products contain surviving candidates after weighting, thresholding
# and regional sorting. Their calorimeter coordinate accessors are unfilled
# constructor zeros, and must not be interpreted as a measured direction.
puppiScanRefTable = puppiScanPFTable.clone(
    src = cms.InputTag("l1tLayer1Barrel", "Puppi"),
    name = cms.string("PuppiScanRef"),
    doc = cms.string("Barrel PUPPI reference after configured weighting, cuts and regional sorting"),
    externalVariables = cms.PSet(),
    variables = dict(
        pt = Var("pt", float, doc = "Weighted PUPPI transverse momentum [GeV]", precision = -1),
        eta = Var("eta", float, doc = "Global output eta; vertex direction for charged candidates", precision = -1),
        phi = Var("phi", float, doc = "Global output phi [rad]; vertex direction for charged candidates", precision = -1),
        caloEta = Var("caloEta()", float, doc = "Unfilled PUPPI accessor: constructor zero, not a calorimeter direction", precision = -1),
        caloPhi = Var("caloPhi()", float, doc = "Unfilled PUPPI accessor: constructor zero, not a calorimeter direction", precision = -1),
        hwEta = Var("hwEta()", int, doc = "Global output eta code, pi/720 per count"),
        hwPhi = Var("hwPhi()", int, doc = "Global output phi code, pi/720 radians per count"),
        puppiWeight = Var("puppiWeight()", float, doc = "Stored PUPPI weight; charged candidates have weight one", precision = -1),
        hwPuppiWeight = Var("hwPuppiWeight()", int, doc = "Stored integer weight accessor, not raw fractional weight bits; use puppiWeight for the weight"),
    ),
)

puppiScanFinalTable = puppiScanRefTable.clone(
    src = cms.InputTag("l1tLayer2Deregionizer", "Puppi"),
    name = cms.string("PuppiScanFinal"),
    doc = cms.string("Full-detector PUPPI reference after the configured Layer-2 deregionizer"),
)

# Keep every row in original collection order, including zero-sumPt vertices.
# The customisation uses the exact barrel vtxCollection InputTag. Offline code
# reproduces Layer-1 selection and its separate 0.05 cm quantisation; hwZ0 here
# is the original unsigned VertexWord bit pattern, NOT a signed PF z0 code.
puppiScanVertexTable = _vtxTable.clone(
    src = cms.InputTag("l1tVertexFinderEmulator", "L1VerticesEmulation"),
    name = cms.string("PuppiScanVertex"),
    doc = cms.string("All barrel Layer-1 input VertexWords in collection order for primary-vertex reconstruction"),
    extension = cms.bool(False),
    variables = dict(
        z0 = Var("z0()", float, doc = "Input vertex z coordinate [cm], before PF requantisation", precision = -1),
        sumPt = Var("pt()", float, doc = "Input vertex track sum pT [GeV] used by the Layer-1 selection", precision = -1),
        hwZ0 = Var("z0Bits()", "uint", doc = "Unsigned 15-bit VertexWord z0 pattern; signed two's complement with 1/512 cm/count"),
    ),
)

puppiScanTask = cms.Task(
    puppiScanPFTable,
    puppiScanRefTable,
    puppiScanFinalTable,
    puppiScanVertexTable,
)


def add_puppi_scan_tables(process, path="p", output="outnano"):
    """Attach the opt-in table task to an existing compatible L1 workflow."""
    for label in ("puppiScanPFTable", "puppiScanRefTable", "puppiScanFinalTable", "puppiScanVertexTable", "puppiScanTask"):
        if hasattr(process, label):
            raise ValueError("add_puppi_scan_tables will not replace process." + label)
    target_path = getattr(process, path)
    target_output = getattr(process, output)
    process.load("DPGAnalysis.Phase2L1TNanoAOD.puppiScan_cff")
    target_path.associate(process.puppiScanTask)
    target_output.outputCommands.append(
        "keep nanoaodFlatTable_puppiScan*Table_*_" + process.name_()
    )
    return process
