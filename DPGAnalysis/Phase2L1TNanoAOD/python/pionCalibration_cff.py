"""Optional pion calibration tables for a configured L1 emulation workflow."""

import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var

# Decoded barrel clusters before the hadronic response correction.
rawBarrelCaloTable = cms.EDProducer(
    "SimpleCandidateFlatTableProducer",
    src = cms.InputTag("l1tLayer1BarrelRaw", "DecodedHadClusters"),
    cut = cms.string(""),
    name = cms.string("RawBarrelCalo"),
    doc = cms.string("Decoded barrel hadronic clusters with hadronic correction disabled"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        pt = Var("pt", float, doc = "Decoded transverse momentum [GeV]", precision = -1),
        eta = Var("eta", float, doc = "Global pseudorapidity", precision = -1),
        phi = Var("phi", float, doc = "Global azimuth [rad]", precision = -1),
    ),
)

# Reuse the decoded calorimeter extension for hardware fields, including EM pt.
# The encoded EM pt must not be replaced by a value inferred from H/E.
rawBarrelCaloExtTable = cms.EDProducer(
    "L1PFDecodedCaloTableProducer",
    src = rawBarrelCaloTable.src,
    name = rawBarrelCaloTable.name,
    cut = cms.string(""),
)

# The same decoded objects with the configured hadronic response correction.
barrelCaloTable = rawBarrelCaloTable.clone(
    src = cms.InputTag("l1tLayer1Barrel", "DecodedHadClusters"),
    name = "BarrelCalo",
    doc = "Decoded barrel hadronic clusters with configured correction",
)
barrelCaloExtTable = rawBarrelCaloExtTable.clone(
    src = barrelCaloTable.src,
    name = barrelCaloTable.name,
)

# Generator particles before the response table's gen-jet matching.
calibGenPiTable = cms.EDProducer(
    "SimpleGenParticleFlatTableProducer",
    src = cms.InputTag("genParticles"),
    cut = cms.string("abs(pdgId) == 211 && statusFlags().isPrompt()"),
    name = cms.string("CalibGenPi"),
    doc = cms.string("Prompt charged generator pions before gen-jet matching"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        pt = Var("pt", float, doc = "Generator transverse momentum [GeV]", precision = -1),
        eta = Var("eta", float, doc = "Generator pseudorapidity", precision = -1),
        phi = Var("phi", float, doc = "Generator azimuth [rad]", precision = -1),
        mass = Var("mass", float, doc = "Generator mass [GeV]", precision = -1),
        pdgId = Var("pdgId", int, doc = "Signed PDG identifier"),
        charge = Var("charge", int, doc = "Electric charge in units of e"),
        vx = Var("vx", float, doc = "Production vertex x coordinate [cm]", precision = -1),
        vy = Var("vy", float, doc = "Production vertex y coordinate [cm]", precision = -1),
        vz = Var("vz", float, doc = "Production vertex z coordinate [cm]", precision = -1),
    ),
)

# Matching, propagation and cone sums require several collections, rather than
# individual-object expressions. All sources intentionally resolve the latest
# matching process. All response columns share the matched pion rows.
# Collection names containing "Calo", "Ecal" or "Hcal" use propagated directions,
# following the original FastPUPPI response calculation.
calibResponseTable = cms.EDProducer(
    "PionCalibrationResponseTableProducer",
    name = cms.string("Response"),
    genJets = cms.InputTag("ak4GenJetsNoNu"),
    genParticles = cms.InputTag("genParticles"),
    objects = cms.PSet(
        L1RawBarrelCalo = cms.VInputTag(cms.InputTag("l1tLayer1BarrelRaw", "DecodedHadClusters")),
        L1Calo = cms.VInputTag(cms.InputTag("l1tLayer1", "Calo")),
    ),
)

pionCalibrationTask = cms.Task(
    rawBarrelCaloTable,
    rawBarrelCaloExtTable,
    barrelCaloTable,
    barrelCaloExtTable,
    calibGenPiTable,
    calibResponseTable,
)


def add_pion_calibration_tables(process, path="p", output="outnano"):
    """Attach the calibration task and keep its tables in the selected output."""
    process.load("DPGAnalysis.Phase2L1TNanoAOD.pionCalibration_cff")
    getattr(process, path).associate(process.pionCalibrationTask)
    getattr(process, output).outputCommands.append(
        "keep nanoaodFlatTable_*Table_*_" + process.name_()
    )
    return process
