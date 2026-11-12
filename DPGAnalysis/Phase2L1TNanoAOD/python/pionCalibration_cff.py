"""Optional single-particle calibration tables for an L1 emulation workflow."""

import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var


def _decoded_calo_tables(src, name, doc):
    """Return a familiar candidate table and its decoded-hardware extension."""
    table = cms.EDProducer(
        "SimpleCandidateFlatTableProducer",
        src = src,
        cut = cms.string(""),
        name = cms.string(name),
        doc = cms.string(doc),
        singleton = cms.bool(False),
        extension = cms.bool(False),
        variables = cms.PSet(
            pt = Var("pt", float, doc = "Decoded transverse momentum [GeV]", precision = -1),
            eta = Var("eta", float, doc = "Global pseudorapidity", precision = -1),
            phi = Var("phi", float, doc = "Global azimuth [rad]", precision = -1),
        ),
    )
    extension = cms.EDProducer(
        "L1PFDecodedCaloTableProducer",
        src = src,
        name = cms.string(name),
        cut = cms.string(""),
    )
    return table, extension


# Barrel decoded hadronic clusters before and after the response correction.
rawBarrelCaloTable, rawBarrelCaloExtTable = _decoded_calo_tables(
    cms.InputTag("l1tLayer1BarrelRaw", "DecodedHadClusters"),
    "RawBarrelCalo",
    "Decoded barrel hadronic clusters with hadronic correction disabled",
)
barrelCaloTable, barrelCaloExtTable = _decoded_calo_tables(
    cms.InputTag("l1tLayer1Barrel", "DecodedHadClusters"),
    "BarrelCalo",
    "Decoded barrel hadronic clusters with configured correction",
)

# HGCal decoded clusters in the tracker and no-tracker regions. The raw Layer-1
# clones have the common HGCal response corrector disabled. Hadronic and EM
# outputs are both retained because their cone-sum ratio defines the fit's EM
# fraction.
rawHGCalTable, rawHGCalExtTable = _decoded_calo_tables(
    cms.InputTag("l1tLayer1HGCalRaw", "DecodedHadClusters"),
    "RawHGCal",
    "Decoded HGCal hadronic clusters with correction disabled",
)
rawHGCalEMTable, rawHGCalEMExtTable = _decoded_calo_tables(
    cms.InputTag("l1tLayer1HGCalRaw", "DecodedEmClusters"),
    "RawHGCalEM",
    "Decoded HGCal electromagnetic clusters with correction disabled",
)
rawHGCalNoTKTable, rawHGCalNoTKExtTable = _decoded_calo_tables(
    cms.InputTag("l1tLayer1HGCalNoTKRaw", "DecodedHadClusters"),
    "RawHGCalNoTK",
    "Decoded no-track HGCal hadronic clusters with correction disabled",
)
rawHGCalNoTKEMTable, rawHGCalNoTKEMExtTable = _decoded_calo_tables(
    cms.InputTag("l1tLayer1HGCalNoTKRaw", "DecodedEmClusters"),
    "RawHGCalNoTKEM",
    "Decoded no-track HGCal electromagnetic clusters with correction disabled",
)
hgcalTable, hgcalExtTable = _decoded_calo_tables(
    cms.InputTag("l1tLayer1HGCal", "DecodedHadClusters"),
    "HGCal",
    "Decoded HGCal hadronic clusters with configured correction",
)
hgcalEMTable, hgcalEMExtTable = _decoded_calo_tables(
    cms.InputTag("l1tLayer1HGCal", "DecodedEmClusters"),
    "HGCalEM",
    "Decoded HGCal electromagnetic clusters with configured correction",
)
hgcalNoTKTable, hgcalNoTKExtTable = _decoded_calo_tables(
    cms.InputTag("l1tLayer1HGCalNoTK", "DecodedHadClusters"),
    "HGCalNoTK",
    "Decoded no-track HGCal hadronic clusters with configured correction",
)
hgcalNoTKEMTable, hgcalNoTKEMExtTable = _decoded_calo_tables(
    cms.InputTag("l1tLayer1HGCalNoTK", "DecodedEmClusters"),
    "HGCalNoTKEM",
    "Decoded no-track HGCal electromagnetic clusters with configured correction",
)

# Generator particles before the response table's gen-jet matching. This table
# remains pion-specific for compatibility; the Response table also accepts the
# mixed particle IDs required by the HGCal fit.
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

# Matching, propagation and cone sums require event-wide access to several
# collections. All sources intentionally resolve the latest matching process.
# The HGCal names preserve the original FastPUPPI convention, which evaluates
# those cones around the generator direction at the vertex.
calibResponseTable = cms.EDProducer(
    "PionCalibrationResponseTableProducer",
    name = cms.string("Response"),
    genJets = cms.InputTag("ak4GenJetsNoNu"),
    genParticles = cms.InputTag("genParticles"),
    particlePdgIds = cms.vint32(11, 22, 111, 130, 211),
    objects = cms.PSet(
        L1RawBarrelCalo = cms.VInputTag(
            cms.InputTag("l1tLayer1BarrelRaw", "DecodedHadClusters"),
        ),
        L1Calo = cms.VInputTag(cms.InputTag("l1tLayer1", "Calo")),
        L1RawHGCal = cms.VInputTag(
            cms.InputTag("l1tLayer1HGCalRaw", "DecodedHadClusters"),
            cms.InputTag("l1tLayer1HGCalNoTKRaw", "DecodedHadClusters"),
        ),
        L1RawHGCalEM = cms.VInputTag(
            cms.InputTag("l1tLayer1HGCalRaw", "DecodedEmClusters"),
            cms.InputTag("l1tLayer1HGCalNoTKRaw", "DecodedEmClusters"),
        ),
        L1HGCal = cms.VInputTag(
            cms.InputTag("l1tLayer1HGCal", "DecodedHadClusters"),
            cms.InputTag("l1tLayer1HGCalNoTK", "DecodedHadClusters"),
        ),
    ),
)

pionCalibrationTask = cms.Task(
    rawBarrelCaloTable,
    rawBarrelCaloExtTable,
    barrelCaloTable,
    barrelCaloExtTable,
    rawHGCalTable,
    rawHGCalExtTable,
    rawHGCalEMTable,
    rawHGCalEMExtTable,
    rawHGCalNoTKTable,
    rawHGCalNoTKExtTable,
    rawHGCalNoTKEMTable,
    rawHGCalNoTKEMExtTable,
    hgcalTable,
    hgcalExtTable,
    hgcalEMTable,
    hgcalEMExtTable,
    hgcalNoTKTable,
    hgcalNoTKExtTable,
    hgcalNoTKEMTable,
    hgcalNoTKEMExtTable,
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
