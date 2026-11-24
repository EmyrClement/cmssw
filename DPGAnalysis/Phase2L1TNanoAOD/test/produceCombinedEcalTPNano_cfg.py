"""Produce legacy and Phase-2 ECAL TP/RCT/GCT tables in one Nano.

The input must contain same-event Phase-2 TPs and regenerated legacy TPs. Use
regenerateLegacyEcalTPs_cfg.py to create such an input from the central no-PU
D128 ECAL RelVal.
"""

from pathlib import Path

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9


options = VarParsing()
options.register(
    "outputFile",
    "combined_ecal_tp_nano.root",
    options.multiplicity.singleton,
    options.varType.string,
    "Nano output file",
)
options.register(
    "maxEvents", 100, options.multiplicity.singleton, options.varType.int, "Number of events"
)
options.register(
    "inputFiles", [], options.multiplicity.list, options.varType.string, "Input EDM files"
)
options.register(
    "legacyTPProcess",
    "LEGACY",
    options.multiplicity.singleton,
    options.varType.string,
    "Process containing simEcalEBTriggerPrimitiveDigis",
)
options.register(
    "phase2TPProcess",
    "DIGI2RAW",
    options.multiplicity.singleton,
    options.varType.string,
    "Process containing simEcalEBTriggerPrimitivePhase2Digis",
)
options.register(
    "genProcess",
    "DIGI2RAW",
    options.multiplicity.singleton,
    options.varType.string,
    "Process containing genParticles",
)
options.parseArguments()

if not options.inputFiles:
    raise ValueError("Specify inputFiles containing both legacy and Phase-2 ECAL TPs")
output = Path(options.outputFile.removeprefix("file:"))
if output.exists():
    raise ValueError(f"Refusing to overwrite {output}")

process = cms.Process("TPNANO", Phase2C26I13M9)
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D128Reco_cff")
process.CaloGeometryBuilder.SelectedCalos = cms.vstring("EcalBarrel", "HCAL")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")
process.load("CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi")
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, "150X_mcRun4_realistic_v1", "")
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
)

legacy_source = cms.InputTag("simEcalEBTriggerPrimitiveDigis", "", options.legacyTPProcess)
phase2_source = cms.InputTag(
    "simEcalEBTriggerPrimitivePhase2Digis", "", options.phase2TPProcess
)

from DPGAnalysis.Phase2L1TNanoAOD.tpNanotables_cff import (
    ecalEBPhase2Table,
    ecalEBTable,
    ecalGCTClusterTable,
    ecalRCTClusterTable,
)
from L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi import (
    l1tPhase2L1CaloEGammaEmulator,
)

process.legacyTPTable = ecalEBTable.clone(
    src=legacy_source,
    name="LegacyEcalEB",
)
process.phase2TPTable = ecalEBPhase2Table.clone(
    src=phase2_source,
    name="Phase2EcalEB",
)

process.legacyClusterProducer = l1tPhase2L1CaloEGammaEmulator.clone(
    ecalTPEB=legacy_source,
    hcalTP=cms.InputTag("simHcalTriggerPrimitiveDigis", "", "DIGI2RAW"),
    usePhase2TPs=cms.bool(False),
    phase2Sample=cms.uint32(0),
    phase2EtLSB=cms.double(0.488),
)
process.phase2ClusterProducer = l1tPhase2L1CaloEGammaEmulator.clone(
    ecalTPEB=phase2_source,
    hcalTP=cms.InputTag("simHcalTriggerPrimitiveDigis", "", "DIGI2RAW"),
    usePhase2TPs=cms.bool(True),
    phase2Sample=cms.uint32(0),
    phase2EtLSB=cms.double(0.488),
)

process.legacyRCTClusterTable = ecalRCTClusterTable.clone(
    src=cms.InputTag("legacyClusterProducer", "RCTClusters"),
    name="LegacyEcalRCTClusters",
)
process.legacyGCTClusterTable = ecalGCTClusterTable.clone(
    src=cms.InputTag("legacyClusterProducer", "GCTClusters"),
    name="LegacyEcalGCTClusters",
)
process.phase2RCTClusterTable = ecalRCTClusterTable.clone(
    src=cms.InputTag("phase2ClusterProducer", "RCTClusters"),
    name="Phase2EcalRCTClusters",
)
process.phase2GCTClusterTable = ecalGCTClusterTable.clone(
    src=cms.InputTag("phase2ClusterProducer", "GCTClusters"),
    name="Phase2EcalGCTClusters",
)

from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.NanoAOD.simpleGenParticleFlatTableProducer_cfi import simpleGenParticleFlatTableProducer

process.genElectronTable = simpleGenParticleFlatTableProducer.clone(
    src=cms.InputTag("genParticles", "", options.genProcess),
    cut=cms.string("abs(pdgId) == 11 && status == 1"),
    name=cms.string("GenElectron"),
    doc=cms.string("Stable generator electrons and positrons"),
    variables=cms.PSet(
        pt=Var("pt", float, doc="Generator transverse momentum [GeV]", precision=-1),
        eta=Var("eta", float, doc="Generator pseudorapidity", precision=-1),
        phi=Var("phi", float, doc="Generator azimuth [rad]", precision=-1),
        pdgId=Var("pdgId", int, doc="PDG identifier"),
    ),
)

process.nano_step = cms.Path(
    process.legacyTPTable
    + process.phase2TPTable
    + process.legacyClusterProducer
    + process.phase2ClusterProducer
    + process.legacyRCTClusterTable
    + process.legacyGCTClusterTable
    + process.phase2RCTClusterTable
    + process.phase2GCTClusterTable
    + process.genElectronTable
)
process.NANOAODoutput = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string(str(output)),
    outputCommands=cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*_*_TPNANO"),
)
process.endjob_step = cms.EndPath(process.NANOAODoutput)
