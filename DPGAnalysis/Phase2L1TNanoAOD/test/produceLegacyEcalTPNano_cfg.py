"""Produce ECAL TP, RCT/GCT cluster and generator tables for legacy TPs.

The input must contain legacy TPs regenerated from the same SIM hits as the
Phase-2 reference. Use regenerateLegacyEcalTPs_cfg.py to make that input.
"""

from pathlib import Path

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9


options = VarParsing()
options.register("outputFile", "legacy_ecal_tp_nano.root", options.multiplicity.singleton,
                 options.varType.string, "Nano output file")
options.register("maxEvents", 100, options.multiplicity.singleton,
                 options.varType.int, "Number of events")
options.register("inputFiles", [], options.multiplicity.list,
                 options.varType.string, "Input EDM files")
options.register(
    "tpProcess",
    "LEGACY",
    options.multiplicity.singleton,
    options.varType.string,
    "Process containing simEcalEBTriggerPrimitiveDigis",
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
    raise ValueError("Specify inputFiles containing regenerated legacy ECAL TPs")
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

tp_source = cms.InputTag("simEcalEBTriggerPrimitiveDigis", "", options.tpProcess)
from DPGAnalysis.Phase2L1TNanoAOD.tpNanotables_cff import (
    ecalEBTable,
    ecalGCTClusterTable,
    ecalRCTClusterTable,
)
from L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi import (
    l1tPhase2L1CaloEGammaEmulator,
)

process.ecalTPTable = ecalEBTable.clone(src=tp_source, name="EcalEB")
process.ecalTPTable.doc = cms.string("One row per ECAL barrel TP; selected sample; encoded units")
process.ecalClusterProducer = l1tPhase2L1CaloEGammaEmulator.clone(
    ecalTPEB=tp_source,
    hcalTP=cms.InputTag("simHcalTriggerPrimitiveDigis", "", "DIGI2RAW"),
    usePhase2TPs=cms.bool(False),
    phase2Sample=cms.uint32(0),
    phase2EtLSB=cms.double(0.488),
)
process.ecalRCTClusterTable = ecalRCTClusterTable.clone(
    src=cms.InputTag("ecalClusterProducer", "RCTClusters"),
)
process.ecalGCTClusterTable = ecalGCTClusterTable.clone(
    src=cms.InputTag("ecalClusterProducer", "GCTClusters"),
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
    process.ecalTPTable
    + process.ecalClusterProducer
    + process.ecalRCTClusterTable
    + process.ecalGCTClusterTable
    + process.genElectronTable
)
process.NANOAODoutput = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string(str(output)),
    outputCommands=cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*_*_TPNANO"),
)
process.endjob_step = cms.EndPath(process.NANOAODoutput)
