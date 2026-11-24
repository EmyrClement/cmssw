"""Regenerate legacy ECAL digis and TPs from no-pileup D128 EB SIM hits.

The output preserves the stored Phase-2 ECAL TPs, HCAL TPs and generator
particles, so the legacy and Phase-2 Nano jobs can be compared event by event.
"""

from pathlib import Path

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9


options = VarParsing()
options.register("outputFile", "ecal_tp_matched_input.root", options.multiplicity.singleton,
                 options.varType.string, "Matched EDM output file")
options.register("maxEvents", 100, options.multiplicity.singleton,
                 options.varType.int, "Number of events")
options.register("inputFiles", [], options.multiplicity.list,
                 options.varType.string, "No-pileup EDM files containing EB SIM hits")
options.register(
    "seed",
    12345,
    options.multiplicity.singleton,
    options.varType.int,
    "Random seed for no-pileup ECAL digitisation",
)
options.parseArguments()

if not options.inputFiles:
    raise ValueError("Specify no-pileup inputFiles containing g4SimHits:EcalHitsEB:SIM")
output = Path(options.outputFile.removeprefix("file:"))
if output.exists():
    raise ValueError(f"Refusing to overwrite {output}")

# Deliberately omit the ECAL development modifiers for the legacy electronics.
process = cms.Process("LEGACY", Phase2C26I13M9)
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D128Reco_cff")
process.CaloGeometryBuilder.SelectedCalos = cms.vstring("EcalBarrel")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, "150X_mcRun4_realistic_v1", "")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.mix.digitizers = cms.PSet(ecal=process.mix.digitizers.ecal.clone())
process.mix.mixObjects = cms.PSet()
process.RandomNumberGeneratorService.mix.initialSeed = options.seed
for label in list(process.aliases_()):
    delattr(process, label)
process.simEcalUnsuppressedDigis = cms.EDAlias(
    mix=cms.VPSet(cms.PSet(type=cms.string("EBDigiCollection")))
)
process.load("SimCalorimetry.EcalEBTrigPrimProducers.ecalEBTriggerPrimitiveDigis_cfi")
process.simEcalEBTriggerPrimitiveDigis.barrelEcalDigis = cms.InputTag(
    "simEcalUnsuppressedDigis", "", "LEGACY"
)

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
    wantSummary=cms.untracked.bool(True),
)
process.legacy_step = cms.Path(process.mix + process.simEcalEBTriggerPrimitiveDigis)
process.EDMoutput = cms.OutputModule(
    "PoolOutputModule",
    fileName=cms.untracked.string(str(output)),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep *_g4SimHits_EcalHitsEB_SIM",
        "keep *_generatorSmeared_*_SIM",
        "keep *_genParticles_*_DIGI2RAW",
        "keep *_simEcalUnsuppressedDigis_*_DIGI2RAW",
        "keep *_simEcalEBTriggerPrimitivePhase2Digis_*_DIGI2RAW",
        "keep *_simHcalTriggerPrimitiveDigis_*_DIGI2RAW",
        "keep *_simEcalUnsuppressedDigis_*_LEGACY",
        "keep *_simEcalEBTriggerPrimitiveDigis_*_LEGACY",
    ),
)
process.endjob_step = cms.EndPath(process.EDMoutput)
