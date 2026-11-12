"""First milestone: the CI pion workflow with Nano output and optional reference tuple."""
from pathlib import Path
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from FastPUPPI.NtupleProducer.runPerformanceNTuple import process, goGun, noPU, goMT, noResp
from DPGAnalysis.Phase2L1TNanoAOD.pionCalibration_cff import add_pion_calibration_tables

options = VarParsing('analysis')
options._tagOrder = []  # Honour the requested output filename exactly.
options.register('writeReference', False, options.multiplicity.singleton, options.varType.bool,
                 'Also write the original FastPUPPI response tuple for validation')
options.register('correctionFile', '', options.multiplicity.singleton, options.varType.string,
                 'Optional barrel hadronic calibration ROOT file')
options.register('hgcalCorrectionFile', '', options.multiplicity.singleton, options.varType.string,
                 'Optional HGCal hadronic calibration ROOT file')
options.maxEvents = 10
options.outputFile = 'pion_calibration_nano.root'
options.parseArguments()
output_files = [options.outputFile]
if options.writeReference:
    output_files.append(options.outputFile + '.reference.root')
for filename in output_files:
    if Path(filename.removeprefix('file:')).exists():
        raise ValueError('Refusing to overwrite ' + filename)
if not options.inputFiles:
    raise ValueError('Supply inputFiles=file:/path/to/pion_EDM.root or an XRootD URL')

# Same production configuration as the reference CI, without source-text edits.
goGun(calib=0)
noPU()
goMT()
process.l1tLayer1BarrelRaw = process.l1tLayer1Barrel.clone(
    gctHadInputConversionParameters = process.l1tLayer1Barrel.gctHadInputConversionParameters.clone(
        gctHadCorrector = '',
    ),
)
process.l1tLayer1HGCalRaw = process.l1tLayer1HGCal.clone(
    hgcalInputConversionParameters = process.l1tLayer1HGCal.hgcalInputConversionParameters.clone(
        corrector = '',
    ),
)
process.l1tLayer1HGCalNoTKRaw = process.l1tLayer1HGCalNoTK.clone(
    hgcalInputConversionParameters = process.l1tLayer1HGCalNoTK.hgcalInputConversionParameters.clone(
        corrector = '',
    ),
)
if options.correctionFile:
    process.l1tLayer1Barrel.gctHadInputConversionParameters.gctHadCorrector = cms.string(
        options.correctionFile
    )
if options.hgcalCorrectionFile:
    process.l1tLayer1HGCal.hgcalInputConversionParameters.corrector = cms.string(
        options.hgcalCorrectionFile
    )
    process.l1tLayer1HGCalNoTK.hgcalInputConversionParameters.corrector = cms.string(
        options.hgcalCorrectionFile
    )
process.extraPFStuff.add(
    process.l1tLayer1BarrelRaw,
    process.l1tLayer1HGCalRaw,
    process.l1tLayer1HGCalNoTKRaw,
)
process.source.fileNames = options.inputFiles
process.maxEvents.input = options.maxEvents
process.outnano.fileName = options.outputFile
process.end += process.outnano  # goGun's respOnly removed this output module.
process = add_pion_calibration_tables(process)
if options.writeReference:
    process.ntuple.objects = process.calibResponseTable.objects.clone()
    process.TFileService.fileName = options.outputFile + '.reference.root'
else:
    noResp()
    del process.TFileService
