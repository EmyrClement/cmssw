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
process.extraPFStuff.add(process.l1tLayer1BarrelRaw)
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
