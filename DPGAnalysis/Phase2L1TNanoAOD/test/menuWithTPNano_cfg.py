"""Wrap the existing generated menu config; preserve its upstream workflow."""
import os
import runpy
from pathlib import Path
from FWCore.ParameterSet.VarParsing import VarParsing
from DPGAnalysis.Phase2L1TNanoAOD.customiseTPNano import add_tp_tables, add_process_tp_tables

options = VarParsing('analysis')
options.register('menuConfig', str(Path(os.environ['CMSSW_BASE'])/'src/v45_rerunL1wTT_cfg.py'),
                 options.multiplicity.singleton, options.varType.string, 'Generated menu cfg')
options.register('includeHLT', False, options.multiplicity.singleton, options.varType.bool,
                 'Add stored HLT tables with HLT suffix')
options.register('sourceProcess', '', options.multiplicity.singleton, options.varType.string,
                 'Add copies pinned to this process, e.g. HLT or GENSIMDIGIRAW')
options.outputFile = 'menu_tp_nano.root'
options.maxEvents = 1
options.parseArguments()
process = runpy.run_path(options.menuConfig)['process']
if options.inputFiles:
    process.source.fileNames = options.inputFiles
process.maxEvents.input = options.maxEvents
process.NANOAODoutput.fileName = options.outputFile
if Path(options.outputFile.removeprefix('file:')).exists():
    raise ValueError('Refusing to overwrite '+options.outputFile)
process = add_tp_tables(process)
if options.includeHLT:
    process = add_process_tp_tables(process)
if options.sourceProcess and not (options.includeHLT and options.sourceProcess == 'HLT'):
    process = add_process_tp_tables(process, options.sourceProcess)
