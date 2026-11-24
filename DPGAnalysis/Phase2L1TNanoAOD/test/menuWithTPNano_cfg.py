"""Wrap a generated menu config and directly attach the standard TP tables."""

import os
import runpy
from pathlib import Path

from FWCore.ParameterSet.VarParsing import VarParsing


options = VarParsing()
options.register(
    "menuConfig",
    str(Path(os.environ["CMSSW_BASE"]) / "src/v45_rerunL1wTT_cfg.py"),
    options.multiplicity.singleton,
    options.varType.string,
    "Generated menu configuration",
)
options.register(
    "outputFile",
    "menu_tp_nano.root",
    options.multiplicity.singleton,
    options.varType.string,
    "Nano output file",
)
options.register(
    "maxEvents", 1, options.multiplicity.singleton, options.varType.int, "Number of events"
)
options.register(
    "inputFiles", [], options.multiplicity.list, options.varType.string, "Input EDM files"
)
options.parseArguments()

process = runpy.run_path(options.menuConfig)["process"]
if options.inputFiles:
    process.source.fileNames = options.inputFiles
process.maxEvents.input = options.maxEvents
process.NANOAODoutput.fileName = options.outputFile
if Path(options.outputFile.removeprefix("file:")).exists():
    raise ValueError(f"Refusing to overwrite {options.outputFile}")

process.load("DPGAnalysis.Phase2L1TNanoAOD.tpNanotables_cff")
process.nanoAOD_step.associate(process.tpTablesTask)
