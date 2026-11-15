"""Optional uncorrected SeededCone jet table for PUPPI JEC derivation."""

import FWCore.ParameterSet.Config as cms

from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nanotables_cff import (
    sc4JetTable as _sc4JetTable,
)


sc4RawJetTable = _sc4JetTable.clone(
    src = cms.InputTag("l1tSC4PFL1PuppiEmulator"),
    name = cms.string("L1puppiJetSC4Raw"),
    doc = cms.string(
        "Uncorrected SeededCone 0.4 PUPPI jet for calibration, origin: Correlator"
    ),
)

puppiJetCalibrationTask = cms.Task(sc4RawJetTable)


def add_puppi_jet_calibration_table(
    process, path = "nanoAOD_step", output = "NANOAODoutput"
):
    """Attach the uncorrected jet table and keep it in the selected output."""
    process.load("DPGAnalysis.Phase2L1TNanoAOD.puppiJetCalibration_cff")
    getattr(process, path).associate(process.puppiJetCalibrationTask)
    getattr(process, output).outputCommands.append(
        "keep nanoaodFlatTable_sc4RawJetTable_*_" + process.name_()
    )
    return process
