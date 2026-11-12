"""cmsDriver customisation for the optional single-particle calibration Nano."""

import FWCore.ParameterSet.Config as cms

from DPGAnalysis.Phase2L1TNanoAOD.pionCalibration_cff import add_pion_calibration_tables


def customisePionCalibrationNano(process):
    """Add raw and corrected calorimeter response tables to a Nano workflow.

    InputTags intentionally have no process name, so they resolve products from
    the current process (normally L1TReRun) before falling back to the input.
    """
    required_modules = (
        "l1tLayer1Barrel",
        "l1tLayer1HGCal",
        "l1tLayer1HGCalNoTK",
        "nanoAOD_step",
        "NANOAODoutput",
    )
    for required in required_modules:
        if not hasattr(process, required):
            raise ValueError(
                "customisePionCalibrationNano requires process." + required
            )

    raw_modules = (
        "l1tLayer1BarrelRaw",
        "l1tLayer1HGCalRaw",
        "l1tLayer1HGCalNoTKRaw",
    )
    for module in raw_modules:
        if hasattr(process, module):
            raise ValueError(
                "customisePionCalibrationNano will not replace process." + module
            )

    process = add_pion_calibration_tables(
        process,
        path="nanoAOD_step",
        output="NANOAODoutput",
    )
    process.l1tLayer1BarrelRaw = process.l1tLayer1Barrel.clone(
        gctHadInputConversionParameters=(
            process.l1tLayer1Barrel.gctHadInputConversionParameters.clone(
                gctHadCorrector=cms.string("")
            )
        )
    )
    process.l1tLayer1HGCalRaw = process.l1tLayer1HGCal.clone(
        hgcalInputConversionParameters=(
            process.l1tLayer1HGCal.hgcalInputConversionParameters.clone(
                corrector=cms.string("")
            )
        )
    )
    process.l1tLayer1HGCalNoTKRaw = process.l1tLayer1HGCalNoTK.clone(
        hgcalInputConversionParameters=(
            process.l1tLayer1HGCalNoTK.hgcalInputConversionParameters.clone(
                corrector=cms.string("")
            )
        )
    )
    process.pionCalibrationTask.add(
        process.l1tLayer1BarrelRaw,
        process.l1tLayer1HGCalRaw,
        process.l1tLayer1HGCalNoTKRaw,
    )
    return process
