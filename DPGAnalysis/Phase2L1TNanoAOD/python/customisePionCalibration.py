"""cmsDriver customisation for the optional pion calibration Nano tables."""

import FWCore.ParameterSet.Config as cms

from DPGAnalysis.Phase2L1TNanoAOD.pionCalibration_cff import add_pion_calibration_tables


def customisePionCalibrationNano(process):
    """Add raw/corrected barrel response tables to a cmsDriver Nano workflow.

    InputTags intentionally have no process name, so they resolve products from
    the current process (normally L1TReRun) before falling back to the input.
    """
    for required in ("l1tLayer1Barrel", "nanoAOD_step", "NANOAODoutput"):
        if not hasattr(process, required):
            raise ValueError(
                "customisePionCalibrationNano requires process." + required
            )
    if hasattr(process, "l1tLayer1BarrelRaw"):
        raise ValueError(
            "customisePionCalibrationNano will not replace an existing "
            "process.l1tLayer1BarrelRaw"
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
    process.pionCalibrationTask.add(process.l1tLayer1BarrelRaw)
    return process
