"""cmsDriver customisation for the optional PUPPI jet calibration table."""

from DPGAnalysis.Phase2L1TNanoAOD.puppiJetCalibration_cff import (
    add_puppi_jet_calibration_table,
)


def customisePuppiJetCalibrationNano(process):
    """Add the uncorrected SeededCone jet table to a Phase-2 L1 Nano workflow."""
    required_modules = ("nanoAOD_step", "NANOAODoutput")
    for required in required_modules:
        if not hasattr(process, required):
            raise ValueError(
                "customisePuppiJetCalibrationNano requires process." + required
            )

    return add_puppi_jet_calibration_table(process)
