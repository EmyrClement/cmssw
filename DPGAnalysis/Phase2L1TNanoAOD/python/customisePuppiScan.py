"""cmsDriver customisation for the optional barrel offline-PUPPI scan inputs."""

import FWCore.ParameterSet.Config as cms
from DPGAnalysis.Phase2L1TNanoAOD.puppiScan_cff import add_puppi_scan_tables


def customisePuppiScanNano(process):
    """Attach PF/PUPPI tables and enable diagnostic maps without changing the algorithm.

    The first offline prototype supports LinearizedPuppi, one emulated input
    vertex and cut-based charged association. The diagnostic maps retain the
    actual alpha sum, region, input slot and charged association used by the
    configured emulation.
    """
    for required in ("l1tLayer1Barrel", "l1tLayer2Deregionizer", "nanoAOD_step", "NANOAODoutput"):
        if not hasattr(process, required):
            raise ValueError("customisePuppiScanNano requires process." + required)

    barrel = process.l1tLayer1Barrel
    if not getattr(barrel, "vtxCollectionEmulation", cms.bool(True)).value():
        raise ValueError("customisePuppiScanNano requires emulated VertexWord inputs")
    if barrel.puAlgo.value() != "LinearizedPuppi":
        raise ValueError("customisePuppiScanNano requires LinearizedPuppi")
    if barrel.nVtx.value() != 1 or barrel.puAlgoParameters.nVtx.value() != 1:
        raise ValueError("customisePuppiScanNano requires one primary vertex")
    if barrel.puAlgoParameters.useMLAssociation.value():
        raise ValueError("customisePuppiScanNano requires cut-based charged association")

    if getattr(barrel.puAlgoParameters, "fakePuppi", cms.bool(False)).value():
        raise ValueError("customisePuppiScanNano requires fakePuppi=False")
    if barrel.puAlgoParameters.nOut.value() != barrel.puAlgoParameters.nIn.value():
        raise ValueError("customisePuppiScanNano requires nOut == nIn for the offline prototype")

    process = add_puppi_scan_tables(process, path="nanoAOD_step", output="NANOAODoutput")
    barrel.producePuppiDiagnostics = cms.bool(True)
    process.puppiScanVertexTable.src = barrel.vtxCollection
    return process
