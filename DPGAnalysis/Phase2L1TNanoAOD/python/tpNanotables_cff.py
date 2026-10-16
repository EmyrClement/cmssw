"""Standard TP input Nano tables, using the menu Nano configuration style.

Defaults resolve the latest matching products
"""
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var


# Legacy ECAL TPs
ecalEBTable = cms.EDProducer(
    'SimpleEcalEBTPFlatTableProducer',
    src = cms.InputTag("simEcalEBTriggerPrimitiveDigis",""),
    cut = cms.string(''),
    name = cms.string('EcalEB'),
    doc = cms.string('One row per ECAL barrel TP; object-defined sample of interest; encoded units'),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        rawId = Var('id().rawId()', 'uint', doc = 'Raw detector ID'),
        ieta = Var('id().ieta()', int, doc = 'Detector eta index'),
        iphi = Var('id().iphi()', int, doc = 'Detector phi index'),
        nSamples = Var('size()', int, doc = 'Number of stored samples'),
        sampleIndex = Var('sampleOfInterest()', int, doc = 'Selected sample index; -1 if unavailable'),
        encodedEt = Var('encodedEt()', int, doc = 'Selected energy code; -1 if unavailable'),
        time = Var('time()', int, doc = 'Selected time code; -1 if unavailable'),
        spike = Var('?sampleOfInterest() >= 0 ? l1aSpike() : -1', int, doc = 'Selected spike flag; -1 if unavailable'), # Explicit guard, as a bool accessor can return -1
        raw = Var('?sampleOfInterest() >= 0 ? encodedEt() + 1024*l1aSpike() + 2048*time() : -1', int, doc = 'Selected sample raw accessor word; -1 if unavailable'), # Reconstruct the exact 16-bit raw word from its 10+1+5 bit fields instead.
    ),
)

# Latest Phase-2 ECAL TPs, with a different sample selection and a different raw word layout.
# Not included in the default task, as not produced by default
ecalEBPhase2Table = cms.EDProducer(
    'SimpleEcalEBPhase2TPFlatTableProducer',
    src = cms.InputTag("simEcalEBTriggerPrimitivePhase2Digis",""),
    cut = cms.string(''),
    name = cms.string('EcalEBPhase2'),
    doc = cms.string('One row per ECAL barrel TP; explicitly selected sample; encoded units'),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        rawId = Var('id().rawId()', 'uint', doc = 'Raw detector ID'),
        ieta = Var('id().ieta()', int, doc = 'Detector eta index'),
        iphi = Var('id().iphi()', int, doc = 'Detector phi index'),
        nSamples = Var('size()', int, doc = 'Number of stored samples'),
        sampleIndex = Var('?(0 >= 0 && 0 < size()) ? 0 : -1', int, doc = 'Selected sample index; -1 if unavailable'),
        encodedEt = Var('?(0 >= 0 && 0 < size()) ? sample(0).encodedEt() : -1', int, doc = 'Selected energy code; -1 if unavailable'),
        time = Var('?(0 >= 0 && 0 < size()) ? sample(0).time() : -1', int, doc = 'Selected time code; -1 if unavailable'),
        spike = Var('?(0 >= 0 && 0 < size()) ? sample(0).l1aSpike() : -1', int, doc = 'Selected spike flag; -1 if unavailable'),
        raw = Var('?(0 >= 0 && 0 < size()) ? sample(0).raw() : -1', int, doc = 'Selected sample raw accessor word; -1 if unavailable'),
    ),
)

# HCAL inputs, barrel
hcalBarrelTable = cms.EDProducer(
    'SimpleHcalTPFlatTableProducer',
    src = cms.InputTag("simHcalTriggerPrimitiveDigis",""),
    cut = cms.string('abs(id().ieta()) <= 16'),
    name = cms.string('HcalBarrel'),
    doc = cms.string('Detector input; raw accessor values; no quality/energy cut'),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        rawId = Var('id().rawId()', 'uint', doc = 'Raw detector ID'),
        ieta = Var('id().ieta()', int, doc = 'Detector eta index'),
        iphi = Var('id().iphi()', int, doc = 'Detector phi index'),
        version = Var('id().version()', int, doc = 'HCAL trigger-tower ID version'),
        depth = Var('id().depth()', int, doc = 'HCAL trigger-tower depth'),
        nSamples = Var('size()', int, doc = 'Number of stored samples'),
        presamples = Var('presamples()', int, doc = 'Number of samples before t0; index of the selected sample'),
        compressedEt = Var('t0().compressedEt()', int, doc = 'Raw energy code; no GeV conversion'),
        raw = Var('t0().raw()', int, doc = 'Raw t0 sample word'),
        fineGrain0 = Var('t0().fineGrain(0)', int, doc = 'Fine-grain bit 0'),
        fineGrain1 = Var('t0().fineGrain(1)', int, doc = 'Fine-grain bit 1'),
        fineGrain2 = Var('t0().fineGrain(2)', int, doc = 'Fine-grain bit 2'),
        fineGrain3 = Var('t0().fineGrain(3)', int, doc = 'Fine-grain bit 3'),
        fineGrain4 = Var('t0().fineGrain(4)', int, doc = 'Fine-grain bit 4'),
        fineGrain5 = Var('t0().fineGrain(5)', int, doc = 'Fine-grain bit 5'),
    ),
)

# Same HCAL fields, with the HF ring selection.
hfTable = hcalBarrelTable.clone(
    name = cms.string("HF"),
    cut = cms.string("abs(id().ieta()) >= 29"),
)

# Track inputs
trackTable = cms.EDProducer(
    'SimpleTPTrackFlatTableProducer',
    src = cms.InputTag("l1tTTTracksFromTrackletEmulation","Level1TTTracks"),
    cut = cms.string(''),
    name = cms.string('Tracks'),
    doc = cms.string('Detector input; raw accessor values; no quality/energy cut'),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        pt = Var('momentum().perp()', float, doc = 'Transverse momentum', precision = -1),
        eta = Var('momentum().eta()', float, doc = 'Pseudorapidity', precision = -1),
        phi = Var('phi()', float, doc = 'Azimuthal angle', precision = -1),
        z0 = Var('z0()', float, doc = 'Track longitudinal impact parameter', precision = -1),
        d0 = Var('d0()', float, doc = 'Track transverse impact parameter', precision = -1),
        rInv = Var('rInv()', float, doc = 'Signed track curvature', precision = -1),
        tanL = Var('tanL()', float, doc = 'Track tan(lambda)', precision = -1),
        chi2 = Var('chi2()', float, doc = 'Track fit chi2', precision = -1),
        chi2Red = Var('chi2Red()', float, doc = 'Reduced track fit chi2', precision = -1),
        chi2BendRed = Var('chi2BendRed()', float, doc = 'Reduced stub-bend chi2', precision = -1),
        trkMVA1 = Var('trkMVA1()', float, doc = 'Track quality MVA1', precision = -1),
        nStubsWord = Var('getNStubs()', int, doc = 'Set bits in packed track-word hit pattern; not reference-vector length'),
        hitPattern = Var('hitPattern()', 'uint', doc = 'Track fit hit pattern'),
        nFitPars = Var('nFitPars()', int, doc = 'Number of track fit parameters'),
    ),
)

# HGCal inputs, clusters
hgcalClusterTable = cms.EDProducer(
    'SimpleTPHGCalClusterFlatTableProducer',
    src = cms.InputTag("l1tHGCalBackEndLayer2Producer","HGCalBackendLayer2Processor3DClustering"),
    cut = cms.string(''),
    name = cms.string('HGCalClusters'),
    doc = cms.string('Detector input; raw accessor values; no quality/energy cut'),
    extension = cms.bool(False),
    minBX = cms.int32(0),
    maxBX = cms.int32(0),
    alwaysWriteBXValue = cms.bool(True),
    variables = cms.PSet(
        pt = Var('pt()', float, doc = 'Transverse momentum', precision = -1),
        eta = Var('eta()', float, doc = 'Pseudorapidity', precision = -1),
        phi = Var('phi()', float, doc = 'Azimuthal angle', precision = -1),
        energy = Var('energy()', float, doc = 'Energy', precision = -1),
        hwPt = Var('hwPt()', int, doc = 'Hardware transverse momentum code'),
        hwEta = Var('hwEta()', int, doc = 'Hardware eta code'),
        hwPhi = Var('hwPhi()', int, doc = 'Hardware phi code'),
        hwQual = Var('hwQual()', int, doc = 'Hardware quality code'),
        hOverE = Var('hOverE()', float, doc = 'Hadronic over electromagnetic energy', precision = -1),
        sigmaEtaEtaTot = Var('sigmaEtaEtaTot()', float, doc = 'Total eta shower width', precision = -1),
        sigmaPhiPhiTot = Var('sigmaPhiPhiTot()', float, doc = 'Total phi shower width', precision = -1),
        sigmaZZ = Var('sigmaZZ()', float, doc = 'Longitudinal shower width', precision = -1),
        sigmaRRTot = Var('sigmaRRTot()', float, doc = 'Total radial shower width', precision = -1),
        zBarycenter = Var('zBarycenter()', float, doc = 'Shower barycenter z coordinate', precision = -1),
        nClusters = Var('size()', int, doc = 'Number of constituent clusters'),
        showerLength = Var('showerLength()', int, doc = 'Shower length in layers'),
        firstLayer = Var('firstLayer()', int, doc = 'First shower layer'),
    ),
)

# HGCal inputs, towers
hgcalTowerTable = cms.EDProducer(
    'SimpleTPHGCalTowerFlatTableProducer',
    src = cms.InputTag("l1tHGCalTowerProducer","HGCalTowerProcessor"),
    cut = cms.string(''),
    name = cms.string('HGCalTowers'),
    doc = cms.string('Detector input; raw accessor values; no quality/energy cut'),
    extension = cms.bool(False),
    minBX = cms.int32(0),
    maxBX = cms.int32(0),
    alwaysWriteBXValue = cms.bool(True),
    variables = cms.PSet(
        pt = Var('pt()', float, doc = 'Transverse momentum', precision = -1),
        eta = Var('eta()', float, doc = 'Pseudorapidity', precision = -1),
        phi = Var('phi()', float, doc = 'Azimuthal angle', precision = -1),
        energy = Var('energy()', float, doc = 'Energy', precision = -1),
        hwPt = Var('hwPt()', int, doc = 'Hardware transverse momentum code'),
        hwEta = Var('hwEta()', int, doc = 'Hardware eta code'),
        hwPhi = Var('hwPhi()', int, doc = 'Hardware phi code'),
        hwQual = Var('hwQual()', int, doc = 'Hardware quality code'),
        etEm = Var('etEm()', float, doc = 'Electromagnetic transverse energy', precision = -1),
        etHad = Var('etHad()', float, doc = 'Hadronic transverse energy', precision = -1),
        rawId = Var('id().rawId()', 'uint', doc = 'Raw detector ID'),
        ieta = Var('id().iEta()', int, doc = 'Detector eta index'),
        iphi = Var('id().iPhi()', int, doc = 'Detector phi index'),
        zside = Var('zside()', int, doc = 'Detector z side'),
        subdet = Var('id().subdet()', int, doc = 'Subdetector identifier'),
        hwEtEm = Var('hwEtEm()', int, doc = 'Hardware electromagnetic transverse energy code'),
        hwEtHad = Var('hwEtHad()', int, doc = 'Hardware hadronic transverse energy code'),
        hwEtRatio = Var('hwEtRatio()', int, doc = 'Hardware energy ratio code'),
    ),
)

# Primary vertex

# Reuse the menu's generic vertex table. Keep all candidates, with no primary cut.
from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nanotables_cff import vtxTable as _menuVtxTable
from DPGAnalysis.Phase2L1TNanoAOD.genVertex_cfi import genVertex as _genVertex

vtxTable = _menuVtxTable

# One event-level generator origin, not a list of PU vertices or an L1 match.
genVertexTable = _genVertex.clone(
    src = cms.InputTag('genParticles', 'xyz0'),
    name = cms.string('GenVtx'),
)

# Phase-2 ECAL is opt-in because standard stored inputs need not contain it.
tpTablesTask = cms.Task(
    ecalEBTable,
    hcalBarrelTable,
    hfTable,
    trackTable,
    hgcalClusterTable,
    hgcalTowerTable,
    vtxTable,
    genVertexTable,
)
