"""Event-level generator origin, using the existing generic Nano singleton producer."""
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.NanoAOD.simpleXYZPointFlatTableProducer_cfi import simpleXYZPointFlatTableProducer

genVertex = simpleXYZPointFlatTableProducer.clone(
    src=cms.InputTag('genParticles', 'xyz0'),
    name=cms.string('GenVtx'),
    doc=cms.string('Generator origin from genParticles:xyz0; event reference, not a nearest-L1 match'),
    variables=cms.PSet(
        x=Var('X', float, doc='Generator origin x [cm]', precision=-1),
        y=Var('Y', float, doc='Generator origin y [cm]', precision=-1),
        z=Var('Z', float, doc='Generator origin z [cm]', precision=-1),
    ),
)
