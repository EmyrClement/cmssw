# Phase2-L1Nano
NanoAOD ntupler for Phase-2 L1 Objects

Initially an independent package here: https://github.com/cms-l1-dpg/Phase2-L1Nano/

## Setup

For more information on the latest L1T Phase 2 software developments in CMSSW see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#Development

Corresponding menu twiki section: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PhaseIIL1TriggerMenuTools#Phase_2_L1_Trigger_objects_based


## Usage

### Via cmsDriver

One can append the L1Nano output to the `cmsDriver` command via the `NANO:@Phase2L1DPGwithGen` autoNANO handle, e.g.:
```bash
cmsDriver.py -s L1,L1TrackTrigger,L1P2GT,NANO:@Phase2L1DPGwithGen
```

Check `PhysicsTools/NanoAOD/python/autoNANO.py` for the way this command is defined.

Note that the step key `Phase2L1DPG` does not include the generator and reco-level objects used for MenuTools studies in the nano!
It is mostly created for workflow tests.

An example `cmsDriver` command for 14x files:
```bash
cmsDriver.py -s L1,L1TrackTrigger,L1P2GT,NANO:@Phase2L1DPGwithGen \
--conditions auto:phase2_realistic_T35 \
--geometry ExtendedRun4D110 \
--era Phase2C17I13M9 \
--eventcontent NANOAOD \
--datatier GEN-SIM-DIGI-RAW-MINIAOD \
--customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2TTOn110.customisePhase2TTOn110 \
--filein /store/mc/Phase2Spring24DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_AllTP_140X_mcRun4_realistic_v4-v1/2560000/11d1f6f0-5f03-421e-90c7-b5815197fc85.root \
--fileout file:output_Phase2_L1T.root \
--python_filename rerunL1_cfg.py \
--inputCommands="keep *, drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT" \
--mc \
-n 10 --nThreads 4 --no_exec
```

### Workflows

Two upgrade workflows are implemented to test/run this nano (`Phase2L1DPG`) after the full DigiTrigger chain after the complete L1:
* `.781` - produces NANO in addition to FEVTDEBUG
* `.782` - produces only NANO

And can be executed as e.g. `runTheMatrix.py --what upgrade -i all --ibeos -l 29634.782`. See [here the readme of `runTheMatrix`](https://github.com/cms-sw/cmssw/tree/master/Configuration/PyReleaseValidation/scripts#interactive-runthematrix-shell) and [here a list of workflows](https://github.com/cms-sw/cmssw/tree/master/Configuration/PyReleaseValidation).

Note that if tou want to include the Gen/Offline references, you need to change the import function in the config after the

## Output

The output file is a nanoAOD file with the output branches in the `Events` tree.

An overview of the corresponding content is shown here: https://alobanov.web.cern.ch/L1T/Phase2/L1Nano/l1menu_nano_V38_1400pre3V9_doc_report.html

Size report: https://alobanov.web.cern.ch/L1T/Phase2/L1Nano/l1menu_nano_V38_1400pre3V9_size_report.html

Example:

```python
'run',
'luminosityBlock',
'event',
'L1tkPhoton_saId',
'L1tkPhoton_hwEta',
'L1tkPhoton_hwIso',
'L1tkPhoton_hwPhi',
```

This can be easily handled with [`uproot/awkward`](https://gitlab.cern.ch/cms-podas23/dpg/trigger-exercise/-/blob/solutions/1_Intro_NanoAwk_Analysis_Solution.ipynb) like this:

```python
f = uproot.open("l1nano.root")
events = f["Events"].arrays() 
```

### P2GT emulator decisions
The GT emulator decisions are stored like this now:
```
L1_pSingleTkMuon22_final # seed name and final for post prescale value
```

## Optional calibration tables

`pionCalibration_cff.py` defines an opt-in task with raw/corrected decoded barrel
and HGCal objects, a prompt-pion diagnostic table, and matched calibration-
particle response quantities. The configured particle IDs are electron, photon, neutral and
charged pion, and K-long. A charged-pion fitter can select its subset through
the stored `Response_mc_id` column.
It does not change the standard menu task. For a cmsDriver workflow containing
`NANO:@Phase2L1DPGwithGen`, append this function to `--customise`:

```text
DPGAnalysis/Phase2L1TNanoAOD/customisePionCalibration.customisePionCalibrationNano
```

The customisation clones `l1tLayer1Barrel`, `l1tLayer1HGCal`, and
`l1tLayer1HGCalNoTK` with their hadronic corrections disabled, adds the clones
and calibration tables to `nanoAOD_step`, and keeps their tables through
`NANOAODoutput`. InputTags intentionally omit the process name and hence select
the latest matching products, normally from `L1TReRun` in a rerun job.
The calibration fitter needs the `Response` table; the standard TP tables are an
independent optional customisation.

For a non-cmsDriver job, attach only the table task to existing compatible
producers with:

```python
from DPGAnalysis.Phase2L1TNanoAOD.pionCalibration_cff import add_pion_calibration_tables
process = add_pion_calibration_tables(process, path="p", output="outnano")
```

`L1PFDecodedCaloTableProducer` extends a generic object table with decoded
hardware quantities. It was moved from FastPUPPI, retaining its plugin name and
columns. Use a FastPUPPI revision with the old registration removed when building
both packages together. `PionCalibrationResponseTableProducer` implements the multi-collection
matching, propagation and cone calculations needed by the response fitter. The historical
plugin name is retained for configuration compatibility. Neither producer
includes FastPUPPI headers.

The HGCal calibration columns combine the tracked and no-track regions for each
raw or corrected response quantity. The object tables keep those regions and
their hadronic/EM collections separate. The HGCal fit uses
`Response_L1RawHGCalEM_pt02 / Response_L1RawHGCal_pt02` for its cone EM fraction
and `Response_L1HGCal_pt02` for corrected closure. The production example accepts
`hgcalCorrectionFile` and applies it to both standard HGCal regions.

Keep reusable producers in `plugins/`, workflow-specific opt-in tasks in named
Python fragments, and example jobs in `test/`. Prefer generic Simple tables with
explicit Python variables whenever possible. A separate package is warranted
when workflows develop distinct dependencies or ownership, rather than simply
because the number of table definitions increases.

`test/pionCalibrationNano_cfg.py` is a hackathon production example using the
FastPUPPI emulation configuration; that example still needs FastPUPPI. The
registered `testPionCalibrationTables` test only checks zero-event construction
and does not import FastPUPPI. Event-level agreement is validated separately
with the pion sample and original analyzer.


## Optional offline PUPPI scan tables

Append this opt-in customisation to the full menu cmsDriver command:

```text
DPGAnalysis/Phase2L1TNanoAOD/customisePuppiScan.customisePuppiScanNano
```

It attaches `puppiScanTask` to `nanoAOD_step`, keeps its tables through
`NANOAODoutput`, and enables the barrel producer's `producePuppiDiagnostics`
ValueMaps. The default menu table task and the emulation settings are unchanged.
The supported prototype uses `LinearizedPuppi`, emulated VertexWords, one primary
vertex, cut-based charged association, `fakePuppi=False` and `nOut == nIn`.
Unsupported modes fail during configuration rather than silently producing an
ambiguous offline contract. The generic table fragment may also be attached
with `add_puppi_scan_tables(process, path="p", output="outnano")` after explicitly
enabling compatible diagnostic products and setting the intended vertex source.

| Table | Input | Meaning |
| --- | --- | --- |
| `PuppiScanPF` | `l1tLayer1Barrel:PF` | Fiducial unweighted PF, with no additional Nano cut |
| `PuppiScanRef` | `l1tLayer1Barrel:Puppi` | Weighted barrel reference after PUPPI cuts and regional sorting |
| `PuppiScanFinal` | `l1tLayer2Deregionizer:Puppi` | Full-detector reference after deregionization |
| `PuppiScanVertex` | Actual barrel `vtxCollection` tag | Every input VertexWord in original order |

All float columns retain full float precision (`precision=-1`). PF `eta/phi`
are the global vertex direction for charged candidates; `caloEta/caloPhi`
are their global calorimeter coordinates. PF `hwEta/hwPhi` are region-local
calorimeter codes. PUPPI `hwEta/hwPhi` are global output codes instead. The PUPPI
products do not fill `caloEta/caloPhi`: those columns preserve their constructor
zeros and are **not** measured directions. Neutrals have zero `z0/hwZ0`;
`hwPuppiWeight` preserves the stored integer accessor, which is not the raw
fractional weight bit pattern; use `puppiWeight` for weight comparisons. The `pt`
column is already weighted in the reference tables, so do not weight it again.
The floating `puppiWeight` retains the 1/256 grid and the pT scale is 0.25 GeV/count.

The PF table also reads four producer-owned integer ValueMaps through generic
`ExtVar` columns: `sumRaw` (the actual fixed-point neutral alpha sum), `region`,
`slot` and `association`. `sumRaw=-1` marks charged or ineligible candidates;
`association=-1` marks neutral or ineligible charged candidates. Neutral input slots start at configured
`nTrack`. These maps preserve the inputs needed to scan the weighting parameters
without reconstructing overlapping regional track lists from fiducial PF rows.
They do not authorize changing alpha geometry, track truncation or association
settings while treating the stored sums as unchanged. The reference tables do
not consume these PF-keyed maps.

Vertex `hwZ0` retains the unsigned 15-bit VertexWord pattern (signed two's
complement, 1/512 cm/count), distinct from the signed PF code (0.05 cm/count).
To recover the primary used by Layer 1, iterate in stored order starting with
`ptsum=0`, replacing the chosen vertex whenever `ptsum == 0` or `sumPt > ptsum`.
Positive-sumPt ties retain the first row; all-zero collections select the last.
An empty collection yields z0=0. Then apply the Layer-1 0.05 cm quantisation
with half-way rounding away from zero. No validity or truth-proximity cut is
applied by that selection. Save the expanded configuration and source identity
alongside every Nano; row indices are event-local, not persistent identities.

The registered `testPuppiScanTables` and `testPuppiScanCustomisation` tests cover
zero-event producer construction, configuration guards and table attachment.
They do not establish event-level closure or agreement with a standalone scan.
