# TP input tables in the menu Nano

The standard definitions live in `python/tpNanotables_cff.py`. They include
legacy ECAL barrel, HCAL barrel/HF, tracks, HGCal clusters/towers, L1 vertices
the GEN origin, and the three ECAL cluster collections. Phase-2 ECAL is an optional table. The task extracts existing
products; it does not schedule TP reconstruction. Default InputTags resolve the
latest matching product.

After constructing the menu Nano configuration, load and associate the task
explicitly:

```python
process.load("DPGAnalysis.Phase2L1TNanoAOD.tpNanotables_cff")
process.nanoAOD_step.associate(process.tpTablesTask)
```

The standard definitions intentionally omit process names and therefore resolve
the latest matching product. For a stored/reproduced comparison, clone only the
required tables in the study configuration, give them distinct names, and pin
their `InputTag` process explicitly. The two ECAL examples below show this
pattern. Other jobs must have a scheduled Nano path and an output module that
keeps the produced `nanoaodFlatTable` products.

## Menu-based example

In the chosen release's src directory with its runtime active:

```sh
cmsRun DPGAnalysis/Phase2L1TNanoAOD/test/menuWithTPNano_cfg.py \
  menuConfig=/absolute/path/v45_rerunL1wTT_cfg.py \
  inputFiles=file:/absolute/path/input.root \
  outputFile=menu_tp_nano.root maxEvents=10
```

The example wraps your cmsDriver-generated menu configuration, preserving its
reconstruction, geometry, conditions and customisations. With no `menuConfig`,
it looks for `v45_rerunL1wTT_cfg.py` in the release src directory; the generated
file is not shipped by this package. Choose a config compatible with the input
and release. The wrapper directly associates the standard task with the menu
Nano path.

These basics have no TPStudies dependency. Optional study provenance recording,
expert extraction adapters and comparisons remain in TPStudies. Save the expanded
configuration, input identities and execution logs for comparisons.

## Standard ECAL cluster outputs

`ecalClusterTablesTask` contains `ecalRCTClusterTable`, `ecalGCTClusterTable`
and `ecalEmulatorClusterTable`, with distinct source labels and table names.
It is included in `tpTablesTask` automatically.
All three input collections must be available when using the standard task.
These tables serialize existing CaloCrystalCluster collections; they do not
schedule clustering or enable support for the new ECAL TP class.
Producer-specific placeholder fields are preserved; pt/calibratedPt are not
an uncalibrated/calibrated pair. See column docs for interpretation.

## Matched legacy and Phase-2 ECAL TP comparison

The following standalone examples target the centrally produced no-pileup D128
single-electron RelVal:

```text
/RelValSingleEFlatPt2To100/CMSSW_20_0_0_patch2-150X_mcRun4_realistic_v1-v1/GEN-SIM-DIGI-RAW
```

The central file contains Phase-2 ECAL TPs but no legacy TP collection. First
regenerate the legacy ECAL digis and TPs from the stored EB simulation hits:

```sh
cmsRun DPGAnalysis/Phase2L1TNanoAOD/test/regenerateLegacyEcalTPs_cfg.py \
  inputFiles=root://cms-xrd-global.cern.ch//store/relval/CMSSW_20_0_0_patch2/RelValSingleEFlatPt2To100/GEN-SIM-DIGI-RAW/150X_mcRun4_realistic_v1-v1/2830000/0a1311b3-ccff-447c-86cd-0631aae15537.root \
  outputFile=file:ecal_tp_matched_input.root \
  maxEvents=100
```

The output retains the original Phase-2 ECAL and HCAL TP collections under
process `DIGI2RAW` and adds the same-event legacy ECAL collection under process
`LEGACY`. The recipe is for no-pileup input; it does not reproduce pileup mixing.

The two Nano jobs intentionally write the same table names and schemas into
separate files, making direct branch-by-branch comparison straightforward:

```sh
cmsRun DPGAnalysis/Phase2L1TNanoAOD/test/producePhase2EcalTPNano_cfg.py \
  inputFiles=file:ecal_tp_matched_input.root \
  outputFile=file:phase2_ecal_tp_nano.root \
  maxEvents=100

cmsRun DPGAnalysis/Phase2L1TNanoAOD/test/produceLegacyEcalTPNano_cfg.py \
  inputFiles=file:ecal_tp_matched_input.root \
  outputFile=file:legacy_ecal_tp_nano.root \
  maxEvents=100
```

Alternatively, schedule both configured RCT/GCT producer instances and write
all variants into one Nano:

```sh
cmsRun DPGAnalysis/Phase2L1TNanoAOD/test/produceCombinedEcalTPNano_cfg.py \
  inputFiles=file:ecal_tp_matched_input.root \
  outputFile=file:combined_ecal_tp_nano.root \
  maxEvents=100
```

The combined file uses distinct `LegacyEcal*` and `Phase2Ecal*` table names.
The separate files deliberately use identical `Ecal*` names and schemas.

Each Nano contains the selected ECAL crystal TP fields, stable generator
electrons/positrons, and RCT and GCT clusters. Both cluster variants consume the
same stored `simHcalTriggerPrimitiveDigis::DIGI2RAW` collection. Legacy ECAL
energy uses 0.125 GeV/count through the digi's object-defined sample of interest.
Phase-2 ECAL energy uses explicit sample 0 and 0.488 GeV/count. Both are converted
to the existing emulator's 0.5 GeV internal representation before the unchanged
RCT/GCT algorithm runs.

The separate examples keep input identity and provenance obvious and also permit
the Phase-2 Nano to run directly on the central file. The combined example is
convenient for same-event comparisons. Use the explicit process-name overrides
for samples whose process names differ from these defaults.
