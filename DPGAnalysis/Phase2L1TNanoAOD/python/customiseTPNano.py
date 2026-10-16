"""Small customisations for the standard menu Nano workflow."""
import FWCore.ParameterSet.Config as cms


def add_tp_tables(process, *, path='nanoAOD_step'):
    """Attach the standard task; the menu output already keeps *Table modules."""
    if not hasattr(process, 'tpTablesTask'):
        process.load('DPGAnalysis.Phase2L1TNanoAOD.tpNanotables_cff')
    getattr(process, path).associate(process.tpTablesTask)
    return process


def add_process_tp_tables(process, source_process='HLT', *, path='nanoAOD_step', output='NANOAODoutput'):
    """Copy the default task, pinning inputs and suffixing names with source_process."""
    task_name = 'tpTables'+source_process+'Task'
    if not source_process.isalnum() or hasattr(process, task_name):
        raise ValueError('Use a nonempty alphanumeric process name, once per process')
    setattr(process, task_name, cms.Task())
    task = getattr(process, task_name)
    for label in sorted(process.tpTablesTask.moduleNames()):
        original = getattr(process, label)
        clone = original.clone(
            src=cms.InputTag(original.src.getModuleLabel(), original.src.getProductInstanceLabel(), source_process),
            name=original.name.value()+source_process,
        )
        setattr(process, label+source_process, clone)
        task.add(getattr(process, label+source_process))
    getattr(process, path).associate(task)
    # Copied module labels no longer end in Table, so explicitly keep them.
    getattr(process, output).outputCommands.append(
        'keep nanoaodFlatTable_*Table'+source_process+'_*_'+process.name_())
    return process


def add_hlt_tp_tables(process, source_process='HLT', **kwargs):
    """Compatibility spelling; prefer add_process_tp_tables for other processes."""
    return add_process_tp_tables(process, source_process, **kwargs)


def customiseTPNano(process):
    return add_tp_tables(process)


def customiseTPNanoWithHLT(process):
    return add_hlt_tp_tables(add_tp_tables(process))
