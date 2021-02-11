import logging
import os
import time
import traceback
from collections import defaultdict
from datetime import datetime

import matplotlib
import numpy as np
import pandas as pd
import psutil

matplotlib.use('Agg')
from matplotlib import pyplot as plt


def process_logs(IS_loc):
    """
    Generate the run report for an IS location
    """
    logging.shutdown()
    logloc = os.path.join(IS_loc, 'log/log.log')
    report_run_stats(logloc, most_recent=False, printToo=True, debug=True, save=True)

def report_run_stats(logloc, save=True, most_recent=True, printToo=True,
                    debug=False, plot=False):
    """
    Overall wrapper of this module.

    1) The text log is parsed into a dataframe with "load_log"
    2) A number of reports are generated with "generate_reports"
    3) These reports are saved / displayed here

    """
    if logloc is None:
        return

    # Load the log
    try:
        Ldb = load_log(logloc, most_recent=most_recent)

    except BaseException as e:
        if debug:
            print('Failed to load log file - {1}'.format('None', str(e)))
            traceback.print_exc()
        return

    # Generate reports
    for run, ldb in Ldb.groupby('run_ID'):
        # Make the multi-processing log
        MLdb = load_multiprocessing_log(ldb)

        # Make the reports
        name2report = generate_reports(ldb, MLdb, debug=debug)

        # Print
        if printToo:
            for name, report in name2report.items():
                print("..:: {0} ::..\n{1}".format(name, report))

        if save == True:
            savelocs = [logloc.replace('log.log', '{0}.runtime_summary.txt'.format(run)),
                        logloc.replace('log.log', 'latest.runtime_summary.txt'.format(run))]
            figloc = logloc.replace('log.log', '{0}.ramProfile.png'.format(run))
        elif save == False:
            continue
        else:
            saveloc = save
            figloc = save + '.mem_usage.png'

        for saveloc in savelocs:
            with open(saveloc, 'w') as o:
                for name, report in name2report.items():
                    o.write("..:: {0} ::..\n{1}\n".format(name, report))

    # # Make the plot
    # if plot:
    #     try:
    #         profile_plot(Ldb, saveloc=figloc)
    #     except BaseException as e:
    #         if debug:
    #             print('Failed to make profile plot - {1}'.format('None', str(e)))
    #             traceback.print_exc()

def load_log(logfile, most_recent=True):
    """
    A bunch of messy logic to turn the text file into a dataframe

    The resulting dataframe has 4 columns:
    * log_type = a description of the information that the log contains
    * time = time (in seconds since epoch) that the log was made
    * parseable_string = a string with arbitrary information
    * run_ID = an identifier of a particular instrain run (time when started)

    There are the following types of logs:
    * program_start
    * program_end
    * checkpoint
        - Made using "log_checkpoint". Used for logging regular checkpoints
    * Special_genes
    * Plotting
    * WorkerLog
    * Failure

    """
    table = defaultdict(list)
    with open(logfile) as o:
        prev_line_2 = None
        prev_line = None
        run_ID = None
        for line in o.readlines():
            line = line.strip()
            linewords = [x.strip() for x in line.split()]

            # load new inStrain run
            if 'inStrain version' in line:
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))

                run_ID = datetime.fromtimestamp(epoch_time).strftime('%Y%m%d_%H%M%S')
                cmd = prev_line_2.strip().split('was:')[1].strip()

                table['log_type'].append('program_start')
                table['time'].append(epoch_time)
                table['parsable_string'].append("version={0}; cmd={1}".format(linewords[5], cmd))
                table['run_ID'].append(run_ID)

            # load inStrain run finish
            elif 'inStrain complete' in line:
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))

                table['log_type'].append('program_end')
                table['time'].append(epoch_time)
                table['parsable_string'].append("loglog={0}".format(linewords[14][:-1]))
                table['run_ID'].append(run_ID)

            # regular checkpoints
            elif (len(linewords) > 3) and ('Checkpoint' == linewords[3]):
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))

                table['log_type'].append('checkpoint')
                table['time'].append(epoch_time)
                table['run_ID'].append(run_ID)
                table['parsable_string'].append("class={0};name={1};status={2};RAM={3}".format(
                        linewords[4], linewords[5], linewords[6], linewords[7]))

            # Special gene multiprocessing reporting
            elif 'SpecialPoint_genes' in line:
                pstring = "scaffold={0};PID={1};status={2};what={3}".format(
                            linewords[1], linewords[3], linewords[5], linewords[4])

                table['log_type'].append('Special_genes')
                table['time'].append(float(linewords[6]))
                table['parsable_string'].append(pstring)
                table['run_ID'].append(run_ID)

            elif "Plotting plot" in line:
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))
                pstring = "plot={0}".format(
                            linewords[5])

                table['log_type'].append('Plotting')
                table['time'].append(epoch_time)
                table['parsable_string'].append(pstring)
                table['run_ID'].append(run_ID)

            # Profile reporting
            elif (len(linewords) > 3) and ('WorkerLog' == linewords[0]):
                pstring = "unit={0};PID={1};status={2};process_RAM={3};command={4}".format(
                        linewords[2], linewords[6], linewords[3], linewords[4], linewords[1])

                table['log_type'].append(linewords[0])
                table['time'].append(float(linewords[5]))
                table['parsable_string'].append(pstring)
                table['run_ID'].append(run_ID)

            # Group log
            elif (len(linewords) > 3) and ('GroupLog' == linewords[0]):
                pstring = "unit={0};PID={1};status={2};process_RAM={3};command={4};task={5}".format(
                    linewords[2], linewords[6], linewords[3], linewords[4], linewords[1], linewords[7])

                table['log_type'].append(linewords[0])
                table['time'].append(float(linewords[5]))
                table['parsable_string'].append(pstring)
                table['run_ID'].append(run_ID)

            # Special Failure
            elif 'FAILURE' in line:
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))
                fail_type = linewords[4]

                if fail_type == 'FilterReads':
                    pstring = "type={0};scaffold={1}".format(fail_type, linewords[5])

                elif fail_type == 'SplitException':
                    pstring = "type={0};scaffold={1};split={2}".format(fail_type, linewords[5], linewords[6])

                elif fail_type == 'MergeError':
                    pstring = "type={0};scaffold={1}".format(fail_type, linewords[5])

                elif fail_type == 'GeneException':
                    pstring = "type={0};scaffold={1}".format(fail_type, linewords[5])

                elif fail_type == 'StbError':
                    pstring = "type={0};scaffold={1};bin={2};stb_type={3}".format(
                                fail_type, linewords[5], linewords[6], linewords[7])

                elif fail_type == 'iRepError':
                    pstring = "type={0};genome={1};mm={2}".format(
                                fail_type, linewords[5], linewords[6])

                else:
                    pstring = "type={0}".format(fail_type)

                table['log_type'].append('Failure')
                table['time'].append(float(epoch_time))
                table['parsable_string'].append(pstring)
                table['run_ID'].append(run_ID)

            # Failture that needs to be captured better
            elif 'Double failure!' in line:
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))

                table['log_type'].append('Failure')
                table['time'].append(epoch_time)
                table['parsable_string'].append("error={0}".format(line.strip()))
                table['run_ID'].append(run_ID)

            prev_line_2 = prev_line
            prev_line = line

    Ldb = pd.DataFrame(table)

    # Filter the log
    if most_recent:
        Ldb = filter_most_recent(Ldb)
        assert len(Ldb['run_ID'].unique()) == 1

    return Ldb

def generate_reports(Ldb, MLdb, debug=False):
    """
    Using the parsed log dataframe, make "name2report".

    The following reports are made:
    * Overall
    * Checkpoints
    (lots more; a lot of this is redundent)

    """
    name2report = {}
    assert len(Ldb['run_ID'].unique()) == 1

    # Make the overall report
    name = 'Overall'
    try:
        report, OVERALL_RUNTIME = gen_overall_report(Ldb)
        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    # Make the main checkpoint report
    name = 'Checkpoints'
    try:
        report = gen_checkpoint_report(Ldb, overall_runtime=OVERALL_RUNTIME.total_seconds(), log_class='main_profile')
        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    # Make the filter reads reaport
    name = "Filter reads report"
    try:
        # Checkpoint log
        report = ''
        report += gen_checkpoint_report(Ldb, log_class='FilterReads')

        # Set up multiprocessing log
        report += '\n'
        report += gen_multiprocessing_report(MLdb, commands=['GetPairedReads'])
        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    # Make the profile RAM report
    name = 'Profile report'
    try:
        # Set up checkpoint log
        report = ''
        report += gen_checkpoint_report(Ldb, log_class='Profile')

        # Set up multiprocessing log
        report += '\n* Profiling splits *\n'
        report += gen_multiprocessing_report(MLdb, commands=['SplitProfile'])
        report += '\n* Merging splits and profiling genes *\n'
        report += gen_multiprocessing_report(MLdb, commands=['MergeProfile'])

        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    # # Make the genes report
    # name = 'Genes paralellization efficiency'
    # try:
    #     report = _gen_genes_report(Ldb)
    #     name2report[name] = report
    # except BaseException as e:
    #     if debug:
    #         print('Failed to make log for {0} - {1}'.format(name, str(e)))
    #         traceback.print_exc()

    # Make the genes report
    name = 'Geneome level report'
    try:
        report = ''
        report += gen_checkpoint_report(Ldb, log_class='GenomeLevel')
        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    name = 'Plotting'
    try:
        report = _gen_plotting_report(Ldb)
        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    name = 'Compare'
    try:
        # Checkpoint log
        report = ''
        report += gen_checkpoint_report(Ldb, log_class='Compare')

        # Set up multiprocessing log
        report += '\n'
        report += gen_multiprocessing_report(MLdb, commands=['Compare'])

        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    # Make the failure report
    name = 'Failures'
    try:
        report = _gen_failures_report(Ldb)
        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    return name2report

def gen_overall_report(Ldb):
    """
    Report on the overall runtime of the program
    """
    report = ''
    if 'program_end' in Ldb['log_type'].tolist():
        start = datetime.fromtimestamp(Ldb[Ldb['log_type'] == 'program_start']['time'].tolist()[0])
        end = datetime.fromtimestamp(Ldb[Ldb['log_type'] == 'program_end']['time'].tolist()[0])
        runtime = end - start

        start = start.strftime('%Y-%m-%d %H:%M:%S')
        end = end.strftime('%Y-%m-%d %H:%M:%S')
        n2a = parse_parsable_string(Ldb[Ldb['log_type'] == 'program_start']['parsable_string'].tolist()[0])

        report += 'InStrain version {3} started at {0} and ended at {1}.\nRuntime = {2}\n'.format(
                    start, end, td_format(runtime), n2a['version'])
        report += 'Command = {0}\n'.format(n2a['cmd'])
    else:
        start = datetime.fromtimestamp(Ldb[Ldb['log_type'] == 'program_start']['time'].tolist()[0])
        end = datetime.fromtimestamp(Ldb['time'].max())
        runtime = end - start

        start = start.strftime('%Y-%m-%d %H:%M:%S')
        end = end.strftime('%Y-%m-%d %H:%M:%S')
        n2a = parse_parsable_string(Ldb[Ldb['log_type'] == 'program_start']['parsable_string'].tolist()[0])

        report += 'InStrain version {3} started at {0} and ended at {1}.\nCOMMAND FAILED AND DID NOT FINISH.\nRuntime = {2}\n'.format(
                    start, end, td_format(runtime), n2a['version'])
        report += 'Command = {0}\n'.format(n2a['cmd'])

    OVERALL_RUNTIME = runtime
    return report, OVERALL_RUNTIME

def load_multiprocessing_log(Odb):
    """
    Generate a log table specifically for "WorkerLog" types
    """
    # Parse the initial datatable
    ldb = Odb[(Odb['log_type'] == 'WorkerLog')]
    table = defaultdict(list)
    for i, row in ldb.iterrows():
        for thing, value in parse_parsable_string(row['parsable_string']).items():
            table[thing].append(value)
        table['time'].append(row['time'])
    Ldb = pd.DataFrame(table)
    Ldb['multi_log_type'] = 'WorkerLog'

    ldb = Odb[(Odb['log_type'] == 'GroupLog')]
    table = defaultdict(list)
    for i, row in ldb.iterrows():
        for thing, value in parse_parsable_string(row['parsable_string']).items():
            table[thing].append(value)
        table['time'].append(row['time'])
    Gdb = pd.DataFrame(table)
    if len(Gdb) > 0:
        Gdb['multi_log_type'] = 'GroupLog'

    return pd.concat([Ldb, Gdb]).reset_index(drop=True)

def gen_multiprocessing_report(MLdb, commands):
    """
    From MLdb
    """
    # Load multiprocessing
    if len(MLdb) > 0:
        mldb = MLdb[MLdb['command'].isin(commands)]
        if len(mldb) > 0:
            ldb, gdb = parse_multiprocessing(mldb)
            return gen_multiprocessing_text(ldb, gdb)

    return ''

def parse_multiprocessing(Odb):
    """
    Ldb looks like:

                              unit    PID status  process_RAM  command  \
    0    N5_271_010G1_scaffold_132  32031  end    188485632.0  Compare
    1    N5_271_010G1_scaffold_75   32032  end    188555264.0  Compare
    2    N5_271_010G1_scaffold_34   32033  end    188387328.0  Compare
    3    N5_271_010G1_scaffold_152  32035  end    187920384.0  Compare
    4    N5_271_010G1_scaffold_135  32036  end    188030976.0  Compare
    ..                         ...    ...  ...            ...      ...
    173  N5_271_010G1_scaffold_96   32035  end    198172672.0  Compare
    174  N5_271_010G1_scaffold_15   32033  end    201838592.0  Compare
    175  N5_271_010G1_scaffold_32   32034  end    199385088.0  Compare
    176  N5_271_010G1_scaffold_10   32031  end    198950912.0  Compare
    177  N5_271_010G1_scaffold_8    32036  end    199233536.0  Compare

                 time
    0    1.597446e+09
    1    1.597446e+09
    2    1.597446e+09
    3    1.597446e+09
    4    1.597446e+09


    """
    table = defaultdict(list)
    Ldb = Odb[Odb['multi_log_type'] == 'WorkerLog']
    Ldb.loc[:, 'time'] = Ldb['time'].astype(float)
    Ldb.loc[:, 'process_RAM'] = Ldb['process_RAM'].astype(float)
    first_time = Ldb['time'].min()

    # Generate this on a per-unit level
    for scaffold, ddb in Ldb.groupby('unit'):
        for cmd, db in ddb.groupby('command'):
            sdb = db[db['status'] == 'start']
            edb = db[db['status'] == 'end']

            table['unit'].append(scaffold)
            table['PID'].append(db['PID'].iloc[0])

            table['start_time'].append(sdb['time'].iloc[0])
            table['adjusted_start'].append(sdb['time'].iloc[0] - first_time)
            table['start_process_RAM'].append(sdb['process_RAM'].iloc[0])

            if len(edb) > 0:
                table['adjusted_end'].append(edb['time'].iloc[0] - first_time)
                table['end_process_RAM'].append(edb['process_RAM'].iloc[0])
                table['end_time'].append(edb['time'].iloc[0])
            else:
                for i in ['adjusted_end', 'end_process_RAM', 'end_time']:
                    table[i].append(np.nan)

            table['runs'].append(len(sdb))
            table['command'].append(cmd)

    db = pd.DataFrame(table)
    db.loc[:, 'runtime'] = [s - e for s, e in zip(db['end_time'], db['start_time'])]
    db.loc[:, 'RAM_usage'] = [s - e for s, e in zip(db['end_process_RAM'], db['start_process_RAM'])]
    WorkerDB = db

    table = defaultdict(list)
    Ldb = Odb[Odb['multi_log_type'] == 'GroupLog']
    if len(Ldb) > 0:
        Ldb.loc[:, 'time'] = Ldb['time'].astype(float)
        Ldb.loc[:, 'process_RAM'] = Ldb['process_RAM'].astype(float)
        first_time = Ldb['time'].min()

        # Generate this on a per-unit level
        for scaffold, ddb in Ldb.groupby('unit'):
            for cmd, db in ddb.groupby('command'):
                sdb = db[db['status'] == 'start']
                edb = db[db['status'] == 'end']

                table['unit'].append(scaffold)
                table['PID'].append(db['PID'].iloc[0])

                table['start_time'].append(sdb['time'].iloc[0])
                table['adjusted_start'].append(sdb['time'].iloc[0] - first_time)
                table['start_process_RAM'].append(sdb['process_RAM'].iloc[0])

                if len(edb) > 0:
                    table['adjusted_end'].append(edb['time'].iloc[0] - first_time)
                    table['end_process_RAM'].append(edb['process_RAM'].iloc[0])
                    table['end_time'].append(edb['time'].iloc[0])
                else:
                    for i in ['adjusted_end', 'end_process_RAM', 'end_time']:
                        table[i].append(np.nan)

                table['runs'].append(len(sdb))
                table['command'].append(cmd)

        db = pd.DataFrame(table)
        db.loc[:, 'runtime'] = [s - e for s, e in zip(db['end_time'], db['start_time'])]
        db.loc[:, 'RAM_usage'] = [s - e for s, e in zip(db['end_process_RAM'], db['start_process_RAM'])]
        GroupDB = db
    else:
        GroupDB = pd.DataFrame()

    return WorkerDB, GroupDB

def gen_multiprocessing_text(rdb, gdb, name='unit'):
    report = ''

    # Overall wall time
    start = datetime.fromtimestamp(rdb['start_time'].min())
    end = datetime.fromtimestamp(rdb['end_time'].max())
    runtime = end - start

    # User time
    parallel_time = rdb['runtime'].sum()
    avg_time = rdb['runtime'].mean()

    # Number of processes used
    PIDs = len(rdb['PID'].unique())

    report += "{0:30}\t{1}\n".format("Wall time", td_format(runtime))
    report += "{0:30}\t{1}\n".format("Total processes used", PIDs)
    report += "{0:30}\t{1:.1f}\n".format("Average number processes used", parallel_time/runtime.total_seconds())
    report += "{0:30}\t{1:.1f}%\n".format("Paralellization efficiency", (parallel_time/runtime.total_seconds()/(PIDs))*100)
    report += "{0:30}\t{1}\n".format("Units profiled", len(rdb['unit'].unique()))

    # Report on splits
    #report += "\n"
    report += "{0:30}\t{1}\n".format("Average time per unit", td_format(None, seconds=rdb['runtime'].mean()))
    report += "{0:30}\t{1}\n".format("Median time per unit", td_format(None, seconds=rdb['runtime'].median()))
    report += "{0:30}\t{1}\n".format("Maximum unit time", td_format(None, seconds=rdb['runtime'].max()))
    report += "{0:30}\t{1}\n".format("Longest running unit", rdb.sort_values('runtime', ascending=False)['unit'].iloc[0])
    report += "{0:30}\t{1}\n".format("Per-process efficiency", sorted(["{0:.1f}".format((d['runtime'].sum()/(rdb['end_time'].max() - d['start_time'].min()))*100) for p, d in rdb.groupby('PID')]))

    # Report on RAM
    #report += "\n"
    report += "{0:35}\t{1}\n".format("{0} per-process strating RAM".format(name),
                ["{0}".format(humanbytes(d['start_process_RAM'].iloc[0])) for p, d in rdb.groupby('PID')])
    report += "{0:35}\t{1}\n".format("{0} per-process final RAM".format(name),
                ["{0}".format(humanbytes(d['start_process_RAM'].iloc[-1])) for p, d in rdb.groupby('PID')])
    report += "{0:35}\t{1}\n".format("{0} per-process minimum RAM".format(name),
                ["{0}".format(humanbytes(d['start_process_RAM'].min())) for p, d in rdb.groupby('PID')])
    report += "{0:35}\t{1}\n".format("{0} per-process maximum RAM".format(name),
                ["{0}".format(humanbytes(d['start_process_RAM'].max())) for p, d in rdb.groupby('PID')])

    # Report on groups
    if len(gdb) > 0:
        report += "{0:30}\t{1}\n".format("Number of groups", len(gdb['unit'].unique()))
        report += "{0:30}\t{1}\n".format("Average time per group", td_format(None, seconds=gdb['runtime'].mean()))
        report += "{0:30}\t{1}\n".format("Median time per group", td_format(None, seconds=gdb['runtime'].median()))

    return report

def gen_checkpoint_report(ldb, overall_runtime=None, log_class=None):
    """
    The prefered way of handling checkpoints

    ldb should have the columns:
        name = name of checkpoint
        status = start or end
        RAM = RAM usage
        time = time in epoch time
    """
    # Parse the parseable strings
    ldb = ldb[ldb['log_type'] == 'checkpoint']
    if len(ldb) > 0:
        for i in ['name', 'class', 'status', 'RAM']:
            ldb.loc[:, i] = [parse_parsable_string(pstring)[i] for pstring in ldb['parsable_string']]

    else:
        return ''

    # Subset to a specific class
    if log_class is not None:
        ldb = ldb[ldb['class'] == log_class]
    if len(ldb) == 0:
        return ''

    # Handle the overall runtime
    if overall_runtime is None:
        overall_runtime = max(ldb[ldb['status'] == 'end']['time'].max() - ldb[ldb['status'] == 'start']['time'].min(), 1)

    report = ''
    order = list(ldb['name'].unique())

    for name in order:
        db = ldb[ldb['name'] == name]
        if len(db) > 2:
            report += '{0} has problems and cannot be reported\n'.format(name)

        elif len(db) == 2:
            start = db[db['status'] == 'start'].iloc[0]['time']
            end = db[db['status'] == 'end'].iloc[0]['time']
            runtime = end - start

            if 'RAM' in ldb.columns:
                startram = int(db[db['status'] == 'start']['RAM'].tolist()[0])
                endram = int(db[db['status'] == 'end']['RAM'].tolist()[0])
                ram_change = endram-startram

                if ram_change > 0:
                    inc_dec = 'increased'
                else:
                    inc_dec = 'decreased'

                report += '{0:20} took {1:15} ({2:4.1f}% of overall)\tRAM went from {5} to {6} ({4} by {3})\n'.format(
                            name, td_format(None, seconds=runtime), (runtime/overall_runtime)*100,
                            humanbytes(ram_change, sign=False), inc_dec,
                            humanbytes(startram, sign=False),
                            humanbytes(endram, sign=False))

            else:
                report += '{0:20} took {1:15} ({2:4.1f}% of overall)\n'.format(name, td_format(runtime), (runtime/overall_runtime)*100)

        elif len(db) == 1:
            try:
                start = db[db['status'] == 'start'].iloc[0]['time']
                report += '{0:20} started at {1} and never finished\n'.format(name, start)
            except:
                report += '{0:20} failed\n'.format(name)

    return report

def _gen_genes_report(Ldb, detailed=False):
    report = ''

    # Set up checkpoint log
    report += gen_checkpoint_report(Ldb, log_class='GeneProfile')
    if report != '':
        report += '\n'

    # Set up paralellization log
    ldb = Ldb[Ldb['log_type'] == 'Special_genes']
    PGdb = _load_genes_logtable(ldb)

    if len(PGdb) == 0:
        return report

    # Generate report on paralellization as a whole
    rdb = PGdb[PGdb['command'] == 'whole']

    start = datetime.fromtimestamp(rdb['start_time'].min())
    end = datetime.fromtimestamp(rdb['end_time'].max())
    runtime = end - start
    parallel_time = rdb['runtime'].sum()
    avg_time = rdb['runtime'].mean()
    PIDs = len(rdb['PID'].unique())

    # Report on paralelized steps as a whole
    report += "{0:30}\t{1}\n".format("Wall time parallelized steps", td_format(runtime))
    report += "{0:30}\t{1}\n".format("Total processes used", PIDs)
    report += "{0:30}\t{1:.1f}\n".format("Average number processes used", parallel_time/runtime.total_seconds())
    report += "{0:30}\t{1:.1f}%\n".format("Paralellization efficiency", (parallel_time/runtime.total_seconds()/(PIDs))*100)
    report += "{0:30}\t{1}\n".format("Scaffolds profiled", len(rdb['scaffold'].unique()))
    report += "{0:30}\t{1}\n".format("Average time per scaffold", td_format(None, seconds=rdb['runtime'].mean()))
    report += "{0:30}\t{1}\n".format("Median time per scaffold", td_format(None, seconds=rdb['runtime'].median()))
    report += "{0:30}\t{1}\n".format("Maximum split scaffold", td_format(None, seconds=rdb['runtime'].max()))
    report += "{0:30}\t{1}\n".format("Longest running scaffold", rdb.sort_values('runtime', ascending=False)['scaffold'].iloc[0])
    report += "{0:30}\t{1}\n".format("Per-process efficiency", sorted(["{0:.1f}".format((d['runtime'].sum()/(rdb['end_time'].max() - d['start_time'].min()))*100) for p, d in rdb.groupby('PID')]))

    # Report on sub-paralalized steps
    for step, db in PGdb.groupby('command'):
        if step == 'whole':
            continue

        step_start = datetime.fromtimestamp(db['start_time'].min())
        step_start = datetime.fromtimestamp(db['end_time'].max())
        step_parallel_time = db['runtime'].sum()

        report += "{0:30}\t{1:.1f}% of parallel runtime\n".format("Step {0}: ".format(step), (step_parallel_time/parallel_time)*100)

    return report

def _gen_failures_report(Ldb):
    report = ''
    ldb = Ldb[Ldb['log_type'] == 'Failure']
    if len(ldb) > 0:
        ldb.loc[:, 'failure_type'] = [parse_parsable_string(p)['type'] for p in ldb['parsable_string']]

        for t, db in ldb.groupby('failure_type'):
            table = defaultdict(list)
            for i, row in db.iterrows():
                for thing, value in parse_parsable_string(row['parsable_string']).items():
                    table[thing].append(value)
                table['time'].append(row['time'])
            fdb = pd.DataFrame(table)

            if t == 'FilterReads':
                report += "The following scaffolds were not in the bam file:\n"
                for s in fdb['scaffold'].unique():
                    report += s + '\n'
                report += '\n'

            elif t == 'SplitException':
                report += "The following splits failed during profiling:\n"
                for i, row in fdb.iterrows():
                    report += "{0} split {1}\n".format(row['scaffold'], row['split'])
                report += '\n'

            elif t == 'MergeError':
                report += "The following scaffolds could not be profiled due to mering errors:\n"
                for i, row in fdb.iterrows():
                    report += "{0}\n".format(row['scaffold'])
                report += '\n'

            elif t == 'GeneException':
                report += "Genes on the following scaffolds could not be profiled due to errors durring profiling:\n"
                for i, row in fdb.iterrows():
                    report += "{0}\n".format(row['scaffold'])
                report += '\n'

            elif t == 'StbError':
                report += "The following scaffolds were in the .stb file given, but not the original .fasta file " \
                        + "used for profiling. They will not be considered in genomeLevel operations:\n"
                for i, row in fdb.iterrows():
                    report += "{0} (intended for genome {1})\n".format(row['scaffold'], row['bin'])
                report += '\n'

            elif t == 'iRepError':
                report += "The following genomes failed to calculate iRep for an unknown reason:\n"
                for i, row in fdb.iterrows():
                    report += "{0} (mm {1})\n".format(row['genome'], row['mm'])
                report += '\n'

            else:
                report += "I dont know how to report {0} failures\n".format(t)
                for i, row in fdb.iterrows():
                    report += str(row) + '\n'
                report += '\n'

    if report == '':
        report = "No failures"
    return report

def _gen_plotting_report(Ldb):
    report = ''
    ldb = Ldb[Ldb['log_type'] == 'Plotting']
    if len(ldb) > 0:
        #ldb['plot'] = [parse_parsable_string(p)['plot'] for p in ldb['parsable_string']]
        ldb.loc[:, 'plot'] = [parse_parsable_string(p)['plot'] for p in ldb['parsable_string']]
        ldb = ldb.sort_values('time').reset_index(drop=True)

        for i, (index, row) in enumerate(ldb.iterrows()):
            if i == len(ldb) - 1:
                break

            if row['plot'] == 'finished':
                break

            start = row['time']
            end = ldb.iloc[i+1]['time']
            report += "Plot {0} took {1}\n".format(row['plot'], td_format(None, seconds=end-start))

    return report

def _load_genes_logtable(ldb):
    table = defaultdict(list)
    for i, row in ldb.iterrows():
        for thing, value in parse_parsable_string(row['parsable_string']).items():
            table[thing].append(value)
        table['time'].append(row['time'])
    Ldb = pd.DataFrame(table)

    if len(Ldb) == 0:
        return Ldb

    table = defaultdict(list)
    Ldb.loc[:, 'time'] = Ldb['time'].astype(float)
    first_time = Ldb['time'].min()
    for scaffold, ddb in Ldb.groupby('scaffold'):
        for cmd, db in ddb.groupby('what'):
            sdb = db[db['status'] == 'start']
            edb = db[db['status'] == 'end']

            table['scaffold'].append(scaffold)
            table['PID'].append(db['PID'].tolist()[0])
            table['start_time'].append(sdb['time'].tolist()[0])
            table['end_time'].append(edb['time'].tolist()[0])
            table['adjusted_start'].append(sdb['time'].tolist()[0] - first_time)
            table['adjusted_end'].append(edb['time'].tolist()[0] - first_time)
            table['runs'].append(len(sdb))
            table['command'].append(cmd)

    db = pd.DataFrame(table)
    db.loc[:, 'runtime'] = [s-e for s,e in zip(db['end_time'], db['start_time'])]

    return db

def td_format(td_object, seconds=None):
    if seconds is None:
        seconds = int(td_object.total_seconds())
    periods = [
        ('year',        60*60*24*365),
        ('month',       60*60*24*30),
        ('day',         60*60*24),
        ('hour',        60*60),
        ('minute',      60),
        ('second',      1)
    ]

    strings=[]
    for period_name, period_seconds in periods:
        if seconds > period_seconds:
            period_value , seconds = divmod(seconds, period_seconds)
            has_s = 's' if period_value > 1 else ''
            strings.append("%s %s%s" % (period_value, period_name, has_s))

    if len(strings) > 0:
        return ", ".join(strings)
    else:
        return "<1 second"

def humanbytes(B, sign=True):
    if B < 0:
        if sign:
            unit = '- '
        else:
            unit = ''
        B = B * -1
    else:
        unit = ''

    B = float(B)
    KB = float(1024)
    MB = float(KB ** 2) # 1,048,576
    GB = float(KB ** 3) # 1,073,741,824
    TB = float(KB ** 4) # 1,099,511,627,776

    if B < KB:
        return '{2}{0} {1}'.format(B,'Bytes' if 0 == B > 1 else 'Byte', unit)
    elif KB <= B < MB:
        return '{1}{0:.2f} KB'.format(B/KB, unit)
    elif MB <= B < GB:
        return '{1}{0:.2f} MB'.format(B/MB, unit)
    elif GB <= B < TB:
        return '{1}{0:.2f} GB'.format(B/GB, unit)
    elif TB <= B:
        return '{1}{0:.2f} TB'.format(B/TB, unit)

def parse_parsable_string(pstring):
    object2string = {}
    linewords = pstring.split(';')
    for word in linewords:
        ws = word.split('=')
        object2string[ws[0].strip()] = ws[1].strip()
    return object2string

def filter_most_recent(Ldb):
    """
    Only keep the most recent run
    """
    ID = Ldb.sort_values('time')['run_ID'].tolist()[-1]
    return Ldb[Ldb['run_ID'] == ID]

def log_fmt_to_epoch(ttime):
    # Old log format with no year
    if len(ttime.split('-')) == 2:
        oldformat = '%m-%d %H:%M'
        datetimeobject = datetime.strptime(ttime,oldformat)
        datetimeobject = datetimeobject.replace(year=datetime.now().year)
    # New log format with year
    else:
        oldformat = '%y-%m-%d %H:%M:%S'
        datetimeobject = datetime.strptime(ttime,oldformat)

    return datetimeobject.timestamp()

def log_checkpoint(log_class, name, status, inc_children=True):
    """
    Log a checkpoint for the program in the debug under "log" status

    Arguments:
        log_class    = where is this log comming from ("main_profile", "GeneProfile", etc.)
        name         = what is the task that you're logging ("load globals", "run loop", etc.)
        stauts       = either the string "start" and "end"
        inc_children = include child processes as well in RAM tally

    Results:
        Make the following log message:
        "Checkpoint class task status RAM"

        When split on tabs, this has the following structure:
        "Checkpoint" = linewords[3]
        class = linewords[4]
        task = linewords[5]
        start/end = linewords[6]
        RAM = linewords[7]
    """
    current_process = psutil.Process(os.getpid())
    mem = current_process.memory_info().rss
    if inc_children:
        for child in current_process.children(recursive=True):
            try:
                mem += child.memory_info().rss
            except:
                pass

    assert status in ['start', 'end'], [log_class, name, status]
    assert len(name.split()) == 1

    logging.debug("Checkpoint {0} {1} {2} {3}".format(log_class, name, status, mem))

def get_worker_log(worker_type, unit, status, inc_children=False):
    """
    Return a string with log information intended to be generated within a worker process to track how long
    actual multiprocessing is taking place

    Arguments:
        worker_type  = The type of worker this is
        unit         = The unit this worker is currently processing
        status       = start / end
        inc_children = include child processes as well in RAM tally

    Returns:
        A string with the following structure:
        "WorkerLog worker_type unit status RAM time PID"

        "WorkerLog" = linewords[0]
        worker_type = linewords[1]
        unit = linewords[2]
        status = linewords[3]
        ram = linewords[4]
        time = linewords[5]
        PID = linewords[6]

    """
    pid = os.getpid()
    current_process = psutil.Process(pid)
    mem = current_process.memory_info().rss
    if inc_children:
        for child in current_process.children(recursive=True):
            try:
                mem += child.memory_info().rss
            except:
                pass

    assert status in ['start', 'end'], status

    return "\nWorkerLog {0} {1} {2} {3} {4} {5}".format(worker_type, unit, status, mem, time.time(), pid)

def get_group_log(worker_type, unit, status, task='overall', inc_children=False):
    """
    Return a string with log information intended to track how long a group of units takes

    Arguments:
        worker_type  = The type of worker this is (Profile, Merge, etc.)
        unit         = The unit this worker is currently processing
        status       = start / end
        task         = what are you doing?
        inc_children = include child processes as well in RAM tally

    Returns:
        A string with the following structure:
        "WorkerLog worker_type unit status RAM time PID"

        "GroupLog" = linewords[0]
        worker_type = linewords[1]
        unit = linewords[2]
        status = linewords[3]
        ram = linewords[4]
        time = linewords[5]
        PID = linewords[6]

    """
    pid = os.getpid()
    current_process = psutil.Process(pid)
    mem = current_process.memory_info().rss
    if inc_children:
        for child in current_process.children(recursive=True):
            try:
                mem += child.memory_info().rss
            except:
                pass

    assert status in ['start', 'end'], status

    return "\nGroupLog {0} {1} {2} {3} {4} {5} {6}".format(worker_type, unit, status, mem, time.time(), pid, task)