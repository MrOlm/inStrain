import os
import sys
import h5py
import time
import shutil
import logging
import traceback
import numpy as np
import pandas as pd
from datetime import datetime
from datetime import timedelta
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.ticker as ticker
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

import seaborn as sns

def process_logs(IS_loc):
    logging.shutdown()
    logloc = os.path.join(IS_loc, 'log/log.log')
    report_run_stats(logloc, most_recent=False, printToo=True, save=True)

def report_run_stats(logloc, save=True, most_recent=True, printToo=True, debug=False, plot=True):
    if logloc == None:
        return

    # Load the log
    try:
        Ldb = load_log(logloc)

        # Filter the log
        if most_recent:
            Ldb = filter_most_recent(Ldb)
    except BaseException as e:
        if debug:
            print('Failed to load log file - {1}'.format('None', str(e)))
            traceback.print_exc()
        return

    # Generate reports
    for run, ldb in Ldb.groupby('run_ID'):
        name2report = generate_reports(ldb, debug=debug)

        # Print
        if printToo:
            for name, report in name2report.items():
                print("..:: {0} ::..\n{1}".format(name, report))

        if save == True:
            saveloc = logloc.replace('log.log', '{0}.runtime_summary.txt'.format(run))
            figloc = logloc.replace('log.log', '{0}.ramProfile.png'.format(run))
        elif save == False:
            continue
        else:
            saveloc = save
            figloc = save + '.ramProfile.png'

        with open(saveloc, 'w') as o:
            for name, report in name2report.items():
                o.write("..:: {0} ::..\n{1}\n".format(name, report))

    # Make the plot
    if plot:
        try:
            profile_plot(Ldb, saveloc=figloc)
        except BaseException as e:
            if debug:
                print('Failed to make profile plot - {1}'.format('None', str(e)))
                traceback.print_exc()

def profile_plot(Ldb, saveloc=None):
    ldb = Ldb[Ldb['log_type'] == 'Special_profile']
    rdb, sys_ram = _load_profile_logtable(ldb)

    plt.scatter(rdb['adjusted_start'], rdb['percent_RAM'])
    plt.xlabel('Runtime (seconds)')
    plt.ylabel('System RAM available at start of thread\n(% of the {0} system)'.format(humanbytes(sys_ram)))

    if saveloc != None:
        plt.gcf().savefig(saveloc, bbox_inches='tight')

def generate_reports(Ldb, debug=False):
    name2report = {}
    assert len(Ldb['run_ID'].unique()) == 1

    # Make the overall report
    name = 'Overall'
    try:
        report, OVERALL_RUNTIME = _gen_overall_report(Ldb)
        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    # Make the checkpoint report
    name = 'Checkpoints'
    try:
        report = _gen_checkpoint_report(Ldb, OVERALL_RUNTIME)
        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    # Make the profile RAM report
    name = 'Profile RAM useage and paralellization efficiency'
    try:
        report = _gen_profileRAM_report(Ldb)
        name2report[name] = report
    except BaseException as e:
        if debug:
            print('Failed to make log for {0} - {1}'.format(name, str(e)))
            traceback.print_exc()

    # Make the genes report
    name = 'Genes paralellization efficiency'
    try:
        report = _gen_genes_report(Ldb)
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

def _gen_overall_report(Ldb):
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

def _gen_checkpoint_report(Ldb, overall_runtime):
    report = ''

    # Set up
    ldb = Ldb[Ldb['log_type'] == 'checkpoint']
    for i in ['name', 'status']:
        ldb[i] = [parse_parsable_string(pstring)[i] for pstring in ldb['parsable_string']]
    try:
        i = 'RAM'
        ldb[i] = [parse_parsable_string(pstring)[i] for pstring in ldb['parsable_string']]
    except:
        pass

    # actually go
    order = list(ldb['name'].unique())
    for name in order:
        db = ldb[ldb['name'] == name]
        if len(db) > 2:
            report += '{0} has problems and cannot be reported\n'.format(name)
        elif len(db) == 2:
            start = datetime.fromtimestamp(db[db['status'] == 'start']['time'].tolist()[0])
            end = datetime.fromtimestamp(db[db['status'] == 'end']['time'].tolist()[0])
            runtime = end - start

            start = start.strftime('%Y-%m-%d %H:%M:%S')
            end = end.strftime('%Y-%m-%d %H:%M:%S')

            if 'RAM' in ldb.columns:
                startram = int(db[db['status'] == 'start']['RAM'].tolist()[0])
                endram = int(db[db['status'] == 'end']['RAM'].tolist()[0])
                report += '{0:20} took {1:15} ({2:3.1f}% of overall)\tRAM use increased by {3}\n'.format(name, td_format(runtime), (runtime/overall_runtime)*100, humanbytes(endram-startram))
            else:
                report += '{0:20} took {1:15} ({2:3.1f}% of overall)\n'.format(name, td_format(runtime), (runtime/overall_runtime)*100)
        elif len(db) == 1:
            start = datetime.fromtimestamp(db[db['status'] == 'start']['time'].tolist()[0])
            start = start.strftime('%Y-%m-%d %H:%M:%S')
            report += '{0:20} started at {1} and never finished\n'.format(name, start)

    return report

def _gen_profileRAM_report(Ldb, detailed=False):
    '''
    Percent_RAM goes down over the run; it reports the percentage of RAM available
    end_system_RAM describes the total about of ram _available_
    '''
    report = ''

    # Set up
    ldb = Ldb[Ldb['log_type'] == 'Special_profile']
    rdb, sys_ram = _load_profile_logtable(ldb)

    if len(rdb) == 0:
        return ''

    # Report on paralellization
    start = datetime.fromtimestamp(rdb['start_time'].min())
    end = datetime.fromtimestamp(rdb['end_time'].max())
    runtime = end - start

    parallel_time = rdb['runtime'].sum()
    avg_time = rdb['runtime'].mean()

    PIDs = len(rdb['PID'].unique())

    report += "{0:30}\t{1}\n".format("Wall time for Profile", td_format(runtime))
    report += "{0:30}\t{1}\n".format("Total number processes used", PIDs)
    report += "{0:30}\t{1:.1f}\n".format("Average number processes used", parallel_time/runtime.total_seconds())
    report += "{0:30}\t{1:.1f}%\n".format("Paralellization efficiency", (parallel_time/runtime.total_seconds()/PIDs)*100)
    report += "{0:30}\t{1}\n".format("Scaffolds profiled", len(rdb['scaffold'].unique()))
    report += "{0:30}\t{1}\n".format("Average time per scaffold", td_format(None, seconds=rdb['runtime'].mean()))
    report += "{0:30}\t{1}\n".format("Median time per scaffold", td_format(None, seconds=rdb['runtime'].median()))
    report += "{0:30}\t{1}\n".format("Maximum scaffold time", td_format(None, seconds=rdb['runtime'].max()))
    report += "{0:30}\t{1}\n".format("Longest running scaffold", rdb.sort_values('runtime', ascending=False)['scaffold'].iloc[0])
    report += "{0:30}\t{1}\n".format("System RAM available", humanbytes(sys_ram))
    report += "{0:30}\t{1:.1f}%\n".format("Starting RAM usage (%)", 100 - rdb['percent_RAM'].iloc[0]) # Percent ram is the amout AVAILABLE
    report += "{0:30}\t{1:.1f}%\n".format("Ending RAM usage (%)", 100 - rdb['percent_RAM'].iloc[-1])
    report += "{0:30}\t{1}\n".format("Peak RAM used", humanbytes(sys_ram - rdb['end_system_RAM'].min()))
    report += "{0:30}\t{1}\n".format("Mimimum RAM used", humanbytes(sys_ram - rdb['start_system_RAM'].max()))

    report += '{0} scaffolds needed to be run a second time\n'.format(
            len(rdb[rdb['runs'] > 1]['scaffold'].unique()))

    return report

def _gen_genes_report(Ldb, detailed=False):
    report = ''

    # Set up
    ldb = Ldb[Ldb['log_type'] == 'Special_genes']
    rdb = _load_genes_logtable(ldb)

    if len(rdb) == 0:
        return ''

    # Report on paralellization
    start = datetime.fromtimestamp(rdb['start_time'].min())
    end = datetime.fromtimestamp(rdb['end_time'].max())
    runtime = end - start

    parallel_time = rdb['runtime'].sum()
    avg_time = rdb['runtime'].mean()

    PIDs = len(rdb['PID'].unique())

    report += 'Gene calling paralellization occured for {0} and involved {1} processes\n'.format(td_format(runtime), PIDs)
    report += 'An average of {0:.1f} processes were used during this time to profile {1} scaffolds\n'.format(
                parallel_time/runtime.total_seconds(), len(rdb['scaffold'].unique()))
    report += 'Scaffolds ran for an average of {0} (median {1}; longest {2})\n'.format(
                td_format(None, seconds=rdb['runtime'].mean()), td_format(None, seconds=rdb['runtime'].median()),
                td_format(None, seconds=rdb['runtime'].max()))
    report += '{0} scaffolds needed to be run a second time\n'.format(
            len(rdb[rdb['runs'] > 1]['scaffold'].unique()))

    return report

def _gen_failures_report(Ldb):
    report = ''
    ldb = Ldb[Ldb['log_type'] == 'Failure']
    for i, row in ldb.iterrows():
        report += "Failure {0}\n".format(ldb['parsable_string'])
    if report == '':
        report = "No failures"
    return report

def _load_profile_logtable(ldb):
    table = defaultdict(list)
    for i, row in ldb.iterrows():
        for thing, value in parse_parsable_string(row['parsable_string']).items():
            table[thing].append(value)
        table['time'].append(row['time'])
    Ldb = pd.DataFrame(table)

    if len(Ldb) == 0:
        return Ldb, None

    table = defaultdict(list)
    Ldb['time'] = Ldb['time'].astype(float)
    Ldb['process_RAM'] = Ldb['process_RAM'].astype(float)
    Ldb['system_RAM'] = Ldb['system_RAM'].astype(float)
    sys_ram = float(Ldb['total_RAM'].tolist()[0])
    first_time = Ldb['time'].min()
    for scaffold, db in Ldb.groupby('scaffold'):
        sdb = db[db['status'] == 'start']
        edb = db[db['status'] == 'end']

        table['scaffold'].append(scaffold)
        table['PID'].append(db['PID'].tolist()[0])

        table['start_time'].append(sdb['time'].tolist()[0])
        table['adjusted_start'].append(sdb['time'].tolist()[0] - first_time)
        table['start_process_RAM'].append(sdb['process_RAM'].tolist()[0])
        table['start_system_RAM'].append(sdb['system_RAM'].tolist()[0])

        if len(edb) > 0:
            table['adjusted_end'].append(edb['time'].tolist()[0] - first_time)
            table['end_process_RAM'].append(edb['process_RAM'].tolist()[0])
            table['end_system_RAM'].append(edb['system_RAM'].tolist()[0])
            table['end_time'].append(edb['time'].tolist()[0])
        else:
            for i in ['adjusted_end', 'end_process_RAM', 'end_system_RAM', 'end_time']:
                table[i].append(np.nan)

        table['runs'].append(len(sdb))

    db = pd.DataFrame(table)
    db['runtime'] = [s-e for s,e in zip(db['end_time'], db['start_time'])]
    db['RAM_usage'] = [s-e for s,e in zip(db['end_process_RAM'], db['start_process_RAM'])]
    db['percent_RAM'] = [(s/sys_ram) * 100 for s in db['end_system_RAM']]

    return db, sys_ram

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
    Ldb['time'] = Ldb['time'].astype(float)
    first_time = Ldb['time'].min()
    for scaffold, db in Ldb.groupby('scaffold'):
        sdb = db[db['status'] == 'start']
        edb = db[db['status'] == 'end']

        table['scaffold'].append(scaffold)
        table['PID'].append(db['PID'].tolist()[0])
        table['start_time'].append(sdb['time'].tolist()[0])
        table['end_time'].append(edb['time'].tolist()[0])
        table['adjusted_start'].append(sdb['time'].tolist()[0] - first_time)
        table['adjusted_end'].append(edb['time'].tolist()[0] - first_time)
        table['runs'].append(len(sdb))

    db = pd.DataFrame(table)
    db['runtime'] = [s-e for s,e in zip(db['end_time'], db['start_time'])]

    return db

def td_format(td_object, seconds=False):
    if seconds == False:
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

def humanbytes(B):
    if B < 0:
        unit = '- '
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
    '''
    Only keep the most recent run
    '''
    ID = Ldb.sort_values('time')['run_ID'].tolist()[-1]
    return Ldb[Ldb['run_ID'] == ID]

def load_log(logfile):
    table = defaultdict(list)
    with open(logfile) as o:
        prev_line_2 = None
        prev_line = None
        run_ID = None
        for line in o.readlines():
            line = line.strip()

            # load new inStrain run
            if 'inStrain version' in line:
                linewords = [x.strip() for x in line.split()]
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))

                run_ID = datetime.fromtimestamp(epoch_time).strftime('%Y%m%d_%H%M%S')
                cmd = prev_line_2.strip().split('was:')[1].strip()

                table['log_type'].append('program_start')
                table['time'].append(epoch_time)
                table['parsable_string'].append("version={0}; cmd={1}".format(linewords[5], cmd))
                table['run_ID'].append(run_ID)

            # load inStrain run finish
            elif 'inStrain complete' in line:
                linewords = [x.strip() for x in line.split()]
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))

                table['log_type'].append('program_end')
                table['time'].append(epoch_time)
                table['parsable_string'].append("loglog={0}".format(linewords[14][:-1]))
                table['run_ID'].append(run_ID)

            # regular checkpoints
            elif 'Checkpoint' in line:
                linewords = [x.strip() for x in line.split()]
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))

                table['log_type'].append('checkpoint')
                table['time'].append(epoch_time)
                table['run_ID'].append(run_ID)

                if len(linewords) == 9:
                    table['parsable_string'].append("status={0};name={1};RAM={2}".format(linewords[5], linewords[4], linewords[8]))
                else:
                    table['parsable_string'].append("status={0};name={1}".format(linewords[5], linewords[4]))

            # Special gene multiprocessing reporting
            elif 'SpecialPoint_genes' in line:
                linewords = [x.strip() for x in line.split()]
                pstring = "scaffold={0};PID={1};status={2}".format(
                            linewords[1], linewords[3], linewords[4])

                table['log_type'].append('Special_genes')
                table['time'].append(float(linewords[5]))
                table['parsable_string'].append(pstring)
                table['run_ID'].append(run_ID)

            # Special profile RAM and multiprocessing reporting
            elif 'RAM. System has' in line:
                linewords = [x.strip() for x in line.split()]
                pstring = "scaffold={0};PID={1};status={2};process_RAM={3};system_RAM={4};total_RAM={5}".format(
                            linewords[0], linewords[2], linewords[3], linewords[7], linewords[11], linewords[13])

                table['log_type'].append('Special_profile')
                table['time'].append(float(linewords[5]))
                table['parsable_string'].append(pstring)
                table['run_ID'].append(run_ID)
                # table['scaffold'].append(linewords[0])
                # table['PID'].append(linewords[2])
                # table['status'].append(linewords[3])
                # table['time'].append(linewords[5])
                # table['process_RAM'].append(linewords[7])
                # table['system_RAM'].append(linewords[11])
                # table['total_RAM'].append(linewords[13])

            # Special Failure
            elif 'FAILURE' in line:
                linewords = [x.strip() for x in line.split()]
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))
                fail_type = linewords[3]
                if fail_type == 'GENES_OUTER':
                    pstring = "type={0}".format(fail_type)
                else:
                    pstring = "type={0}".format(fail_type)

                table['log_type'].append('Failure')
                table['time'].append(float(epoch_time))
                table['parsable_string'].append(pstring)
                table['run_ID'].append(run_ID)

            # Failture that needs to be captured better
            elif 'Double failure!' in line:
                linewords = [x.strip() for x in line.split()]
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))

                table['log_type'].append('Failure')
                table['time'].append(epoch_time)
                table['parsable_string'].append("error={0}".format(line.strip()))
                table['run_ID'].append(run_ID)

            prev_line_2 = prev_line
            prev_line = line

    Ldb = pd.DataFrame(table)
    #Ldb = add_run_IDs(Ldb)

    return Ldb

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
