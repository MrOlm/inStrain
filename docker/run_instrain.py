#!/usr/bin/env python

import os
import sys
import glob
import shlex
import shutil
import logging
import subprocess
from argparse import ArgumentParser

from s3_utils import download_file, upload_file, download_folder, upload_folder
from job_utils import generate_working_dir, delete_working_dir#, setup_logger


def download_reference(s3_path, working_dir):
    """
    Downloads reference folder that has been configured to run with Isaac
    :param s3_path: S3 path that the folder resides in
    :param working_dir: working directory
    :return: local path to the folder containing the reference
    """

    reference_folder = os.path.join(working_dir, 'reference')

    try:
        os.mkdir(reference_folder)
    except Exception as e:
        pass

    download_folder(s3_path, reference_folder)

    # Update sorted reference
    update_sorted_reference(reference_folder)

    return reference_folder


def download_fastq_files(fastq1_s3_path, fastq2_s3_path, working_dir):
    """
    Downlodas the fastq files
    :param fastq1_s3_path: S3 path containing FASTQ with read1
    :param fastq2_s3_path: S3 path containing FASTQ with read2
    :param working_dir: working directory
    :return: local path to the folder containing the fastq
    """
    fastq_folder = os.path.join(working_dir, 'fastq')

    try:
        os.mkdir(fastq_folder)
    except Exception as e:
        pass

    local_fastq1_path = download_file(fastq1_s3_path, fastq_folder)
    local_fastq2_path = download_file(fastq2_s3_path, fastq_folder)

    # Isaac requires the fastqs to be symlinked as lane1_read1.fastq.gz and lane1_read2.fastq.gz
    os.symlink(local_fastq1_path, os.path.join(fastq_folder, 'lane1_read1.fastq.gz'))
    os.symlink(local_fastq2_path, os.path.join(fastq_folder, 'lane1_read2.fastq.gz'))

    return fastq_folder


def upload_bam(bam_s3_path, local_folder_path):
    """
    Uploads results folder containing the bam file (and associated output)
    :param bam_s3_path: S3 path to upload the alignment results to
    :param local_folder_path: local path containing the alignment results
    """

    upload_folder(bam_s3_path, local_folder_path)


def run_isaac(reference_dir, fastq_folder_path, memory, cmd_args, working_dir):
    """
    Runs Isaac3
    :param reference_dir: local path to directory containing reference
    :param fastq_folder_path: local path to directory containing fastq files
    :param memory: memory in GB
    :param cmd_args: additional command-line arguments to pass in
    :param working_dir: working directory
    :return: path to results
    """

    # Maps to Isaac's folder structure and change working directory
    os.chdir(working_dir)
    bam_folder = os.path.join(working_dir, 'Projects/default/default')

    try:
        os.mkdir(bam_folder)
    except Exception as e:
        pass

    cmd = 'isaac-align -r %s -b %s --base-calls-format fastq-gz -o %s -m %s %s' % \
          (os.path.join(reference_dir, 'sorted-reference.xml'), fastq_folder_path, working_dir, memory, cmd_args)
    print ("Running: %s" % cmd)
    subprocess.check_call(shlex.split(cmd))

    return bam_folder


def update_sorted_reference(reference_dir):
    """
    Updates sorted-reference.xml to map to the correct directory path. Since each analysis occurs in subfolder of the
    working directory, it will change each execution
    :param reference_dir: Reference directory
    """
    with open(os.path.join(reference_dir, 'sorted-reference.xml'), 'r') as infile:
        sorted_reference = infile.read()

    with open(os.path.join(reference_dir, 'sorted-reference.xml'), 'w') as outfile:
        outfile.write(sorted_reference.replace('/scratch', reference_dir))

def download_data(args, working_dir, tmp_dir):
    cmd_string = ''

    # Download the .bam
    if args.bam != None:
        logging.info("Downloading bam to {0}".format(tmp_dir))
        download_file(args.bam, tmp_dir)
        cmd_string = cmd_string + os.path.join(tmp_dir, os.path.basename(args.bam)) + ' '

        try:
            logging.info("Downloading bam index".format(tmp_dir))
            download_file(args.bam + '.bai', tmp_dir)
        except:
            pass

    # Downlaod the .fasta
    if args.fasta != None:
        logging.info("Downloading fasta to {0}".format(tmp_dir))
        download_file(args.fasta, tmp_dir)
        cmd_string = cmd_string + os.path.join(tmp_dir, os.path.basename(args.fasta)) + ' '

    # Download other files
    for f, name in zip([args.genes, args.stb], ['-g', '-s']):
        if f is not None:
            logging.info("{0} is {1}; downloading".format(name, f))
            download_file(f, tmp_dir)
            cmd_string = cmd_string + name + ' ' + os.path.join(tmp_dir, os.path.basename(f)) + ' '

    # 2) Unzip if need be
    to_unzip = glob.glob(tmp_dir + '/*.gz')
    if len(to_unzip) > 0:
        for g in to_unzip:
            cmd = 'gzip -d {0}'.format(g)
            subprocess.check_call(shlex.split(cmd))

    # Get the work directory
    wd_loc = os.path.join(working_dir, args.wd_name)
    cmd_string = cmd_string + ' -o ' + wd_loc

    return cmd_string, wd_loc

def download_checkm(arags):
    pass

def run_inStrain(data_prefix, command, cmd_args, timeout=0):
    '''
    Actually call inStrain

    Args:
        genomes_string: string with genomes, INCLUDING -g argument
        wd_loc: location of output work directory
        command: the string to put right after "command"
        cmd_args: things to put after the command
    '''
    if timeout == 0:
        cmd = 'inStrain {0} {1} {2}'.format(command, data_prefix, cmd_args)
    else:
        cmd = 'timeout {3} inStrain {0} {1} {2}'.format(command, data_prefix, cmd_args, timeout)
    logging.info ("Running: %s" % cmd)
    subprocess.check_call(shlex.split(cmd))

def upload_results(args, working_dir, wd_loc):
    if args.light_upload:
        # delete the raw_data folder
        data_folder = str(glob.glob(wd_loc + '/raw_data/')[0])
        logging.info("deleting {0}".format(data_folder))
        shutil.rmtree(data_folder)

    upload_folder(args.results_directory, working_dir)

def parse_arguments():
    argparser = ArgumentParser()

    file_path_group = argparser.add_argument_group(title='File paths')
    file_path_group.add_argument('--bam', type=str, help='s3 path to the .bam file')
    file_path_group.add_argument('--fasta', type=str, help='s3 path to the .fasta file')
    file_path_group.add_argument('--genes', type=str, help='s3 path to the genes file')
    file_path_group.add_argument('--stb', help='s3 path to the stb file', nargs='*', default=[]))
    file_path_group.add_argument('--results_directory', type=str, help='s3 path to the folder to put the results in')

    run_group = argparser.add_argument_group(title='Run command args')
    run_group.add_argument('--wd_name', type=str, help='Name of the output directory', default='instrain_output')
    run_group.add_argument('--timeout', type=int, help='Kill job after this many seconds and upload what was gotten so far. Input 0 here to have no limit', default=0)
    run_group.add_argument('--command', type=str, help='The command that should go after the command inStrain', default='profile')
    run_group.add_argument('--cmd_args', type=str, help='A string (as long as you want) that will be put after the inStrain command, .bam, fasta, and output directory', default=' ')
    run_group.add_argument('--light_upload', help='By default it will upload the /raw_data folder; this will make it not', default=False, action='store_true')

    argparser.add_argument('--working_dir', type=str, default='/mnt/scratch')

    args = argparser.parse_args()
    return args


def main():
    args = parse_arguments()

    # Make a working directory to be uploaded at the end
    working_dir = generate_working_dir(args.working_dir)
    log_loc = os.path.join(working_dir, 'docker_log.log')
    setup_logger(log_loc)
    logging.info('Here are the args: {0}'.format(args))

    # Get a temporary directory for scratch things
    temp_dir = generate_working_dir('/mnt/temp')

    # Download the required things
    logging.info('downloading data')
    data_prefix, wd_loc = download_data(args, working_dir, temp_dir)

    logging.info('running inStrain command')
    try:
        run_inStrain(data_prefix, args.command, args.cmd_args, args.timeout)
    except:
        pass

    logging.info('uploading results')
    upload_results(args, working_dir, wd_loc)

    logging.info('cleaning up')
    delete_working_dir(working_dir)
    delete_working_dir(temp_dir)

### TO BE MOVED TO COMMON UTILS
def setup_logger(loc):
    ''' set up logger such that DEBUG goes only to file, rest go to file and console '''

    # Cancel if a logger already exists:
    if logging.getLogger('').handlers:
        return

    # set up logging everything to file
    logging.basicConfig(level=logging.DEBUG,
                       format='%(asctime)s %(levelname)-8s %(message)s',
                       datefmt='%m-%d %H:%M',
                       filename=loc)

    # set up logging of INFO or higher to sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    console.setFormatter(formatter)

    logging.getLogger('').addHandler(console)

    logging.debug("!"*80)
    logging.debug("***Logger started up at {0}***".format(loc))
    logging.debug("Command to run was: {0}\n".format(' '.join(sys.argv)))
    logging.debug("!"*80 + '\n')
###

if __name__ == '__main__':
    main()
