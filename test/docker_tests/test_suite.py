#!/usr/bin/env python
'''
Run tests
'''

import io
import os
import sys
import glob
import boto3
import shlex
import docker
import shutil
import pickle
import logging
import warnings

import importlib
import subprocess
import numpy as np
warnings.filterwarnings("ignore")
import pandas as pd
from subprocess import call
from Bio import SeqIO
from collections import defaultdict
from inStrain._version import __version__

def estimate_cost(downfolder, upfolder):
    '''
    downfolder and upfolder can be a string, list, or none
    '''
    n2c = {}
    for var, name in zip([downfolder, upfolder], ['down', 'up']):
        if type(var) == type('test'):
            cost = gb_to_cost(s3_folder_size(var))
        elif type(var) == type([]):
            cost = gb_to_cost(sum([s3_folder_size(x) for x in var]))
        elif downfolder == None:
            cost = 0
        n2c[name] = cost

    downcost = n2c['down']
    upcost = n2c['up']
    print("$$$$$$$$$$$$$$\nIt cost {0:.2f}¢ to run this test ({1}¢ for download and {2}¢ for upload)\n$$$$$$$$$$$$$$".format(
        downcost + upcost, downcost, upcost))

def s3_folder_size(folder):
    '''
    return total size of objects in s3 folder in Gb
    '''
    bucket = folder.split('/')[2]
    key = '/'.join(folder.split('/')[3:])

    total_size = bytesto(sum([float(o['Size']) for o in get_matching_s3_objects(bucket, key)]), 'g')
    return total_size

def bytesto(bytes, to, bsize=1024):
    """convert bytes to megabytes, etc.
       sample code:
           print('mb= ' + str(bytesto(314575262000000, 'm')))
       sample output:
           mb= 300002347.946
    """
    a = {'k' : 1, 'm': 2, 'g' : 3, 't' : 4, 'p' : 5, 'e' : 6 }
    r = float(bytes)
    for i in range(a[to]):
        r = r / bsize
    return(r)

def gb_to_cost(gb):
    '''
    return cost in cents
    '''
    return gb * 0.09

def sync_test_data():
    cmd = 'aws s3 sync ./s3_test_data/ s3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/'
    subprocess.check_call(shlex.split(cmd))
    #print(s3_folder_size('s3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/'))

def clear_s3_results():
    s3 = boto3.resource('s3')
    bucket = s3.Bucket('czbiohub-microbiome')
    bucket.objects.filter(Prefix="Sonnenburg_Lab/Software/docker_testing/s3_results/").delete()

def load_s3_results():
    return [f for f in get_matching_s3_keys('czbiohub-microbiome', 'Sonnenburg_Lab/Software/docker_testing/s3_results/')]

def get_s3_results_folder():
    return 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/s3_results/'

def load_random_test_dir():
    loc = os.path.join(str(os.getcwd()), \
        'test_backend/testdir/')
    return loc

def get_credentials():
    loc = os.path.expanduser('~/.aws/credentials')
    return loc

def get_accessible_test_data():
    loc = os.path.join(str(os.getcwd()), \
        'accessible_test_data/')
    return loc

def read_s3_file(key, bucketname='czbiohub-microbiome'):
    s3 = boto3.resource('s3')
    obj = s3.Object(bucketname, key)
    body = obj.get()['Body'].read()
    return body.decode("utf-8").strip()

def get_matching_s3_objects(bucket, prefix="", suffix=""):
    """
    Generate objects in an S3 bucket.

    :param bucket: Name of the S3 bucket.
    :param prefix: Only fetch objects whose key starts with
        this prefix (optional).
    :param suffix: Only fetch objects whose keys end with
        this suffix (optional).
    """
    s3 = boto3.client("s3")
    paginator = s3.get_paginator("list_objects_v2")

    kwargs = {'Bucket': bucket}

    # We can pass the prefix directly to the S3 API.  If the user has passed
    # a tuple or list of prefixes, we go through them one by one.
    if isinstance(prefix, str):
        prefixes = (prefix, )
    else:
        prefixes = prefix

    for key_prefix in prefixes:
        kwargs["Prefix"] = key_prefix

        for page in paginator.paginate(**kwargs):
            try:
                contents = page["Contents"]
            except KeyError:
                return

            for obj in contents:
                key = obj["Key"]
                if key.endswith(suffix):
                    yield obj


def get_matching_s3_keys(bucket, prefix="", suffix=""):
    """
    Generate the keys in an S3 bucket.

    :param bucket: Name of the S3 bucket.
    :param prefix: Only fetch keys that start with this prefix (optional).
    :param suffix: Only fetch keys that end with this suffix (optional).
    """
    for obj in get_matching_s3_objects(bucket, prefix, suffix):
        yield obj["Key"]

def check_s3_file(floc):
    '''
    Return True if exists and False if it does not
    '''
    bucket = floc.split('/')[2]
    prefix = '/'.join(floc.split('/')[3:])

    found = False
    for key in get_matching_s3_keys(bucket, prefix):
        if prefix in key:
            found = True
    return found

def run_docker(image, cmd, simulate_aegea=True, overwrite_accessible=False):
    '''
    Load the image and run the command
    '''
    # Run the Docker
    cred_loc = get_credentials()
    test_loc = get_accessible_test_data()
    output_loc = load_random_test_dir()
    mt_loc = os.getcwd() + '/mnt'

    if overwrite_accessible:
        cmd = 'cp /root/accessible_testing_data/* ./;' + cmd

    if simulate_aegea == True:
        AEGEA_JUNK = ['/bin/bash', '-c', 'for i in "$@"; do eval "$i"; done; cd /', 'aegea.util.aws.batch', 'set -a', 'if [ -f /etc/environment ]; then source /etc/environment; fi', 'if [ -f /etc/default/locale ]; then source /etc/default/locale; else export LC_ALL=C.UTF-8 LANG=C.UTF-8; fi', 'export AWS_DEFAULT_REGION=us-west-2', 'set +a', 'if [ -f /etc/profile ]; then source /etc/profile; fi', 'set -euo pipefail', 'sed -i -e "s|/archive.ubuntu.com|/us-west-2.ec2.archive.ubuntu.com|g" /etc/apt/sources.list', 'apt-get update -qq', 'apt-get install -qqy --no-install-suggests --no-install-recommends httpie awscli jq lsof python3-virtualenv > /dev/null', 'python3 -m virtualenv -q --python=python3 /opt/aegea-venv', '/opt/aegea-venv/bin/pip install -q argcomplete requests boto3 tweak pyyaml', '/opt/aegea-venv/bin/pip install -q --no-deps aegea==3.4.3']
        cmd = AEGEA_JUNK + [cmd]

    elif simulate_aegea == 'semi':
        AEGEA_JUNK = ['/bin/bash', '-c', 'for i in "$@"; do eval "$i"; done; cd /', 'aegea.util.aws.batch', 'set -a', 'if [ -f /etc/environment ]; then source /etc/environment; fi', 'if [ -f /etc/default/locale ]; then source /etc/default/locale; else export LC_ALL=C.UTF-8 LANG=C.UTF-8; fi', 'export AWS_DEFAULT_REGION=us-west-2', 'set +a', 'if [ -f /etc/profile ]; then source /etc/profile; fi', 'set -euo pipefail', 'sed -i -e "s|/archive.ubuntu.com|/us-west-2.ec2.archive.ubuntu.com|g" /etc/apt/sources.list']
        cmd = AEGEA_JUNK + [cmd]

    else:
        quick_junk = ["/bin/bash", "-c", 'for i in "$@"; do eval "$i"; done; cd /', __name__]
        cmd = quick_junk + [cmd]
        #cmd = [cmd]
    # FULL_CMD="docker run -v {2}:/root/.aws/credentials -v {3}:/root/accessible_testing_data/ -v {4}:/root/accessible_results/ {0} /bin/bash -c \"{1}\"".format(image, cmd, cred_loc, test_loc, output_loc)
    # print(FULL_CMD)
    # call(FULL_CMD, shell=True)


    BASE_CMD = shlex.split("docker run -v {2}:/root/.aws/credentials -v {3}:/root/accessible_testing_data/ -v {4}:/root/accessible_results/ {0}".format(image, cmd, cred_loc, test_loc, output_loc, mt_loc))
    FULL_CMD = BASE_CMD + cmd
    print(FULL_CMD)
    print(' '.join(FULL_CMD))
    call(FULL_CMD)
    #FULL_CMD="docker run -v {2}:/root/.aws/credentials -v {3}:/root/accessible_testing_data/ -v {4}:/root/accessible_results/ {0} {1}".format(image, cmd, cred_loc, test_loc, output_loc)

class test_instrain():
    def setUp(self):
        self.test_dir = load_random_test_dir()
        self.IMAGE = "mattolm/instrain:latest"
        self.BAM_S3 = 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.SAM_S3 = 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sam'
        self.GENES_S3 = 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/N5_271_010G1_scaffold_min1000.fa.genes.fna'
        self.STB_S3 = 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/N5_271_010G1.maxbin2.stb'
        self.FASTA_S3 = 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/N5_271_010G1_scaffold_min1000.fa'
        self.GENOME_LIST = 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/genomelist.txt'
        self.IS_1 = 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.forRC.IS/'
        self.IS_2 = 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.forRC.IS/'
        self.IS_FIF = 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/is_locs.txt'
        self.scafflist = 's3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/test_data/scaffList.txt'

        self.tearDown()

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)
        clear_s3_results()

    def run(self):
        self.setUp()
        self.test0()
        self.tearDown()

        self.setUp()
        self.test1()
        self.tearDown()

        self.setUp()
        self.test2()
        self.tearDown()

        self.setUp()
        self.test3()
        self.tearDown()

        self.setUp()
        self.test4()
        self.tearDown()

        self.setUp()
        self.test5()
        self.tearDown()

        self.setUp()
        self.test6()
        self.tearDown()

        self.setUp()
        self.test7()
        self.tearDown()


    def test0(self):
        '''
        make sure dependencies are working; make sure the right version of inStrain is in there
        '''
        # Set up command

        CMD = "./prepare.sh; conda activate work; inStrain profile --version > version.txt; aws s3 cp version.txt {0}/"\
                .format(get_s3_results_folder())

        # Run command
        run_docker(self.IMAGE, CMD, simulate_aegea=True)

        # Upload results
        output_files = load_s3_results()
        reported_version = read_s3_file(output_files[0])
        correct_version = "inStrain version {0}".format(__version__)
        assert reported_version == correct_version, [correct_version, reported_version]

        # Estimate cost
        estimate_cost(None, get_s3_results_folder())

    def test1(self):
        '''
        Full on basic test
        '''
        # Set up command

        CMD = "./prepare.sh; conda activate work; pip install inStrain --upgrade; time ./run_instrain.py --bam {0} --fasta {1} --results_directory {2} --wd_name {3} --cmd_args='--skip_plot_generation'".format(self.BAM_S3, self.FASTA_S3, get_s3_results_folder(), 'test')
        # SETUP_CMD = "cp /root/accessible_testing_data/run_instrain.py /mnt/;"
        # CMD = SETUP_CMD + CMD

        # Run command
        run_docker(self.IMAGE, CMD, simulate_aegea=True)

        # Set up intended output
        OUTPUT = ['docker_log.log', 'log.log', 'test_genome_info.tsv', 'scaffold_2_mm_2_read_2_snvs.pickle']

        # Estimate cost
        dls = [self.BAM_S3, self.FASTA_S3, self.BAM_S3 + 'bai']
        estimate_cost(dls, get_s3_results_folder())

        output_files = load_s3_results()
        basenames = [os.path.basename(b) for b in output_files]
        for o in OUTPUT:
            have = o in basenames
            assert have, [o, basenames]

    def test2(self):
        '''
        Test with .sam, genes, and .stb file
        '''
        # Set up command

        CMD = "./prepare.sh; conda activate work; ls; ./run_instrain.py --bam {0} --fasta {1} --results_directory {2} --wd_name {3} --cmd_args='--skip_plot_generation' --genes {4} --stb {5}".format(self.SAM_S3, self.FASTA_S3, get_s3_results_folder(), 'test', self.GENES_S3, self.STB_S3)
        # SETUP_CMD = "cp /root/accessible_testing_data/run_instrain.py /mnt/;"
        # CMD = SETUP_CMD + CMD

        # Run command
        run_docker(self.IMAGE, CMD, simulate_aegea=True)

        # Set up intended output
        OUTPUT = ['docker_log.log', 'log.log', 'test_gene_info.tsv', 'scaffold_2_mm_2_read_2_snvs.pickle']

        # Estimate cost
        dls = [self.BAM_S3, self.FASTA_S3, self.BAM_S3 + 'bai', self.GENES_S3, self.STB_S3]
        estimate_cost(dls, get_s3_results_folder())

        output_files = load_s3_results()
        basenames = [os.path.basename(b) for b in output_files]
        for o in OUTPUT:
            have = o in basenames
            assert have, [o, basenames]

    def test3(self):
        '''
        Test the timeout functionality
        '''
        # Set up command

        CMD = "./prepare.sh; conda activate work; ls; time ./run_instrain.py --bam {0} --fasta {1} --results_directory {2} --wd_name {3} --timeout 5 --cmd_args='--skip_plot_generation'".format(self.BAM_S3, self.FASTA_S3, get_s3_results_folder(), 'test')
        # SETUP_CMD = "cp /root/accessible_testing_data/run_instrain.py /mnt/;"
        # CMD = SETUP_CMD + CMD

        # Run command
        run_docker(self.IMAGE, CMD, simulate_aegea=True)

        # Set up intended output
        OUTPUT = ['docker_log.log', 'log.log']#, 'test_genomeWide_scaffold_info.tsv', 'scaffold_2_mm_2_read_2_snvs.pickle']
        MISSING_OUT =  ['test_genomeWide_scaffold_info.tsv', 'scaffold_2_mm_2_read_2_snvs.pickle']

        # Estimate cost
        dls = [self.BAM_S3, self.FASTA_S3, self.BAM_S3 + 'bai']
        estimate_cost(dls, get_s3_results_folder())

        output_files = load_s3_results()
        basenames = [os.path.basename(b) for b in output_files]
        for o in OUTPUT:
            have = o in basenames
            assert have, [o, basenames]

        for o in MISSING_OUT:
            have = o in basenames
            assert not have, o

    def test4(self):
        '''
        Test quick_profile
        '''
        CMD = "./prepare.sh; conda activate work;  pip install instrain --upgrade --pre; inStrain -h; time ./run_instrain.py --bam {0} --fasta {1} --results_directory {2} --wd_name {3} --command quick_profile".format(self.BAM_S3, self.FASTA_S3, get_s3_results_folder(), 'test')
        #SETUP_CMD = "cp /root/accessible_testing_data/run_instrain.py /mnt/;"
        #CMD = SETUP_CMD + CMD

        # Run command
        run_docker(self.IMAGE, CMD, simulate_aegea='semi')

        # Set up intended output
        OUTPUT = ['docker_log.log', 'coverm_raw.tsv']#, 'test_genomeWide_scaffold_info.tsv', 'scaffold_2_mm_2_read_2_snvs.pickle']
        MISSING_OUT =  ['test_genomeWide_scaffold_info.tsv', 'scaffold_2_mm_2_read_2_snvs.pickle']

        # Estimate cost
        dls = [self.BAM_S3, self.FASTA_S3, self.BAM_S3 + 'bai']
        estimate_cost(dls, get_s3_results_folder())

        output_files = load_s3_results()
        basenames = [os.path.basename(b) for b in output_files]
        for o in OUTPUT:
            have = o in basenames
            assert have, [o, basenames]

        for o in MISSING_OUT:
            have = o in basenames
            assert not have, o

    def test5(self):
        '''
        Test beta version installation
        '''
        # Set up command

        CMD = "./prepare.sh; conda activate work; mkdir inStrain; git clone git://github.com/MrOlm/inStrain.git; cd inStrain; git checkout v1.3.0; pip install . --upgrade; cd ../; inStrain -h".format(get_s3_results_folder())

        # Run command
        run_docker(self.IMAGE, CMD, simulate_aegea=True)

        # Estimate cost
        estimate_cost(None, get_s3_results_folder())

    def test6(self):
        '''
        Test compare
        '''
        # Set up command
        # CMD = "./prepare.sh; conda activate work; cp /root/accessible_testing_data/run_instrain.py ./; ./run_instrain.py --IS {0} {1} --results_directory {2} --wd_name {3} --command compare --cmd_args='--store_mismatch_locations' --scaffolds {4}".format(
        #     self.IS_1, self.IS_2, get_s3_results_folder(), 'testRC', self.scafflist)
        CMD = "./prepare.sh; conda activate work; ./run_instrain.py --IS {0} {1} --results_directory {2} --wd_name {3} --command compare --cmd_args='--store_mismatch_locations' --scaffolds {4}".format(
            self.IS_1, self.IS_2, get_s3_results_folder(), 'testRC', self.scafflist)

        # Run command
        run_docker(self.IMAGE, CMD, simulate_aegea='semi')

        # Estimate cost
        estimate_cost(None, get_s3_results_folder())

        # Make sure done
        output_files = load_s3_results()
        basenames = [os.path.basename(b) for b in output_files]

        # Set up intended output
        OUTPUT = ['docker_log.log', 'testRC_comparisonsTable.tsv', 'testRC_pairwise_SNP_locations.tsv']

        for o in OUTPUT:
            have = o in basenames
            assert have, [o, basenames]

    def test7(self):
        '''
        Test compare using FOF (file of file
        '''
        # Set up command
        # CMD = "./prepare.sh; conda activate work; cp /root/accessible_testing_data/* ./; ./run_instrain.py --IS {0} --results_directory {2} --wd_name {3} --command compare --cmd_args='--store_mismatch_locations' --scaffolds {4}".format(
        #     self.IS_FIF, None, get_s3_results_folder(), 'testRC', self.scafflist)
        CMD = "./prepare.sh; conda activate work; ./run_instrain.py --IS {0} {1} --results_directory {2} --wd_name {3} --command compare --cmd_args='--store_mismatch_locations' --scaffolds {4}".format(
            self.IS_1, self.IS_2, get_s3_results_folder(), 'testRC', self.scafflist)

        # Run command
        run_docker(self.IMAGE, CMD, simulate_aegea='semi')

        # Estimate cost
        estimate_cost(None, get_s3_results_folder())

        # Make sure done
        output_files = load_s3_results()
        basenames = [os.path.basename(b) for b in output_files]

        # Set up intended output
        OUTPUT = ['docker_log.log', 'testRC_comparisonsTable.tsv', 'testRC_pairwise_SNP_locations.tsv']

        for o in OUTPUT:
            have = o in basenames
            assert have, [o, basenames]


if __name__ == '__main__':
    sync_test_data()

    test_instrain().run()

    print('everything is working swimmingly!')
