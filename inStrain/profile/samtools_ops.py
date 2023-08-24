"""
Methods that make external system calls to samtools
"""

import os
import sys
import pysam
import logging
from subprocess import call
from Bio import bgzf
from subprocess import PIPE, run, Popen

def prepare_bam_fie(bam, processes=6, max_reads=None):
    """
    Make this a .bam file

    Args:
        bam = location of .bam file
    """
    if bam[-4:] == '.sam':
        logging.info("You gave me a sam- I'm going to make it a .bam now")

        bam = _sam_to_bam(bam, processes=processes)
        bam = _sort_index_bam(bam, processes=processes)

    elif bam[-4:] == '.bam':
        # If there's an index, assume its sorted
        if (os.path.exists(bam + '.bai')) | (os.path.exists(bam[:-4] + '.bai')):
            pass

        # If not sorted, sort and index
        elif not _is_sorted_bam(bam):
            bam = _sort_index_bam(bam, processes=processes, rm_ori=False)

        # If sorted, just index
        else:
            logging.info("Indexing {0}".format(bam))
            _index_bam(bam, processes=processes)

    if max_reads is not None:
        total = calc_number_reads(bam)
        ratio = int(max_reads) / total

        if ratio < 1:
            seed = str(42) + f"{ratio:.3f}"[1:]

            new_bam = bam + '.sub'

            logging.info(f"Subsetting with seed {seed}")
            bam = subset_bam(bam, seed, new_bam, p=processes)

            _index_bam(bam, processes=processes)

    if os.stat(bam).st_size == 0:
        logging.error("Failed to generated a sorted .bam file! Make sure you have " +
                      "samtools version 1.6 or greater.")
        sys.exit()

    # Do a quick sanity check
    try:
        pysam.AlignmentFile(bam).mapped
    except ValueError:
        logging.error("It seems like the .bam file could not be indexed!" +
                      "Make sure you have samtools version 1.6 or greater.")
        sys.exit()

    return bam


def _sam_to_bam(sam, processes=6):
    """
    From the location of a .sam file, convert it to a bam file and return the location
    """
    if sam[-4:] != '.sam':
        print('Sam file needs to end in .sam')
        sys.exit()

    bam = sam[:-4] + '.bam'
    logging.info("Converting {0} to {1}".format(sam, bam))
    cmd = ['samtools', 'view', '-S',  '-@', str(processes), '-b', sam, '>', bam]
    print(' '.join(cmd))
    call(' '.join(cmd), shell=True)

    return bam

def _index_bam(bam, processes=6):
    # Try standard way
    cmd = ['samtools', 'index', bam, bam + '.bai', '-@', str(processes)]
    print(' '.join(cmd))
    call(cmd)
    exit_code = call(cmd)

    if exit_code != 0:
        # Try non-standard way
        cmd = ['samtools',  'index', '-@', str(processes), bam]
        print(' '.join(cmd))
        call(cmd)

def _sort_index_bam(bam, processes=6, rm_ori=False):
    """
    From a .bam file, sort and index it. Remove original if rm_ori
    Return path of sorted and indexed bam
    """
    if bam[-4:] != '.bam':
        logging.error('Bam file needs to end in .bam')
        sys.exit()

    if 'sorted.bam' not in bam:
        logging.info("sorting {0}".format(bam))
        sorted_bam = bam[:-4] + '.sorted.bam'
        cmd = ['samtools', 'sort', bam, '-o', sorted_bam, '-@', str(processes)]
        print(' '.join(cmd))
        call(cmd)

    else:
        sorted_bam = bam
        rm_ori = False

    logging.info("Indexing {0}".format(sorted_bam))
    _index_bam(sorted_bam, processes=processes)


    if rm_ori:
        logging.info("Deleting {0}".format(bam))
        os.remove(bam)

    return sorted_bam


def _is_sorted_bam(bam):
    """
    Checks if a BAM file is sorted by coordinate.
    """
    with bgzf.BgzfReader(bam, "rb") as fin:
        bam_header = fin.readline().strip()
        return b"SO:coordinate" in bam_header

def calc_number_reads(bam):
    """
    Return the number of reads in the .bam file
    """
    cmd = ['samtools', 'idxstats', bam]
    logging.info(' '.join(cmd))

    # run idxstats
    result = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    assert result.returncode == 0, result

    # sum number of reads
    total_reads = 0
    for line in result.stdout.split('\n'):
        if len(line.split()) != 4:
            continue
        total_reads += (int(line.split()[2]) + int(line.split()[3]))

    return total_reads

def subset_bam(ori_bam, substr, new_bam, p=6):
    cmd = f"samtools view -@ {p} -bs {substr} {ori_bam} > {new_bam}"
    #cmd = f"sambamba view -s {substr} -t {p} -f bam {ori_bam} > {new_bam}"
    logging.info(cmd)
    result = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
    assert result.returncode == 0, result
    return new_bam