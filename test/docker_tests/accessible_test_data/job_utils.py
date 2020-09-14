import os
import sys
import shutil
import uuid
import logging


def generate_working_dir(working_dir_base):
    """
    Creates a unique working directory to combat job multitenancy
    :param working_dir_base: base working directory
    :return: a unique subfolder in working_dir_base with a uuid
    """

    working_dir = os.path.join(working_dir_base, str(uuid.uuid4()))
    try:
        os.makedirs(working_dir, exist_ok=True)
    except Exception as e:
        print(e)
        return working_dir_base
    return working_dir


def delete_working_dir(working_dir):
    """
    Deletes working directory
    :param working_dir:  working directory
    """

    try:
        shutil.rmtree(working_dir)
    except Exception as e:
        print ('Can\'t delete %s' % working_dir)

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
