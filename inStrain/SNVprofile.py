#!/usr/bin/env python

import os
import sys
import json
import time
import h5py
import pickle
import logging
import argparse
import numpy as np
import pandas as pd

if __name__ != '__main__':
    from ._version import __version__
    import inStrain.profileUtilities

import inStrain.profileUtilities

class SNVprofile:
    '''
    This object holds the profile of a single .fasta / .bam pair
    '''

    firstLevels = ['output', 'raw_data', 'log', 'figures']

    def __init__(self, location):
        # Parse and store the location
        self.location = os.path.abspath(location)

        # Make the file structure
        self._make_fileStructure()

        # Initialize attributes
        self._initialize_attributes()

    def store(self, name, value, type, description):
        '''
        Store an attribute in the attributes.tsv file

        Args:
            name = name of attribute
            value = actual value of the attribute
            type = type of attribute
            description = one-line description of the attribute
        '''
        # Prepare for storage
        if type == 'value':
            adb = pd.DataFrame({'value':value, 'type':type, 'description':description}, index=[name])

        elif type == 'dictionary':
            fileloc = self._get_fileloc(name) + '.json'
            self._store_dictionary(value, fileloc)
            adb = pd.DataFrame({'value':fileloc, 'type':type, 'description':description}, index=[name])

        elif type == 'list':
            fileloc = self._get_fileloc(name) + '.txt'
            self._store_list(value, fileloc)
            adb = pd.DataFrame({'value':fileloc, 'type':type, 'description':description}, index=[name])

        elif type == 'numpy':
            fileloc = self._get_fileloc(name) + '.npz'
            self._store_numpy(value, fileloc)
            adb = pd.DataFrame({'value':fileloc, 'type':type, 'description':description}, index=[name])

        elif type == 'pandas':
            fileloc = self._get_fileloc(name) + '.csv.gz'
            self._store_pandas(value, fileloc)
            adb = pd.DataFrame({'value':fileloc, 'type':type, 'description':description}, index=[name])

        elif type == 'special':
            fileloc = self._store_special(value, name)
            adb = pd.DataFrame({'value':fileloc, 'type':type, 'description':description}, index=[name])

        elif type == 'pickle':
            fileloc = self._get_fileloc(name) + '.pickle'
            self._store_pickle(value, fileloc)
            adb = pd.DataFrame({'value':fileloc, 'type':type, 'description':description}, index=[name])

        else:
            print(logging.error('I dont know how to save a {0} type, so Im just going to pickle it'.format(type)))
            fileloc = self._get_fileloc(name) + '.pickle'
            self._store_pickle(value, fileloc)
            adb = pd.DataFrame({'value':fileloc, 'type':'pickle', 'description':description}, index=[name])

        # See if it's already in there
        Adb = self._get_attributes_file()
        if name in Adb.index:

            # Make sure it's the same type and description; otherwise crash
            for thing in ['type', 'description']:
                if (Adb.loc[name, thing] != eval(thing)):
                    logging.error('WILL NOT OVERWRITE {0}; {1} arent the same ({2} vs {3}))'.format(
                            name, thing, Adb.loc[name, thing], eval(thing)))
                    return

            # Overwrite
            Adb.at[name, 'value'] = adb['value'].tolist()[0]

        # Add the value if not
        else:
            Adb = Adb.append(adb)

        self._store_attributes_file(Adb)

    def get(self, name, **kwargs):
        '''
        Get an attribute from the attributes.tsv file
        '''
        Adb = self._get_attributes_file()
        if name not in Adb.index:
            #logging.error("{0} is not in this IS object".format(name))
            return None

        type = Adb.loc[name, 'type']

        if type == 'value':
            return Adb.loc[name, 'value']

        filename = self.get_location('raw_data') + '/' + os.path.basename(Adb.loc[name, 'value'])

        if type == 'dictionary':
            return self._load_dictionary(filename)

        elif type == 'list':
            return self._load_list(filename)

        elif type == 'numpy':
            return self._load_numpy(filename)

        elif type == 'pandas':
            return self._load_pandas(filename)

        elif type == 'special':
            return self._load_special(filename, name, **kwargs)

        elif type == 'pickle':
            return self._load_pickle(filename)

        else:
            logging.error('I dont know how to load a {0} type!'.format(type))

    def get_location(self, name):
        if name == 'output':
            return os.path.join(self.location, 'output/')
        if name == 'raw_data':
            return os.path.join(self.location, 'raw_data/')
        if name == 'log':
            return os.path.join(self.location, 'log/')
        if name == 'figures':
            # To maintain backwards compatibility
            loc = os.path.join(self.location, 'figures/')
            if not os.path.exists(loc):
                os.makedirs(loc)
            return loc
        else:
            logging.error('I dont know the location of {0}!'.format(name))

    def verify(self):
        '''
        Run a series of tests to make sure everything is kosher

        1) No duplicates in the attributes
        2) Everything in the attributes exists
        '''
        pass

    def get_nonredundant_RC_table(self):
        '''
        Get an RC table with just one line per scaffold
        '''
        scdb = self.get('comparisonsTable')
        if (scdb is None) or (len(scdb) == 0):
            return pd.DataFrame()
        else:
            return scdb.sort_values('mm')\
                    .drop_duplicates(subset=['scaffold', 'name1', 'name2'], keep='last')\
                    .sort_index().drop(columns=['mm'])

    def get_nonredundant_scaffold_table(self):
        '''
        Get a scaffold table with just one line per scaffold
        '''
        scdb = self.get('cumulative_scaffold_table')
        if (scdb is None) or (len(scdb) == 0):
            return pd.DataFrame()
        else:
            return scdb.sort_values('mm')\
                    .drop_duplicates(subset=['scaffold'], keep='last')\
                    .sort_index().drop(columns=['mm'])

    def get_nonredundant_snv_table(self):
        '''
        Get a SNP table with just one line per scaffold
        '''
        scdb = self.get('cumulative_snv_table')
        if (scdb is None) or (len(scdb) == 0):
            return pd.DataFrame()
        else:
            return scdb.sort_values('mm')\
                    .drop_duplicates(subset=['scaffold', 'position'], keep='last')\
                    .sort_index().drop(columns=['mm'])

    def get_nonredundant_linkage_table(self):
        '''
        Get a SNP table with just one line per scaffold
        '''
        scdb = self.get('raw_linkage_table')
        if (scdb is None) or (len(scdb) == 0):
            return pd.DataFrame()
        else:
            return scdb.sort_values('mm')\
                    .drop_duplicates(subset=['scaffold', 'position_A', 'position_B'], keep='last')\
                    .sort_index().drop(columns=['mm'])

    def get_clonality_table(self, nonredundant=True):
        '''
        Get a clonality table
        '''
        clonT = self.get('clonT')
        if clonT is None:
            return pd.DataFrame()

        else:
            dbs = []
            scaff2clonT = clonT
            for scaff, clonT in scaff2clonT.items():
                db = inStrain.profileUtilities._clonT_to_table(clonT)
                db['scaffold'] = scaff
                dbs.append(db)

            # The dropna is necessary because clonT doesn't have the "run_up_NaN"
            Cdb = pd.concat(dbs).dropna().reset_index(drop=True)

            if nonredundant:
                Cdb = Cdb.sort_values('mm').dropna()\
                        .drop_duplicates(subset=['scaffold', 'position'], keep='last')\
                        .sort_index().drop(columns=['mm'])

            return Cdb

    def __str__(self):
        '''
        String representation of attributes file
        '''
        string = str(self._get_attributes_file())
        return string

    def _make_fileStructure(self):
        '''
        Make the top level file structure
        '''
        location = self.location

        if not os.path.exists(location):
            os.makedirs(location)

        for l in SNVprofile.firstLevels:
            loc = location + '/' + l
            if not os.path.exists(loc):
                os.makedirs(loc)

    def _initialize_attributes(self):
        '''
        Make an "attributes.tsv" file, and store the location in
        '''
        # This is a new IS object
        if not os.path.exists(os.path.join(self.location, 'raw_data/attributes.tsv')):

            table = {'value':[], 'type':[], 'description':[]}
            Adb = pd.DataFrame(table)
            self._store_attributes_file(Adb)

            # Store the location
            self.store('location', self.location, 'value', 'Location of SNVprofile object')

            # Store the version
            self.store('version', __version__, 'value', 'Version of inStrain')

            # Make a README.txt file
            rm_loc = self._get_fileloc('_README.txt')
            with open(rm_loc, 'w') as o:
                o.write("The data in this folder can be easily accessed using the inStrain python API.\nFor information on how this is done, see the inStrain documentaion at https://instrain.readthedocs.io/en/latest/\n")


        # This is an existing IS object
        else:
            Adb = self._get_attributes_file()

            if not same_versions(Adb.loc['version', 'value'], __version__):
                logging.error("Warning! Your inStrain folder is from version {0}, while the installed version is {1}.\nIf you experience weird behavior, this might be why".format(
                        Adb.loc['version', 'value'], __version__))

            if self.location != self.get('location'):
                self.store('location', self.location, 'value', 'Location of SNVprofile object')

    def _get_attributes_file(self):
        '''
        Load the attributes file
        '''

        Aloc = str(os.path.join(self.location, 'raw_data/attributes.tsv')).strip()
        # This crazy mumbo-jumbo is to prevent a failure when trying to read while its being written
        for i in range(0,100):
            while True:
                try:
                    Adb = pd.read_csv(Aloc, sep='\t', index_col='name')
                except:
                    logging.error("Cannot load {0}; try {1}".format(Aloc, i))
                    time.sleep(1)
                break

        return Adb

    def _store_attributes_file(self, Adb):
        '''
        Store the attributes files
        '''
        Adb.to_csv(os.path.join(self.location, 'raw_data/attributes.tsv'),
                                    sep='\t', index_label='name')

    def _get_covt_keys(self):
        '''
        A special method for getting the keys to covT
        '''
        Adb = self._get_attributes_file()

        scaffs = set()
        filename = self.get_location('raw_data') + '/' + os.path.basename(Adb.loc['covT', 'value'])
        f = h5py.File(filename, 'r')
        for thing in list(f.keys()):
            scaff, mm = thing.split('::')
            scaffs.add(scaff)
        return scaffs


    # def _store_special(self, obj, name):
    #     '''
    #     store special things
    #     '''
    #     if name in ['covT', 'snpsCounted', 'clonT']:
    #         '''
    #         The idea here is to save a series of numpy arrays.
    #
    #         Each scaffold has two arrays in a row; one ending in ::covs and one
    #         ending in ::mms. The covs is a 2d array of mm and counts, and mms
    #         is just a list of the true mm levels, so that you can translate
    #         the index of ::covs back to the original
    #         '''
    #         fileloc = self._get_fileloc(name) + '.npz'
    #
    #         to_save = {}
    #         for scaff, covt in obj.items():
    #             if (scaff == 'N5_271_010G1_scaffold_0') & (name == 'snpsCounted'):
    #                 print('storing')
    #
    #             covs = []
    #             mms = []
    #             for mm, cov in covt.items():
    #                 mms.append(mm)
    #                 covs.append(cov)
    #
    #             to_save[scaff + '::covs'] = np.array(covs)
    #             to_save[scaff + '::mms'] = np.array(mms)
    #
    #         np.savez_compressed(fileloc, **to_save)
    #
    #     else:
    #         logging.error("I dont know how to store {0}! Ill just pickle it".format(name))
    #         fileloc = self._get_fileloc(name) + '.pickle'
    #         self._store_pickle(obj, fileloc)
    #
    #     return fileloc
    #
    def _load_special_old(self, location, name):
        '''
        store special things; this is for loading deprecated things
        '''
        if name in ['covT', 'snpsCounted', 'clonT', 'clonTR']:
            base = np.load(location)
            files = base.files


            covT = {}
            # iterate in groups of 2
            for covs, mms in [files[i:i+2] for i in range(0, len(files), 2)]:
                scaff = covs[:-6]
                assert covs[:-6] == mms[:-5], [covs[:-6],  mms[:-5]]

                covT[scaff] = {}
                mm_list = base[mms]
                cov_array = base[covs]

                for i, mm in enumerate(mm_list):
                    covT[scaff][mm] = cov_array[i]

            return covT

        else:
            logging.error("I dont know how to load {0}! Ill just try pickle ".format(name))
            return self._load_pickle(location)

    def _store_special(self, obj, name):
        '''
        store special things using hd5
        '''
        if name in ['covT', 'snpsCounted', 'clonT', 'clonTR']:
            '''
            The idea here is to save the series'
            '''
            fileloc = self._get_fileloc(name) + '.hd5'
            f = h5py.File(fileloc, "w")
            for scaff, clon in obj.items():
                for mm, arr in clon.items():
                    #arr.to_hdf(fileloc, key="{0}::{1}".format(scaff, mm), mode='a')
                    dset = f.create_dataset("{0}::{1}".format(scaff, mm),
                            data=np.array([arr.values, arr.index]),
                            compression="gzip")# convert from series to 2d array

        elif name in ['scaff2pair2mm2SNPs', 'scaff2pair2mm2cov']:
            fileloc = self._get_fileloc(name) + '.hd5'
            if name == 'scaff2pair2mm2cov':
                convert = True
            else:
                convert = False
            store_scaff2pair2mm2SNPs(obj, fileloc, convert=convert)


        else:
            logging.error("I dont know how to store {0}! Ill just pickle it".format(name))
            fileloc = self._get_fileloc(name) + '.pickle'
            self._store_pickle(obj, fileloc)

        return fileloc

    def _load_special(self, location, name, **kwargs):
        '''
        store special things

        with the kwarg "scaffolds", only load certain scaffolds
        '''
        scaffolds = kwargs.get('scaffolds', [])
        pairs = kwargs.get('pairs', [])

        if name in ['covT', 'snpsCounted', 'clonT', 'clonTR']:
            if location[-4:] == '.npz':
                return self._load_special_old(location, name)

            scaff2mm = {}
            f = h5py.File(location, 'r')
            for thing in list(f.keys()):
                scaff, mm = thing.split('::')
                if scaffolds != []:
                    if scaff not in scaffolds:
                        continue

                #pd.read_hdf(location, key=thing)
                dset = list(f[thing])
                mm = int(mm)

                if scaff not in scaff2mm:
                    scaff2mm[scaff] = {}

                scaff2mm[scaff][mm] = pd.Series(data = dset[0],
                                                index = np.array(dset[1].astype('int'))) # convert from 2d array to series

            return scaff2mm

        elif name in ['scaff2pair2mm2SNPs', 'scaff2pair2mm2cov']:
            return load_scaff2pair2mm2SNPs(location, scaffolds=scaffolds, pairs=pairs)

        else:
            logging.error("I dont know how to load {0}! Ill just try pickle ".format(name))
            return self._load_pickle(location)


    def _store_list(self, obj, location):
        '''
        Store a list using txt
        '''
        assert type(obj) == type(list())

        with open(location, "w") as f:
            for s in obj:
                f.write(str(s) +"\n")

    def _load_list(self, location):
        '''
        Load a list using txt
        '''
        l  = []
        with open(location, "r") as f:
            for line in f:
                l.append(line.strip())
        return l

    def _store_numpy(self, obj, location):
        np.savez_compressed(location, obj)

    def _load_numpy(self, location):
        loaded = np.load(location, allow_pickle=True)['arr_0']
        return loaded

    def _store_pandas(self, obj, location):
        obj.to_csv(location)

    def _load_pandas(self, location):
        return pd.read_csv(location, index_col=0)

    def _store_pickle(self, obj, location):
        f = open(location, 'wb')
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    def _load_pickle(self, location):
        f = open(location, 'rb')
        tmp_dict = pickle.load(f)
        f.close()
        return tmp_dict

    def _store_dictionary(self, dic, location):
        '''
        Store a dictionary using json
        '''
        assert type(dic) == type({})

        with open(location, 'w') as fp:
            json.dump(dic, fp, default=_json_helper)

    def _load_dictionary(self, location):
        '''
        Store a dictionary using json
        '''
        with open(location, 'r') as fp:
            return(json.load(fp))

    def _get_fileloc(self, name):
        return os.path.join(self.location, 'raw_data/{0}'.format(name))

def same_versions(v1, v2):
    '''
    If the major or minor version are the same, return True; else, False
    '''
    minor1 = int(v1.split('.')[1])
    minor2 = int(v2.split('.')[1])

    major1 = int(v1.split('.')[0])
    major2 = int(v2.split('.')[0])

    return (minor1 == minor2) & (major1 == major2)

def store_scaff2pair2mm2SNPs(obj, fileloc, convert=False):
    '''
    Store as an hdf5 object

    If convert, declare it as a numpy array from a set
    '''
    f = h5py.File(fileloc, "w")
    for scaff, pair2mm2SNPs in obj.items():
        for pair, mm2SNPs in pair2mm2SNPs.items():
            for mm, arr in mm2SNPs.items():
                if convert:
                    dset = f.create_dataset("{0}::{1}::{2}".format(scaff, pair, mm),
                            data=np.fromiter(arr, int, len(arr)), compression="gzip")
                else:
                    dset = f.create_dataset("{0}::{1}::{2}".format(scaff, pair, mm),
                            data=np.array(arr),compression="gzip")

def load_scaff2pair2mm2SNPs(location, scaffolds=[], pairs=[]):
    scaff2pair2mm2SNPs = {}
    f = h5py.File(location, 'r')

    for thing in list(f.keys()):
        scaff, pair, mm = thing.split('::')

        if scaffolds != []:
            if scaff not in scaffolds:
                continue

        if pairs != []:
            if pair not in pairs:
                continue

        dset = f[thing].value
        mm = int(mm)

        if scaff not in scaff2pair2mm2SNPs:
            scaff2pair2mm2SNPs[scaff] = {}

        if pair not in scaff2pair2mm2SNPs[scaff]:
            scaff2pair2mm2SNPs[scaff][pair] = {}

        scaff2pair2mm2SNPs[scaff][pair][mm] = dset # convert from 2d array to series

    return scaff2pair2mm2SNPs

def _json_helper(o):
    if isinstance(o, np.int64): return int(o)
    raise TypeError

class SNVprofile_old:
    '''
    The class holds the profile of a single .fasta / .bam pair
    '''

    def __init__(self, **kwargs):
        '''
        initialize all attributes to None
        '''

        self.ATTRIBUTES = (\
        'filename', # Filename of this object

        'fasta_loc',
        'scaffold2length', # Dictionary of scaffold 2 length
        'bam_loc',

        'scaffold_list', # 1d list of scaffolds, in same order as counts_table
        'counts_table', #1d numpy array of 2D counts tables for each scaffold
        'raw_snp_table', # Contains raw SNP information on a mm level
        'raw_ANI_table', # Contains raw ANI information on a mm level
        'raw_coverage_table', # Contains raw coverage information on a mm level
        'raw_linkage_table', # Contains raw linkage information on a mm level

        'cumulative_scaffold_table', # Cumulative coverage on mm level. Formerly "scaffoldTable.csv"
        'cumulative_snv_table', # Cumulative SNP on mm level. Formerly "snpLocations.pickle"

        # THE FOLLOWING ARE HEAVY-RAM OBJECTS
        'scaff2covT',
        'scaff2basesCounted',
        'scaff2snpsCounted',
        )
        for att in self.ATTRIBUTES:
            setattr(self, att, kwargs.get(att, None))
        self.version = __version__

    def __str__(self):
        string = '\n'.join(["{0} - {1}".format(att, type(getattr(self, att)))\
            for att in self.ATTRIBUTES])
        return string

    def get_nonredundant_scaffold_table(self):
        '''
        Get a scaffold table with just one line per scaffold
        '''
        if (self.cumulative_scaffold_table is None) or (len(self.cumulative_scaffold_table) == 0):
            return pd.DataFrame()
        else:
            return self.cumulative_scaffold_table.sort_values('mm')\
                    .drop_duplicates(subset=['scaffold'], keep='last')\
                    .sort_index().drop(columns=['mm'])

    def get_nonredundant_linkage_table(self):
        '''
        Get a SNP table with just one line per scaffold
        '''
        if (self.raw_linkage_table is None) or (len(self.raw_linkage_table) == 0):
            return pd.DataFrame()
        else:
            return self.raw_linkage_table.sort_values('mm')\
                    .drop_duplicates(subset=['scaffold', 'position_A', 'position_B'], keep='last')\
                    .sort_index().drop(columns=['mm'])

    def get_nonredundant_snv_table(self):
        '''
        Get a SNP table with just one line per scaffold
        '''
        if (self.cumulative_snv_table is None) or (len(self.cumulative_snv_table) == 0):
            return pd.DataFrame()
        else:
            return self.cumulative_snv_table.sort_values('mm')\
                    .drop_duplicates(subset=['scaffold', 'position'], keep='last')\
                    .sort_index().drop(columns=['mm'])

    def get_clonality_table(self, nonredundant=True):
        '''
        Get a clonality table
        '''
        if not hasattr(self, 'clonT'):
            return pd.DataFrame()

        elif self.clonT is None:
            return pd.DataFrame()

        else:
            dbs = []
            scaff2clonT = self.clonT
            for scaff, clonT in scaff2clonT.items():
                db = inStrain.profileUtilities._clonT_to_table(clonT)
                db['scaffold'] = scaff
                dbs.append(db)

            # The dropna is necessary because clonT doesn't have the "run_up_NaN"
            Cdb = pd.concat(dbs).dropna().reset_index(drop=True)

            if nonredundant:
                Cdb = Cdb.sort_values('mm').dropna()\
                        .drop_duplicates(subset=['scaffold', 'position'], keep='last')\
                        .sort_index().drop(columns=['mm'])

            return Cdb

    def store(self):
        '''
        Store self. MUST have the attribute "filename" set
        '''
        if self.filename is None:
            print("Cant save this SNVprofile- no filename!")
            return
        elif self.filename[-7:] == '.pickle':
            self.filename = self.filename[:-7]

        f = open(self.filename + ".pickle", 'wb')
        pickle.dump(self.__dict__, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    def load(self, basename):
        '''
        Load self from the basename
        '''
        if basename[-7:] == '.pickle':
            basename = basename[:-7]

        #print("loading {0}".format(basename + '.pickle'))
        f = open(basename + '.pickle', 'rb')
        tmp_dict = pickle.load(f)
        f.close()

        self.__dict__.clear()
        self.__dict__.update(tmp_dict)

        return self

def convert_SNVprofile(pickle_loc):
    '''
    From the old SNVprofile type (a pickle), make a new one
    '''
    # Load the old version
    oIS = SNVprofile_old()
    oIS.load(pickle_loc)

    # Get a new version set up
    nIS = SNVprofile(pickle_loc.replace('.pickle', ''))

    for attr in oIS.__dict__.keys():
    #for attr in oIS.ATTRIBUTES:
        #print('storing {0}'.format(attr))
        if attr in ['ATTRIBUTES', 'filename']:
            continue

        elif getattr(oIS, attr) is None:
            pass

        elif attr == 'version':
            nIS.store('old_version', getattr(oIS, attr), 'value', 'Version of old SNVprofile this was made from')

        elif attr == 'fasta_loc':
            # THIS IS INCORRECTLY SAVED AT TH MOMEMENT! - M.O. 6/6/19
            #nIS.store('fasta_loc', oIS.fasta_loc, 'value', 'Location of .fasta file')
            nIS.store('fasta_loc', 'unk', 'value', 'Location of .fasta file')

        elif attr == 'scaffold2length':
            nIS.store('scaffold2length', oIS.scaffold2length, 'dictionary', 'Dictionary of scaffold 2 length')

        elif attr == 'bam_loc':
            nIS.store('bam_loc', getattr(oIS, attr), 'value', 'Location of .bam file')

        elif attr == 'scaffold_list':
            nIS.store('scaffold_list', getattr(oIS, attr), 'list', '1d list of scaffolds, in same order as counts_table')

        elif attr == 'counts_table':
            nIS.store('counts_table', getattr(oIS, attr), 'numpy', '1d numpy array of 2D counts tables for each scaffold')

        elif attr == 'raw_snp_table':
            nIS.store('raw_snp_table', getattr(oIS, attr), 'pandas', 'Contains raw SNP information on a mm level')

        elif attr == 'raw_ANI_table':
            assert getattr(oIS, attr) == None
            #nIS.store('raw_ANI_table', getattr(oIS, attr), 'pandas', 'Contains raw ANI information on a mm level')

        elif attr == 'raw_coverage_table':
            assert getattr(oIS, attr) == None

        elif attr == 'raw_linkage_table':
            nIS.store('raw_linkage_table', getattr(oIS, attr), 'pandas', 'Contains raw linkage information on a mm level')

        elif attr == 'cumulative_scaffold_table':
            nIS.store('cumulative_scaffold_table', getattr(oIS, attr), 'pandas', "Cumulative coverage on mm level. Formerly scaffoldTable.csv")

        elif attr == 'cumulative_snv_table':
            nIS.store('cumulative_snv_table', getattr(oIS, attr), 'pandas', "Cumulative SNP on mm level. Formerly snpLocations.pickle")
        #
        # '''
        # THOSE DON'T WORK BECAUSE THE BASEWISE SHRINK ISN'T MEANT TO WORK WHEN YOU HAVE A SCAFFOLD FISRT!!
        # '''
        #

        elif attr == 'covT':
            new = {s:inStrain.profileUtilities.shrink_basewise(cov, 'coverage') for s, cov in getattr(oIS, attr).items()}
            nIS.store('covT', new,
                    'special', "Scaffold -> mm -> position based coverage")
            # nIS.store('covT', inStrain.profileUtilities.shrink_basewise(getattr(oIS, attr), 'coverage'),
            #         'special', "Scaffold -> mm -> position based coverage")

        elif attr == 'snpsCounted':
            new = {s:inStrain.profileUtilities.shrink_basewise(cov, 'snpCounted') for s, cov in getattr(oIS, attr).items()}
            nIS.store('snpsCounted', new,
                    'special', "Scaffold -> mm -> position based True/False on if a SNPs is there")
            # nIS.store('snpsCounted', inStrain.profileUtilities.shrink_basewise(getattr(oIS, attr), 'snpsCounted'),
            #         'special', "Scaffold -> mm -> position based True/False on if a SNPs is there")

        elif attr == 'clonT':
            new = {s:inStrain.profileUtilities.shrink_basewise(cov, 'clonality') for s, cov in getattr(oIS, attr).items()}
            nIS.store('clonT', new,
                    'special', "Scaffold -> mm -> position based clonality")

        elif attr == 'clonTR':
            new = {s:inStrain.profileUtilities.shrink_basewise(cov, 'clonality') for s, cov in getattr(oIS, attr).items()}
            nIS.store('clonTR', new,
                    'special', "Scaffold -> mm -> position based clonality")
            # nIS.store('clonT', inStrain.profileUtilities.shrink_basewise(getattr(oIS, attr), 'clonality'),
            #         'special', "Scaffold -> mm -> position based clonality")

        elif attr == 'read_report':
            nIS.store('read_report', getattr(oIS, attr), 'pandas', "Report on reads")

        else:
            logging.error('I dont know how to store {0}!'.format(attr))
            print(type(getattr(oIS, attr)))
            break

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description= """
#         A quick way to convert old SNP profiles (v0.3.x) to new SNP profiles (v0.4.x)""")
#
#     # Required positional arguments
#     parser.add_argument('input', help="an on inStrain object pickle")
#     args = parser.parse_args()
#
#     from _version import __version__
#
#     convert_SNVprofile(args.input)
