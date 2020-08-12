'''
This file contains the API-like elements that are called from the outside
'''

import inStrain.profile.profile_controller

def profile_bam(bam, Fdb, sR2M, **kwargs):
    '''
    Profile a bam file with inStrain

    arguments:
        bam = location of .bam file
        Fdb = dictionary listing fasta locations to profile
        sR2M = dictionary of scaffold -> read pair -> number of mm
    '''
    return inStrain.profile.profile_controller.BamProfileController(
                                                bam, Fdb, sR2M, **kwargs).main()
