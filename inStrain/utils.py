#!/usr/bin/env python

import sys
import subprocess
import matplotlib
import shutil
import seaborn
import numpy
import Bio
import pandas
import pysam

def humanbytes(B):
   '''
   Return the given bytes as a human friendly KB, MB, GB, or TB string
   https://stackoverflow.com/questions/12523586/python-format-size-application-converting-b-to-kb-mb-gb-tb/37423778
   '''
   B = float(B)
   KB = float(1024)
   MB = float(KB ** 2) # 1,048,576
   GB = float(KB ** 3) # 1,073,741,824
   TB = float(KB ** 4) # 1,099,511,627,776

   if B < KB:
      return '{0} {1}'.format(B,'Bytes' if 0 == B > 1 else 'Byte')
   elif KB <= B < MB:
      return '{0:.2f} KB'.format(B/KB)
   elif MB <= B < GB:
      return '{0:.2f} MB'.format(B/MB)
   elif GB <= B < TB:
      return '{0:.2f} GB'.format(B/GB)
   elif TB <= B:
      return '{0:.2f} TB'.format(B/TB)

def gen_dependency_report():
   """
   Return a string listing a number of dependencies and their versions
   """
   PROGRAMS = ['samtools', 'coverm']
   MODULES = [matplotlib, seaborn, numpy, Bio, pysam, pandas]

   ret_string = f'{"$"*10} DEPENDENCY REPORT {"$"*10}\n'
   ret_string += f'PYTHON\nRunning python v{sys.version}\n'

   ret_string += '\nPROGRAM DEPENDENCIES\n'
   for program in PROGRAMS:
      loc, works, version = find_program(program)
      works_message = {True: 'all good', False: '! NOT WORKING !'}[works]
      ret_string += f'{program:.<30} {works_message:10} (version={version}) (location = {loc})\n'.format(program, works_message, loc)

   ret_string += '\nPYTHON DEPENDENCIES\n'
   for module in MODULES:
      v = str(module.__version__)
      l = str(module.__file__)
      s = f"{l.split('/')[-2]:.<30} all good   (version={v}) (location = {l})\n"
      ret_string += s

   ret_string += f'{"$"*40}\n'

   return ret_string


def find_program(dep):
   """
   return location of progrgam, works = True/False (based on calling the help), and the version
   """
   # find the location of the program
   loc = shutil.which(dep)

   # make sure the help on the program works
   works = False
   v = 'na'
   if loc != None:
      try:
         v = subprocess.check_output([loc, '--version'], stderr=subprocess.STDOUT)
         v = str(v, 'utf-8').strip().split('\n')[0]
         works = True
      except:
         pass

   return loc, works, v
