#!/usr/bin/env python3
# coding: utf-8

import re
import fileinput
import argparse
from math import ceil

parser = argparse.ArgumentParser()
parser.add_argument('choice', choices=['run', 'tspec'], default=None, help='switch for run.in or tspec.dat')
parser.add_argument('-c', '--comment', nargs='+', default=[], type=int, help='linenumbers to comment out')
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='adds verbosity output')

args = parser.parse_args()
if args.verbose:
    for arg in vars(args):
        print(arg, getattr(args, arg))

if args.choice == 'run':
    rfile = 'run.in'
elif args.choice == 'tspec':
    rfile = 'nonhel_k80/tspec.dat'

clines = []
if args.comment:
    clines = [l+1 for l in args.comment]

for l, line in enumerate(fileinput.input(rfile, inplace=1, backup='.bak')):
    if l in args.comment:
        line = '!' + line.rstrip()
    else:
        line = re.sub('tmax=', 'tmax=100', line.rstrip())
        line = re.sub('dspec=[0-9]*\.[0-9]*', 'dspec=1.0', line.rstrip())
    if args.choice == 'tspec' and l == 0:
        line = '   {0:0<8.1f}               {1}'.format(ceil(float(line.split()[0])),line.split()[1])
    else:
        line = line.rstrip()
    print(line)
