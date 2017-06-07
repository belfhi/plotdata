#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import pencil as pc
from scipy.optimize import curve_fit
from helperfuncs import *
import matplotlib as mpl
mpl.use('pgf')
import sys
from mpl_toolkits.axes_grid1 import Grid
import argparse
from os import mkdir, path, listdir
from os.path import join, isdir, isfile
from fractions import Decimal, Fraction
from scipy.integrate import simps
