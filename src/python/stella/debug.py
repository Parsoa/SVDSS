from __future__ import print_function

import io
import os
import re
import pwd
import sys
import copy
import json
import time
import argparse
import operator
import traceback

from stella import (
    config,
)

from stella.logger import *
print = pretty_print

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

def debug_breakpoint():
    c = config.Configuration()
    if c.debug:
        print(magenta('BEGIN DEBUG *********************************************************************************************'))
        print(magenta('Debug Breakpoint. Press Enter to continue...'))
        print(magenta('END DEBUG ***********************************************************************************************'))
        s = input()

def debug_log(*args):
    c = config.Configuration()
    if c.debug:
        print(magenta('BEGIN DEBUG *********************************************************************************************'))
        print(*args)
        print(magenta('END DEBUG ***********************************************************************************************'))

def debug_terminate():
    c = config.Configuration()
    if c.debug:
        print(magenta('DEBUG TERMINATE *********************************************************************************************'))
        exit()

def debug_print(*args):
    c = config.Configuration()
    if c.debug:
        print(*args)
