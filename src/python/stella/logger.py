import os
import json
import time
import functools

import colorama

colorama.init()

# ============================================================================================================================ #
# STDIO Wrappers/Helpers
# ============================================================================================================================ #

LOG_LEVEL = 1

def colorize(*vargs):
    s = ''.join(functools.reduce(lambda x, y: x + str(y) + ' ', vargs))
    return s.strip() + colorama.Fore.WHITE

def red(*args):
    return colorize(colorama.Fore.RED, *args)

def cyan(*args):
    return colorize(colorama.Fore.CYAN, *args)

def blue(*args):
    return colorize(colorama.Fore.BLUE, *args)

def white(*args):
    return colorize(colorama.Fore.WHITE, *args)

def green(*args):
    return colorize(colorama.Fore.GREEN, *args)

def yellow(*args):
    return colorize(colorama.Fore.YELLOW, *args)

def magenta(*args):
    return colorize(colorama.Fore.MAGENTA, *args)

def pretty_print(*args):
    def inner(*vargs):
        return ''.join(functools.reduce(lambda x, y: x + str(y) +  colorama.Fore.WHITE + ' ', vargs))
    print(inner(colorama.Fore.WHITE, *args))

def user_print(*args):
    print(white(*args))

def user_print_error(*args):
    print(red(*args))

def user_print_warning(*args):
    print(yellow(*args))

def system_print(*args):
    if 0 >= LOG_LEVEL:
        print(white(*args))

def system_print_high(*args):
    if 1 >= LOG_LEVEL:
        print(white(*args))

def system_print_warning(*args):
    if 2 >= LOG_LEVEL:
        print(yellow(*args))

def system_print_error(*args):
    if 3 >= LOG_LEVEL:
        print(red(*args))

def json_print(d):
    print(json.dumps(d, sort_keys = True, indent = 4, separators = (',', ': ')))

def jsonify(s):
    return json.dumps(s, sort_keys = True, indent = 4, separators = (',', ': '))

def compactify(l):
    return '[' + ','.join(map(lambda t: str(t), l)) + ']'
