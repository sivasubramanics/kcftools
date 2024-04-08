!#/usr/bin/env python3

"""
This module is used to parse the command line arguments and return specific to the plugin
"""


class ArgumentParser:
    def __init__(self, args):
        self.command = args.command
        self.inputs = parse_args(args)

    def parse_args(args):
        if self.command == "