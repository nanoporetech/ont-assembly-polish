#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from collections import OrderedDict

from wub.util import parse as parse_util
from wub.util import misc
from wub.vis import report

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Plot accuracies as bar plot.')
parser.add_argument(
    '-i', metavar='input', type=str, help="Input specifier string: tag:file,....", default="")
parser.add_argument(
    '-t', metavar='type', type=str, help="Pickles type: dnadiff or lastal", choices=('dnadiff', 'lastal'))
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF.", default="plot_accuracies.pdf")


def get_accuracy_dnadiff(res_file):
    """Get accuracy from dandiff pickle file."""
    res = misc.pickle_load(res_file)
    props = res['Alignments']['1-to-1']['AvgIdentity']
    if props.ref != props.query:
        raise Exception('Mismatch in calcuated dnadiff identity')
    return props.ref


def get_accuracy_lastal(res_file):
    """Get accuracy from dandiff pickle file."""
    res = misc.pickle_load(res_file)
    return res['Accuracy']


if __name__ == '__main__':
    args = parser.parse_args()

    plotter = report.Report(args.r)
    input_pickles = parse_util.args_string__to_dict(args.i)

    if args.t == 'dnadiff':
        get_accuracy = get_accuracy_dnadiff
    elif args.t == 'lastal':
        get_accuracy = get_accuracy_lastal

    accuracies = OrderedDict((tag, get_accuracy(res_file))
                             for tag, res_file in input_pickles.iteritems())
    plotter.plot_bars_simple(accuracies, title="Accuracies as measured by {}".format(args.t), xlab="Method", ylab="Accuracy", auto_limit=True)

    plotter.close()
