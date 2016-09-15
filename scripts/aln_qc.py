#!/usr/bin/env python

import argparse
import sys
import os
import pysam

from collections import defaultdict
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # import seaborn as sns
warnings.resetwarnings()


class Report:

    """
    Wrapper class for plotting figures using matplotlib.
    """

    def __init__(self, pdf):
        self.pdf = pdf
        self.pages = PdfPages(pdf)

    def plot_arrays(self, a, b, title="", xlab="", ylab=""):
        """
        Dot plot of two arrays.
        """
        fig = plt.figure()
        plt.plot(a, b, 'b.')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)

    def plot_heatmap(self, z, title="", xlab="", ylab=""):
        """
        Plot heatmap.
        """
        fig = plt.figure()
        p = plt.contourf(z)
        plt.colorbar(p, orientation='vertical')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        plt.close(fig)
        self.pages.savefig(fig)

    def plot_hash_line(self, d, title="", xlab="", ylab=""):
        """
        Plot heatmap.
        """
        fig = plt.figure()

        plt.plot(d.keys(), d.values(), '-')

        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        plt.close(fig)
        self.pages.savefig(fig)

    def plot_heatmap(self, z, title="", xlab="", ylab=""):
        """
        Plot heatmap.
        """
        fig = plt.figure()
        p = plt.contourf(z)
        plt.colorbar(p, orientation='vertical')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        plt.close(fig)
        self.pages.savefig(fig)

    def plot_hists(self, d1, d2, l1="", l2="", title="", xlab="", ylab="", bins1=50, bins2=50):
        """
        Plot two histograms of data.
        """
        fig = plt.figure()
        plt.hist(d1, bins=bins1, label=l1, alpha=0.7)
        plt.hist(d2, bins=bins2, label=l2, alpha=0.7)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.legend(loc='upper right')
        plt.title(title)
        self.pages.savefig(fig)

    def close(self):
        """
        Close backend PDF.
        """
        self.pages.close()


def mean_qscore(scores):
    """ Returns the phred score corresponding to the mean of the probabilities
    associated with the phred scores provided.
    :param scores: Iterable of phred scores.
    :returns: Phred score corresponding to the average error rate, as
        estimated from the input phred scores.
    """
    if len(scores) == 0:
        return 0.0
    sum_prob = 0.0
    for val in scores:
        sum_prob += pow(10, -0.1 * val)
    mean_prob = sum_prob / len(scores)
    return -10.0 * np.log10(mean_prob)


def parse_pileups(bam):
    st = defaultdict(lambda: defaultdict(list))
    samfile = pysam.AlignmentFile(bam, "rb")
    for pileupcolumn in samfile.pileup():
        # print pileupcolumn.reference_name, pileupcolumn.reference_pos, pileupcolumn.nsegments
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # print pileupcolumn.reference_name, pileupcolumn.reference_pos,
                # pileupread.alignment.query_qualities[pileupread.query_position]
                st[pileupcolumn.reference_name][pileupcolumn.reference_pos].append(
                    pileupread.alignment.query_qualities[pileupread.query_position])
    samfile.close()
    return st


def find_max_pos(d):
    return max(d.keys())


def find_max_qual(d):
    return max(itertools.chain.from_iterable(d.itervalues()))


def mean_qual_per_pos(d):
    mq = {}
    for pos, quals in d.iteritems():
        mq[pos] = mean_qscore(quals)
    return mq


def stats_to_matrix(rst):
    max_pos = find_max_pos(rst)
    max_qual = find_max_qual(rst)
    mat = np.zeros((max_qual + 1, max_pos + 1), dtype=float)
    for pos, quals in rst.iteritems():
        for q in quals:
            mat[q][pos] += 1
    return mat


def parse_reads(bam):
    res = defaultdict(list)
    bf = pysam.AlignmentFile(bam, "rb")
    for r in bf.fetch(until_eof=True):
        if r.is_unmapped:
            res['unaligned_quals'].append(mean_qscore(r.query_qualities))
            res['unaligned_lengths'].append(len(r.query_qualities))
        else:
            res['aligned_quals'].append(mean_qscore(r.query_qualities))
            res['alignment_lengths'].append(r.query_alignment_length)
            res['aligned_lengths'].append(len(r.query_qualities))
    return res


def ref_qual_qc(st, report):
    for ref, stats in st.iteritems():
        mat = stats_to_matrix(stats)
        report.plot_heatmap(mat, title="Quality values across {}".format(
            ref), xlab="Position", ylab="Quality bin")
        mq = mean_qual_per_pos(stats)
        report.plot_hash_line(mq, title="Mean quality values across {}".format(
            ref), xlab="Position", ylab="Mean quality")


def read_qual_qc(st, report):
    report.plot_hists(st['unaligned_quals'], st['aligned_quals'], l1="Unaligned", l2="Aligned",
                      title="Distribution of mean quality values in fractions", xlab="Mean base quality", ylab="Count")
    report.plot_hists(st['unaligned_lengths'], st[
                      'aligned_lengths'], l1="Unaligned", l2="Aligned", title="Distribution of read lengths in fractions")
    report.plot_arrays(st['alignment_lengths'], st['aligned_quals'], xlab="Alignment length",
                       ylab="Mean base quality", title="Alignment length vs. mean base quality")


def parse_arguments():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(
        description='Produce plot of quality values across reference.')
    parser.add_argument(
        '-b', metavar='bam_file', type=str, help="BAM file.", required=True)
    parser.add_argument(
        '-r', metavar='report_pdf', type=str, help="Output PDF.", default="aln_qc.pdf")

    args = parser.parse_args()
    return args

args = parse_arguments()
R = Report(args.r)

ist = parse_reads(args.b)
read_qual_qc(ist, R)

# Really regret this, but have to parse the data again:
st = parse_pileups(args.b)
ref_qual_qc(st, R)

R.close()
