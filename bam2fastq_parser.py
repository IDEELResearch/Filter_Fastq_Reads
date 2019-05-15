#!/usr/bin/python2

"""
This entire script largely relies on
dpryan79's answer to this stackexchange question
https://bioinformatics.stackexchange.com/questions/2149/how-to-safely-and-efficiently-convert-subset-of-bam-to-fastq
https://github.com/dpryan79/Answers/blob/master/bioinfoSE_2149/convert.py
with very minor modifications
"""

import pysam
import argparse
import sys

def parseargs():
    parser = argparse.ArgumentParser(description = "Extract a bam file mapping quality and exclude reads with a MQ above a specified threshold.")
    parser.add_argument("-b", "--BAMfile",
                        help = "BAM file to extract/parse")
    parser.add_argument("-m", "--MQthreshold",
                        type = int,
                        help = "Mapping Quality Threshold to Exclude On. Reads with a MQ of this or greater than are excluded.")
    parser.add_argument("-o", "--outputPrefix",
                        help = "The prefix for the newly created fastq files. They will be name prefix_1.fastq and prefix_2.fastq")
    return parser

def filterAlignment(bam, mq):
    """
    Only include those reads that fall below our speficied map quality.
    Note, we are not explicitly doing something about multimappers,
     which may actually align well
    to the first reference genome (here) but just are hypervariable.
    """
    if bam.mapping_quality > mq:
        return True
    else:
        return False


def inBuffer(b, buf):
    if b.query_name in buf:
        return True
    return False


def revComp(s):
    d = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    s = [d[c] for c in s]
    return ''.join(s[::-1])


def storeInBuffer(b, buf):
    s = b.query_sequence
    q = ''.join([chr(c+33) for c in b.query_qualities])
    if b.is_read1:
        rNum = 1
    else:
        rNum = 2
    if b.is_reverse:
        s = revComp(s)
        q = q[::-1]
    buf[b.query_name] = (s, q, rNum)


def dumpAlignments(b, buf, of1, of2):
    b2 = buf[b.query_name]

    # Write the mate
    of = of1
    if b2[2] == 2:
        of = of2
    of.write("@{}\n{}\n+\n{}\n".format(b.query_name, b2[0], b2[1]))

    # Write the read
    s = b.query_sequence
    q = ''.join([chr(c+33) for c in b.query_qualities])
    if b.is_reverse:
        s = revComp(s)
        q = q[::-1]
    of = of1
    if b.is_read2:
        of = of2
    of.write("@{}\n{}\n+\n{}\n".format(b.query_name, s, q))

    # Remove from the buffer
    del buf[b.query_name]


def main(args):
    args = parseargs().parse_args(args)

    bam = pysam.AlignmentFile(args.BAMfile, "rb")
    mqnum = args.MQthreshold
    of1 = open("{}_1.fastq".format(args.outputPrefix), "w")
    of2 = open("{}_2.fastq".format(args.outputPrefix), "w")

    buf = dict()
    for b in bam:
        if filterAlignment(b, mqnum):
            continue
        if inBuffer(b, buf):
            dumpAlignments(b, buf, of1, of2)
        else:
            storeInBuffer(b, buf)

    bam.close()
    of1.close()
    of2.close()


if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]
    main(args)
