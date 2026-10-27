#!/usr/bin/env python3
"""
Author : Pierre Levy <levy.pierre@yahoo.fr>
Date   : 2022-05-18
Purpose: Bowtie2 indexing for TCR lib
"""

import argparse
import subprocess
import os
import re
import sys
import logging

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='TCR library screen processing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('lib',
                        metavar='str',
                        help='fasta file containing library tcr sequences')

    return parser.parse_args()

    # --------------------------------------------------

def exec_command(cmd):
    logger = logging.getLogger()
    logger.info(f'{cmd}\n')
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, executable='/bin/bash')
    output, error = p.communicate()
    if p.returncode != 0:
        for line in output.decode("utf-8").split("\n") if output else "":
            logger.error(line.rstrip())
        for line in error.decode("utf-8").split("\n") if error else "":
            logger.error(line.rstrip())
        sys.exit()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    lib = os.path.abspath(args.lib)

    # Set up the logger
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

    # bowtie2 indexing (if not done before)
    if not os.path.exists(f"{os.path.dirname(lib)}/bowtie2_indices/tcrlib.1.bt2"):
	    logger.info("****** Bowtie2: library fasta indexing ******\n\n")
	    cmd = f"mkdir {os.path.dirname(lib)}/bowtie2_indices && \
	            bowtie2-build {lib} {os.path.dirname(lib)}/bowtie2_indices/tcrlib"
	    exec_command(cmd)

# --------------------------------------------------
if __name__ == '__main__':
    main()