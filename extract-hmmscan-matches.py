#! /usr/bin/env python

import sys
import argparse
import screed


def tblout_to_names(tblout):
    names = []
    with open(tblout, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                line = line.strip().split()
                names.append(line[2])
    return names


def get_fasta_matches(tblout, fasta):
    names = set(tblout_to_names(tblout))
    found = set()

    for record in screed.open(args.fasta):
        shortname = record.name.split()[0]
        if shortname in names:
            sys.stdout.write(f'>{record.name}\n{record.sequence}\n')
            found.add(shortname)

    sys.stderr.write(f'found {len(found)} of {len(names)} ({len(names-found)} missing)\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('tblout')
    p.add_argument('--fasta')
    args = p.parse_args()
    #assert args.output, "must specify location for output configfile using '-o'"
    sys.exit(get_fasta_matches(args.tblout, args.fasta))
