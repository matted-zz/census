#!/usr/bin/env python
#
# Census read duplicate histogram generator.
#
# See usage instructions in the README or http://github.com/matted/census/.
#
# Matt Edwards
# Copyright 2014 MIT
# Released under the MIT license
#
#

import pysam, sys, numpy, bisect, re, argparse
from collections import defaultdict

VERSION = "0.9"

def covered(regions, pos):
    targ = bisect.bisect_left(regions, (pos,0))
    if targ == len(regions):
        targ -= 1
    return regions[targ][0] <= pos and regions[targ][1] >= pos or regions[targ-1][0] <= pos and regions[targ-1][1] >= pos

def parse_bed(fin_bed):
    bed = {}
    # fin_bed = open(fname)
    total = 0
    sizes = []
    D = 0

    for line in fin_bed:
        line = line.strip().split()
        chromo = line[0]
        if not chromo in bed:
            bed[chromo] = [(0,0)]
        bed[chromo].append((int(line[1]) - D, int(line[2]) + D))
        total += bed[chromo][-1][1] - bed[chromo][-1][0]
        sizes.append(bed[chromo][-1][1] - bed[chromo][-1][0])

    for chromo, regions in bed.items():
        bed[chromo] = sorted(regions)

    sys.stderr.write("total bases excluded: %d\n" % total)
    sys.stderr.write("total chromosomes included: %d\n" % len(bed))
    
    fin_bed.close()
    return bed

def parse(s, regexp):
    match = re.match(regexp, s)
    if match is None:
        # print >>sys.stderr, "warning: read name regular expression couldn't parse the flowcell position"
        return BAD_COORDS
    else:
        return tuple(map(int, match.groups())) 

def flowcell_dist(x,y):
    if x[0] != y[0] or x[1] != y[1] or x[0] == -1 or y[0] == -1:
        return float("inf")
    return ((x[2]-y[2])**2.0 + (x[3]-y[3])**2.0)**0.5

parser = argparse.ArgumentParser(description="Histogram generator for Census library complexity package.")
parser.add_argument("excluded_regions.bed", type=argparse.FileType('r'))
parser.add_argument("sorted_reads.bam", type=argparse.FileType('rb'))
parser.add_argument("-v", "--version", action="version", version=VERSION)

parser.add_argument("-s", "--single_ended", action="store_true", dest="unpaired", help="Include only single-ended reads, instead of only paired-end reads where both ends map.")
parser.add_argument("-q", "--mapq", type=int, default=1, help="Minimum read mapping quality for a read or read pair to be included.  Default is 1.")
parser.add_argument("-d", "--mindist", type=int, default=100, help="Maximum distance in reported flowcell coordinates for reads to be considered optical duplicates.  Default is 100.")
parser.add_argument("-r", "--regexp", default="[\w\.\-\_ ]+:([\d]):([\d]+):([\d]+):([\d]+).*", help="Regular expression for finding flowcell coordinates from read names.  Default is [\w\.\-\_ ]+:([\d]):([\d]+):([\d]+):([\d]+).*")

args = parser.parse_args()

samf = pysam.Samfile(getattr(args, "sorted_reads.bam").name, 'rb')
bed = parse_bed(getattr(args, "excluded_regions.bed"))

hits = {}
records = 0.01
valid = 0

DEBUG = True
MIN_DIST = args.mindist
PAIRED = not args.unpaired
MIN_MAPQ = args.mapq

maxX = -1
maxY = -1

histo = defaultdict(int)
dupsKilled = 0
badChrKilled = 0
badMapqKilled = 0
cnvsKilled = 0
current = (-1,-1)
BAD_COORDS = (-1,-1,-1,-1)
parseFailure = False

for read in samf:
    if read.is_secondary or read.is_unmapped or read.is_read2 or read.mate_is_unmapped:
        continue
    if read.rname == read.mrnm and read.pos == read.mpos:
        # print "skipping as mirrored read:", read.qname, read.is_paired, read.is_proper_pair, read.mate_is_unmapped, read.isize
        continue

    if records > 0 and records % 1000000 == 0:
        if DEBUG: sys.stderr.write("processed %d reads, max x: %d, max y: %d\n" % (records, maxX, maxY))

    records += 1

    if PAIRED:
        if read.mrnm == -1 or not read.is_read1:
            continue
    else:
        if read.is_read1 or read.is_read2:
            continue
        read.mrnm = read.rname # hack to treat it as the paired case for optical dup detection, etc.
        read.mpos = read.pos # same

    valid += 1

    if read.mapq < MIN_MAPQ:
        badMapqKilled += 1
        continue
       
    refname = samf.getrname(read.rname)
    mateRefname = samf.getrname(read.mrnm) 

    if len(bed) > 0 and (not refname in bed or not mateRefname in bed):
        badChrKilled += 1
        continue

    if len(bed) > 0 and (covered(bed[refname], read.pos) or covered(bed[mateRefname], read.mpos)):
        cnvsKilled += 1
        continue

    if (samf.getrname(read.rname), read.pos) != current:
        for key, names in hits.items():
            if len(names) > 1 and BAD_COORDS not in names: # if we couldn't parse flowcell coordinates, don't try any optical duplicate detection
                legit = set(names)
                if len(legit) > 500:
                    if DEBUG: sys.stderr.write("skipping a cluster (%s) because it has too many hits\n" % key)
                    continue
                for index1, name1 in enumerate(names):
                    for index2 in range(index1+1, len(names)):
                        dist = flowcell_dist(name1, names[index2])
                        if dist < MIN_DIST:
                            legit.discard(names[index2])
                if False and len(legit) >= 4:
                    print("%s %d %s" % (key, len(legit), str(sorted(legit))))
                histo[len(legit)] += 1
                dupsKilled += len(names) - len(legit)
            else:
                histo[len(names)] += 1

        hits.clear()
        current = (samf.getrname(read.rname), read.pos)

    key = (samf.getrname(read.rname), read.pos, samf.getrname(read.mrnm) if read.mrnm != -1 else "", read.mpos, read.isize)

    if not key in hits:
        hits[key] = []
    coords = parse(read.qname, args.regexp)
    if coords == BAD_COORDS:
        parseFailure = True
    hits[key].append(coords)
    if coords[-2] > maxX: maxX = coords[-2]
    if coords[-1] > maxY: maxY = coords[-1]
    
# process hits for the last set (duplicate code...)
for key, names in hits.items():
    if len(names) > 1 and BAD_COORDS not in names:
        legit = set(names)
        if len(legit) > 500:
            if DEBUG: sys.stderr.write("skipping a cluster (%s) because it has too many hits\n" % key)
            continue
        for index1, name1 in enumerate(names):
            for index2 in range(index1+1, len(names)):
                dist = flowcell_dist(name1, names[index2])
                if dist < MIN_DIST:
                    legit.discard(names[index2])
        if False and len(legit) >= 4:
            print("%s %d %s" % (key, len(legit), str(sorted(legit))))
        histo[len(legit)] += 1
        dupsKilled += len(names) - len(legit)
    else:
        histo[len(names)] += 1
      
samf.close()

if parseFailure:
    sys.stderr.write("WARNING: read name regular expression couldn't parse the flowcell position, so no optical duplicate detection was performed\n")

print("# %d considered read pairs, %d valid read pairs, %d optical dup read pairs (%.2f pct.), %d bad chr read pairs (%.2f pct.), %d CNV read pairs (%.2f pct.), %d bad mapq pairs (%.2f pct.)" % (records, valid, dupsKilled, 100.0*dupsKilled/records, badChrKilled, 100.0*badChrKilled/records, cnvsKilled, 100.0*cnvsKilled/records, badMapqKilled, 100.0*badMapqKilled/records))

for hit, count in sorted(histo.items()):
    print("%d\t%d" % (hit, count))
