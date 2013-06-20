#!/usr/bin/python2.6
import pysam, sys, numpy, bisect
from collections import defaultdict

def covered(regions, pos):
    targ = bisect.bisect_left(regions, (pos,0))
    if targ == len(regions):
        targ -= 1
    return regions[targ][0] <= pos and regions[targ][1] >= pos or regions[targ-1][0] <= pos and regions[targ-1][1] >= pos

def parse_bed(fname):
    bed = {}
    fin_bed = open(fname)
    total = 0
    sizes = []
    D = 0

    for line in fin_bed:
        line = line.strip().split()
        chromo = line[0]
        if len(chromo) > 5 or chromo.startswith("chrM") or chromo.startswith("chrX") or chromo.startswith("chrY"):
            continue
        if not bed.has_key(chromo):
            bed[chromo] = [(0,0)]
        bed[chromo].append((int(line[1]) - D, int(line[2]) + D))
        total += bed[chromo][-1][1] - bed[chromo][-1][0]
        sizes.append(bed[chromo][-1][1] - bed[chromo][-1][0])

    for chromo, regions in bed.iteritems():
        bed[chromo] = sorted(regions)

    print >>sys.stderr, "total bases excluded:", total
    print >>sys.stderr, "chromosomes total:", len(bed)
    if len(sizes) > 0:
        print >>sys.stderr, "minimum size block:", min(sizes), numpy.average(sizes), numpy.median(sizes), numpy.std(sizes)
    
    fin_bed.close()
    return bed

# if covered(bed[chromo1], pos1) or covered(bed[chromo2], pos2): continue

def parse(s):
    ret = s.split("/")[0]
    ret = ret.split("#")[0]
    # ret = ret.split(" ")[1] # for SRA reads.
    # ret = s
    # ret = ret.split(":")[1:5]
    # ret = ret.split(":")[3:7] # for DNAse-seq reads, unsure of source
    ret = ret.split(":")[-4:] # for DNAse-seq reads, unsure of source
    # ret[-1] = ret[-1].split("_")[0].split("#")[0] # should regexp it...

    return tuple(map(int, ret))

def flowcell_dist(x,y):
    if x[0] != y[0] or x[1] != y[1]:
        return float("inf")
    return ((x[2]-y[2])**2.0 + (x[3]-y[3])**2.0)**0.5

print >>sys.stderr, "trying to open [%s]" % sys.argv[2]
samf = pysam.Samfile(sys.argv[2], 'rb')
bed = parse_bed(sys.argv[1])

hits = {}
records = 0

DEBUG = True


maxX = -1
maxY = -1

histo = defaultdict(int)
dupsKilled = 0
badChrKilled = 0
badMapqKilled = 0
cnvsKilled = 0
current = (-1,-1)

for read in samf:
    # print dir(read)
    if read.is_secondary or read.is_unmapped or read.is_read2 or read.mate_is_unmapped:
        continue
    if read.rname == read.mrnm and read.pos == read.mpos or read.mrnm == -1:
        # print "skipping as mirrored read:", read.qname, read.is_paired, read.is_proper_pair, read.mate_is_unmapped, read.isize
        continue
    # if records > 15000000:
    #    break

    if records % 1000000 == 0:
        if DEBUG: print >>sys.stderr, "processed %d reads, max x: %d, max y: %d" % (records, maxX, maxY)

    records += 1

    if read.mapq <= 0:
        badMapqKilled += 1
        continue
       
    refname = samf.getrname(read.rname)
    mateRefname = samf.getrname(read.mrnm) 
    if len(bed) > 0 and (not bed.has_key(refname) or not bed.has_key(mateRefname)):
        badChrKilled += 1
        continue

    if len(bed) > 0 and (covered(bed[refname], read.pos) or covered(bed[mateRefname], read.mpos)):
        cnvsKilled += 1
        continue

    if (samf.getrname(read.rname), read.pos) != current:
        ### process hits
        for key, names in hits.iteritems():
            if len(names) > 1:
                # print key, names
                legit = set(names)
                if len(legit) > 500:
                    if DEBUG: print >>sys.stderr, "skipping a cluster because it has too many hits"
                    if DEBUG: print >>sys.stderr, key, names
                    continue
                for index1, name1 in enumerate(names):
                    for index2 in xrange(index1+1, len(names)):
                        dist = flowcell_dist(name1, names[index2])
                        # print "\t", name1, names[index2], dist
                        if dist < 100:
                            legit.discard(names[index2])
                # print "\t", len(legit)
                if False and len(legit) >= 4:
                    print key, len(legit), sorted(legit)
                histo[len(legit)] += 1
                dupsKilled += len(names) - len(legit)
            else:
                histo[len(names)] += 1
        ### end processing
        hits.clear()
        current = (samf.getrname(read.rname), read.pos)

    key = (samf.getrname(read.rname), read.pos, samf.getrname(read.mrnm) if read.mrnm != -1 else "", read.mpos, read.isize)
    if not hits.has_key(key):
        hits[key] = []
    coords = parse(read.qname)
    hits[key].append(coords)
    if coords[-2] > maxX: maxX = coords[-2]
    if coords[-1] > maxY: maxY = coords[-1]
    
### process hits: code duplication!
for key, names in hits.iteritems():
    if len(names) > 1:
        # print key, names
        legit = set(names)
        if len(legit) > 500:
            if DEBUG: print >>sys.stderr, "skipping a cluster because it has too many hits"
            if DEBUG: print >>sys.stderr, key, names
            continue
        for index1, name1 in enumerate(names):
            for index2 in xrange(index1+1, len(names)):
                dist = flowcell_dist(name1, names[index2])
                # print "\t", name1, names[index2], dist
                if dist < 200:
                    legit.discard(names[index2])
        # print "\t", len(legit)
        if False and len(legit) >= 4:
            print key, len(legit), sorted(legit)
        histo[len(legit)] += 1
        dupsKilled += len(names) - len(legit)
    else:
        histo[len(names)] += 1
### end processing
      
samf.close()

# print len(hits)
#histo = defaultdict(int)
#for key in hits.iterkeys():
#    histo[len(hits[key])] += 1

print sys.argv[2]

print "after: ", sorted(histo.iteritems())
print "removed %d (%.2f%%) optical duplicate read pairs" % (dupsKilled, 100.0*dupsKilled/records)
print "removed %d (%.2f%%) bad chr read pairs" % (badChrKilled, 100.0*badChrKilled/records)
print "removed %d (%.2f%%) CNV read pairs" % (cnvsKilled, 100.0*cnvsKilled/records)
print "removed %d (%.2f%%) bad mapq read pairs" % (badMapqKilled, 100.0*badMapqKilled/records)

fout = open("%s.popCNV_filt_histo.txt" % sys.argv[2], 'w')
print >>fout, "#\t%d considered read pairs\t%d valid read pairs\t%d optical dup read pairs\t%.2f\tpct. optical dups\t%d bad chr read pairs\t%.2f\tpct. bad chrs\t%d CNV read pairs\t%.2f\tpct. CNVs" % (records, records-dupsKilled-badChrKilled-cnvsKilled, dupsKilled, 100.0*dupsKilled/records, badChrKilled, 100.0*badChrKilled/records, cnvsKilled, 100.0*cnvsKilled/records)
for hit, count in sorted(histo.iteritems()):
    print >>fout, "%d\t%d" % (hit, count)
fout.close()
