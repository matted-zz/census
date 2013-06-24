#!/usr/bin/env python
#
# Census library complexity estimation routines.
#
# See usage instructions in the README or http://github.com/matted/census/.
#
# Matt Edwards
# Copyright 2013 MIT
# Released under the MIT license
#
#

import sys, numpy, random, math, scipy.optimize, scipy.stats, argparse
from collections import defaultdict

VERSION = "0.8"

def log_trunc(x):
    return numpy.log(max(x, 1e-300))

def K(L, left, right, out=False, cdf=None):
    if out:
        print "K:", scipy.stats.poisson.cdf(right, L), scipy.stats.poisson.cdf(left - 1, L), L, left, right
    if cdf is None:
        return scipy.stats.poisson.cdf(right, L) - scipy.stats.poisson.cdf(left - 1, L)
    else:
        return cdf(right, *L) - cdf(left - 1, *L)

# simplified form that uses the sufficient statistic xavg:
def loglik(L, xavg, left, right):
    # multiply by M (xcount) if we care about the normalization factor afterwards.
    return -(-log_trunc(K(L, left, right)) - L + xavg * log_trunc(L))

# full form that uses the whole histogram:
def loglik2(L, x, left, right, dist=scipy.stats.nbinom):
    LL = 0
    cdf = dist.cdf
    pmf = dist.pmf

    logK = log_trunc(cdf(right, *L) - cdf(left - 1, *L))

    for index in numpy.arange(left, right+1):
        if index >= len(x):
            break
        LL += (log_trunc(pmf(index, *L)) - logK) * x[index]
    return -LL

# target is goal duplication fraction to accept
def find_max_coverage(target, L, CV, RL=100, G=2.634e9):
    def func(cov):
        return (dupfrac(L, CV, cov, RL, G) - target)**2.0
    ret = scipy.optimize.brent(func)
    if func(ret) != func(ret):
        return 1000.0 # float("nan")
    return ret

def estimate(counts, left, right, out=True, verbose=False):
    xsum = sum([count*index if count>=0 and index>=left and index<=right else 0 for index,count in enumerate(counts)])
    xcount = sum([count if count>=0 and index>=left and index<=right else 0 for index,count in enumerate(counts)])
    if xcount <= 0:
        return
    xavg = 1.0 * xsum / xcount

    if out:
        print "Input data summary:"
        print "\tSmallest included event count (inclusive):\t%d" % left
        print "\tLargest included event count (inclusive):\t%d" % right
        print "\tTotal reads:\t%d" % xsum
        print "\tTotal unique reads (molecules):\t%d (%.1f%%)" % (xcount, 100.0*xcount/xsum)
        print "\tAverage reads per molecule:\t", xavg

    # Fit Poisson model.
    ret = scipy.optimize.fminbound(loglik, 0.001, 200.0, args=(xavg, left, right), maxfun=1000, disp=0)

    # Fit negative binomial model.
    # Iterate over several initialization choices to find the best one (or at least a better one):
    bestNegLL = float("inf")
    for init in [0.2, 1.0, 100.0, 1000.0]:
    # for init in [1000.0]:
        ret2 = scipy.optimize.fmin(loglik2, [init, 1.0/(ret/init + 1.0)], args=(counts, left, right), maxiter=1000, disp=False)
        # print ">>", init, loglik2(ret2, counts, left, right), ret2
        if loglik2(ret2, counts, left, right) < bestNegLL:
            bestNegLL = loglik2(ret2, counts, left, right)
            bestRet2 = ret2

    ret2 = bestRet2

    # Fit log-series distribution.
    ret3 = scipy.optimize.fmin(loglik2, [(1.0 - ret2[1]) / (2.0 - ret2[1])], args=(counts, left, right, scipy.stats.logser), maxiter=1000, disp=False)

    alpha = ret2[0]
    beta = (1.0-ret2[1])/ ret2[1]
    if out:
        print "Estimated parameters:"
        print "\tPoisson mean:\t%.4f" % ret
        print "\tNB mean:\t%.4f" % (alpha*beta)
        print "\tNB k:\t%.4f" % (1.0/alpha)
        print "\tNB beta:\t%.4f" % (beta)
        print "\tNB CV:\t%.4f" % (1.0/numpy.sqrt(alpha))

    MLlambda = ret
    libsize = 1.0 * xcount / K(MLlambda, left, right, False)
    libsize2 = 1.0 * xcount / max(1e-12, K(ret2, left, right, False, cdf=scipy.stats.nbinom.cdf))
    libsize3 = 1.0 * xcount / K(ret3, left, right, False, cdf=scipy.stats.logser.cdf)

    fish_alpha = xsum * (1.0 - ret3[0]) / ret3[0]
    # print "LSD alpha estimate (biased because of right truncation):", fish_alpha

    print "Library size estimates:"
    print "\tUnique molecules in 1T reads using Poisson\t %.3f" % (libsize * (1.0 - scipy.stats.poisson.pmf(0, 1e12/libsize)))
    print "\tUnique molecules in 1T reads using NB:\t%.3f" % (libsize2 * (1.0 - scipy.stats.nbinom.pmf(0, ret2[0], 1.0 - (1e12/libsize2)/(ret2[0]+(1e12/libsize2)))))
    print "\tUnique molecules in 1T reads using LSD:\t%.3f" % (fish_alpha * numpy.log(1.0 + 1e12 / fish_alpha))

    # print "maximum useful coverage with NB model (assuming 100bp reads onto mappable hg19):\t%.2f" % find_max_coverage(0.1, libsize2, numpy.sqrt(ret2[0]))
    # if True or ret2[0] >= 100.0: print "1x coverage predicted dupfrac:", dupfrac(libsize2, numpy.sqrt(ret2[0]), 1.0)

    print "\tPoisson predicted library size:\t%.3f" % libsize
    print "\tNB predicted library size:\t%.3f" % libsize2
    print "Applications of library size estimates:"
    print "\tPredicted fraction of library observed, Poisson model:\t%.3f" % K(MLlambda, left, right)
    print "\tPredicted fraction of library observed, NB model:\t%.3f" % (1.0*xcount/libsize2)
    
    print "Maximum useful reads for given experimental efficiency cutoffs:"
    print "\tUniq. Frac.\tReads (Poi.)\tReads (NB)\tReads (LSD)"

    for target_unique_frac in [0.5, 0.8, 0.95]:
        max_reads_poisson = scipy.optimize.brent(lambda reads : (libsize * (1.0 - scipy.stats.poisson.pmf(0, max(1, reads) / libsize)) / max(1, reads) - target_unique_frac)**2.0, brack=(1e5, 1e9))
        max_reads_NB = scipy.optimize.brent(lambda reads : (libsize2 * (1.0 - scipy.stats.nbinom.pmf(0, ret2[0], 1.0 - beta*max(1, reads)/xsum/(1.0+beta*max(1, reads)/xsum))) / max(1, reads) - target_unique_frac)**2.0, brack=(1e5, 1e9))
        max_reads_LSD = scipy.optimize.brent(lambda reads : (fish_alpha * numpy.log(1.0 + max(1, reads) / fish_alpha) / max(1, reads) - target_unique_frac)**2.0, brack=(1e5, 1e9))
        
        print "\t%.2f:\t\t%.1f\t%.1f\t%.1f" % (target_unique_frac, max_reads_poisson, max_reads_NB, max_reads_LSD)

    if verbose:
        if out:
            print "Observed and predicted count table:"
            print "\tHits\tObs.\tPoi.\tNB\tLSD"
        printed = False
        total = 0.0
        cutoff = None
        for index in numpy.arange(0, min(len(counts), right+1)):
            count = counts[index]
            temp = libsize * scipy.stats.poisson.pmf(index, MLlambda)
            temp2 = scipy.stats.poisson.pmf(index, MLlambda)
            temp3 = scipy.stats.nbinom.pmf(index, *ret2)
            temp4 = libsize2 * scipy.stats.nbinom.pmf(index, *ret2)
            temp5 = libsize3 * scipy.stats.logser.pmf(index, *ret3)
            if total >= 0.99 and not printed:
                printed = True
                cutoff = index
            total += temp2
            if count >= 0 and (True or count >= 10 or index==0):
                if out:
                    print "\t%d:\t%d\t%.0f\t%.0f\t%.1f" % (index, count, temp, temp4, temp5)

    ### 
    print "Prediction of observed unique molecules with given read counts:"
    print "\tReads\tUnique (Poi.)\tUnique (NB)\tUnique (LSD)"
    for more_reads in [xsum] + range(int(50e6), int(500e6+1), int(100e6)) + [1e9] + [2e9]:
        L = 1.0*(more_reads)/libsize
        new_beta = beta * L / (xsum / libsize) # scale beta to stay proportional to lambda

        # poisson estimate
        temp1 = libsize * (1.0 - scipy.stats.poisson.pmf(0, L))
        # NB estimate
        temp2 = libsize2 * (1.0 - scipy.stats.nbinom.pmf(0, ret2[0], 1.0 - new_beta/(1.0+new_beta)))
        # LSD estimate
        temp3 = fish_alpha * numpy.log(1.0 + more_reads / fish_alpha)

        print "\t%dM:\t%.1f\t%.1f\t%.1f" % (more_reads/1000000, temp1, temp2, temp3)
    ###

    return libsize, xcount

def map_to_histo_vec(temp):
    counts = numpy.zeros((max([-1]+temp.keys())+1))
    for k,v in temp.iteritems():
        counts[k] = v
    return counts

# alternate optimization procedure (unused for now):
def estimate_EM(counts, left, right):
    xsum = sum([count*index if count>=0 and index>=left and index<=right else 0 for index,count in enumerate(counts)])
    xcount = sum([count if count>=0 and index>=left and index<=right else 0 for index,count in enumerate(counts)])
    if xcount <= 0:
        return
    xavg = xsum / xcount

    if False:
        print "Total reads:\t%d" % xsum
        print "Total unique reads (molecules):\t%d" % xcount
        print "Average reads per molecule:\t", xavg

    L = xavg
    for its in xrange(20):
        # unobs = numpy.exp(-L)/(1.0-numpy.exp(-L)) * xcount
        unobs = (1.0 / K(L, left, right, False) - 1.0) * xcount
        L = 1.0 * xcount / (xcount + unobs) * xavg # need to account for right truncation here...
        print "\t", L, unobs
    unobs = numpy.exp(-L)/(1.0-numpy.exp(-L)) * xcount
    print "EM results:\tlambda: %.3f, library size: %.2f" % (L, unobs+xcount)

    return unobs+xcount

def sample_histo_fast(histo, frac=0.5):
    choices = []
    for index, count in enumerate(histo):
        choices.extend(range(len(choices),len(choices)+int(count))*index)
    numpy.random.shuffle(choices)
    sampled = choices[:scipy.stats.binom.rvs(len(choices), frac)]
    hits = {}
    for sample in sampled:
        if not hits.has_key(sample):
            hits[sample] = 0
        hits[sample] += 1
    histo = {}
    for hit, count in hits.iteritems():
        if not histo.has_key(count):
            histo[count] = 0
        histo[count] += 1
    return map_to_histo_vec(histo)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Census, library complexity estimation.")
    parser.add_argument("count_histogram.txt", type=argparse.FileType('r'), help="File for duplicate count histogram, or - for stdin.")
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    parser.add_argument("-l", "--mincount", type=int, default=1, help="Minimum duplicate count to use in estimation.  Default is 1.")
    parser.add_argument("-r", "--maxcount", type=int, default=10, help="Maximum duplicate count to use in estimation.  Default is 10.")
    parser.add_argument("-s", "--subsample", type=float, default=1, help="Fraction of counts to use (float), useful for testing.  Default is 1 (no downsampling).")
    
    args = parser.parse_args()

    temp = {}
    f = getattr(args, "count_histogram.txt")
    left = args.mincount
    right = args.maxcount

    for line in f:
        if line[0] == "#": continue
        index, count = map(int, line.strip().split()[:2])
        temp[index] = count
        counts = map_to_histo_vec(temp)

    if counts[0] != 0:
        print "true 0 counts, which we are ignoring in the estimation procedure:", counts[0]
    counts[0] = 0

    if args.subsample < 1.0 and args.subsample > 0.0:
        subsample = sample_histo_fast(counts, args.subsample)
        print >>sys.stderr, "randomly downsampling to %.4f of reads" % args.subsample
    else:
        subsample = counts

    libsize, unique_old = estimate(subsample, left, right, out=True, verbose=True)
