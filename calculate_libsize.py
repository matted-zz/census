#!/usr/bin/env python
import sys, numpy, random, math
import scipy.optimize, scipy.stats
from collections import defaultdict

def log_trunc(x):
    return numpy.log(max(x, 1e-300))

def K(L, left, right, out=False, cdf=None):
    if out:
        print "K:", scipy.stats.poisson.cdf(right, L), scipy.stats.poisson.cdf(left - 1, L), L, left, right
    if cdf is None:
        return scipy.stats.poisson.cdf(right, L) - scipy.stats.poisson.cdf(left - 1, L)
    else:
        return cdf(right, *L) - cdf(left - 1, *L)

# simplified form that uses the sufficient statistic xavg
def loglik(L, xavg, left, right):
    # did I drop some extra terms, required for normalization? YES, multiply by M (xcount)
    return -(-log_trunc(K(L, left, right)) - L + xavg * log_trunc(L))

# fitting new distribution, propto p^x / x.  needs suff1 = sum_i (i*x_i) and suff2 = sum_i logi * x_i
def loglik3(p, suff1, suff2, n):
    # print suff1, suff2, n
    # + n*numpy.log(p/(1.0-p))
    return -(numpy.log(p)*suff1 - suff2 - n*numpy.log(numpy.log(1.0/(1.0-p))))
    
def loglik_with_size(params, xavg, left, right, observed):
    L = params[0]
    k = K(L, left, right)
    M = observed / k
    # explanation: we previously dropped a factor of "n" because it
    # didn't affect the maximization.  now it does, so we add it back.
    # we also include the likelihood term for the library size
    # explicitly.  we play some tricks: we want max_{L,M}
    # loglik(L,M,x), which we can break into max_L f(L,x) * max_M
    # g(M,x).  we know that the best M for a given L is observed/k and
    # the log likelihood of that portion is 1/sqrt(2*pi*M*k*(1-k))
    # (from normal approximation to binomial). that gives us the new
    # term in the objective function.
    # print observed*(-numpy.log(k) - L + xavg * numpy.log(L)), -0.5*(numpy.log(M) + numpy.log(k) + numpy.log(1.0-k))

    ret = -(observed*(-numpy.log(k) - L + xavg * numpy.log(L)) - 0.5*(numpy.log(M) + numpy.log(k) + numpy.log(1.0-k)))
    return ret

# full form that uses the whole histogram
def loglik2(L, x, left, right, dist=scipy.stats.nbinom):
    LL = 0
    cdf = dist.cdf
    pmf = dist.pmf
    
    # if len(L) > 1: L[0] = max(L[0], 1.001)
    # if K(L, left, right, False, cdf=scipy.stats.nbinom.cdf) < 0.001:
    #    return numpy.inf

    logK = log_trunc(cdf(right, *L) - cdf(left - 1, *L))

    for index in numpy.arange(left, right+1):
        if index >= len(x):
            break
        LL += (log_trunc(pmf(index, *L)) - logK) * x[index]
    return -LL

def loglik2prior(L, x, left, right, dist=scipy.stats.nbinom, priorB=1.0):
    A = 5.0
    B = priorB
    # LL = numpy.log(scipy.stats.gengamma.pdf(1.0/L[0], A, 1.0, scale=1.0/B))
    LL = (A-1.0)*log_trunc(1.0/L[0]) - 1.0/L[0]*B - log_trunc(scipy.special.gamma(A)) + A*log_trunc(B)

    cdf = dist.cdf
    pmf = dist.pmf
    
    # if len(L) > 1: L[0] = max(L[0], 1.001)
    # if K(L, left, right, False, cdf=scipy.stats.nbinom.cdf) < 0.001:
    #    return numpy.inf

    logK = log_trunc(max(cdf(right, *L) - cdf(left - 1, *L), 1e-100))

    for index in numpy.arange(left, right+1):
        if index >= len(x):
            break
        LL += (log_trunc(max(pmf(index, *L), 1e-300)) - logK) * x[index]
    return -LL

def dupfrac(L, CV, C, RL=100, G=2.634e9):
    k = CV**2.0
    R = C*G/RL
    Rprime = L*(1.0 - scipy.stats.nbinom.pmf(0, 1.0/k, 1.0 - R/L/(1.0/k+R/L))) # number of UNIQUE molecules in the read sample.  should use this for coverage calculation, not the number of reads.

    return 1.0 - Rprime / R

    # return (1.0 - scipy.stats.nbinom.pmf(0, 1.0/k, 1.0 - R/L/(1.0/k+R/L)) - scipy.stats.nbinom.pmf(1, 1.0/k, 1.0 - R/L/(1.0/k+R/L))) / (1.0 - scipy.stats.nbinom.pmf(0, 1.0/k, 1.0 - R/L/(1.0/k+R/L)))
    # print "distr:", scipy.stats.nbinom.pmf(0, 1.0/k, 1.0 - R/L/(1.0/k+R/L)), scipy.stats.nbinom.pmf(1, 1.0/k, 1.0 - R/L/(1.0/k+R/L)), scipy.stats.nbinom.pmf(2, 1.0/k, 1.0 - R/L/(1.0/k+R/L)), 1.0 - scipy.stats.nbinom.pmf(0, 1.0/k, 1.0 - R/L/(1.0/k+R/L)) - scipy.stats.nbinom.pmf(1, 1.0/k, 1.0 - R/L/(1.0/k+R/L)) - scipy.stats.nbinom.pmf(2, 1.0/k, 1.0 - R/L/(1.0/k+R/L))
    # return 1.0 - scipy.stats.nbinom.pmf(1, 1.0/k, 1.0 - R/L/(1.0/k+R/L)) / (1.0 - scipy.stats.nbinom.pmf(0, 1.0/k, 1.0 - R/L/(1.0/k+R/L)))

# target is goal duplication fraction to accept
def find_max_coverage(target, L, CV, RL=100, G=2.634e9):
    def func(cov):
        return (dupfrac(L, CV, cov, RL, G) - target)**2.0
    ret = scipy.optimize.brent(func)
    if func(ret) != func(ret):
        return 1000.0 # float("nan")
    return ret

def estimate(counts, left, right, out=True, verbose=False, prior=1.0):
    xsum = sum([count*index if count>=0 and index>=left and index<=right else 0 for index,count in enumerate(counts)])
    xcount = sum([count if count>=0 and index>=left and index<=right else 0 for index,count in enumerate(counts)])
    if xcount <= 0:
        return
    xavg = xsum / xcount

    if out:
        print "Smallest included event count (inclusive):\t%d" % left
        print "Largest included event count (inclusive):\t%d" % right
        print "Total reads:\t%d" % xsum
        print "Total unique reads (molecules):\t%d (%.1f%%)" % (xcount, 100.0*xcount/xsum)
        print "Average reads per molecule:\t", xavg
        # print "Count variance:\t%.4f, \tMax variance:\t%.4f" % (sum([count*(index-xavg)**2.0 if count>=0 and index>=left and index<=right else 0 for index,count in enumerate(counts)]) / xcount, xavg+xavg**2.0)
    ret = scipy.optimize.fminbound(loglik, 0.01, 20.0, args=(xavg, left, right), maxfun=1000, disp=0)
    # ret2 = scipy.optimize.fmin(loglik2, [0.05, 0.9], args=(counts, left, right), maxiter=1000, disp=False)
    # original initialization:
    # ret2 = scipy.optimize.fmin(loglik2, [0.2, 0.99], args=(counts, left, right), maxiter=1000, disp=False)
    # better initialization (try to set overdisp=0, mean=poisson mean)
    ret2 = scipy.optimize.fmin(loglik2, [1000.0, 1.0/(ret/1000.0 + 1.0)], args=(counts, left, right), maxiter=1000, disp=False)
    # TODO: do multiple and compare?
    # ret2 = scipy.optimize.fmin(loglik2prior, [1.5, 0.8], args=(counts, left, right, scipy.stats.nbinom, prior), maxiter=1000, disp=False)
    # print "using prior!"
    ret3 = scipy.optimize.fmin(loglik2, [(1.0 - ret2[1]) / (2.0 - ret2[1])], args=(counts, left, right, scipy.stats.logser), maxiter=1000, disp=False)

    # ret3 = scipy.optimize.fmin(loglik3, [0.1], args=(xsum, sum([numpy.log(index)*count if index>=left and index<=right else 0 for index,count in enumerate(counts)]), xcount), maxiter=1000)
    # ret3 = scipy.optimize.fmin(loglik2, [5.0, 0.5], args=(counts, left, right), maxiter=1000)
    # ret3 = scipy.optimize.fmin_tnc(loglik2, [1.5, 0.8], args=(counts, left, right), approx_grad=True)
    
    # ret3 = scipy.optimize.fmin(loglik_with_size, [xavg], args=(xavg, left, right, xcount), maxiter=1000)
    # print "new method results:", ret3

    if False: # new for debugging/thinking
        import pylab
        L = numpy.arange(0.001, 4, 0.01)
        LL = [-loglik(l, xavg, left, right) for l in L]
        # LL2 = [-loglik2([1000.0,1.0 - l/(1000.0+l)], counts, left, right) for l in L]
        LL2 = [-loglik2([l, 0.784], counts, left, right) for l in L]
        pylab.plot(L, LL)
        pylab.show()

    alpha = ret2[0]
    beta = (1.0-ret2[1])/ ret2[1]
    if out:
        # print "ret:", ret, "LL:", loglik2(ret2, counts, left, right, scipy.stats.poisson)
        # print "ret2:", ret2, "LL:", loglik2prior(ret2, counts, left, right), 
        print "Best-fit mean: %.2f, stdev: %.2f," % (alpha*beta, numpy.sqrt(alpha)*beta), "lambda: %.4f, k: %.4f, CV: %.4f" % (alpha*beta, 1.0/alpha, 1.0/numpy.sqrt(alpha))
        # print "ret3:", ret3,  "LL:", loglik2(ret2, counts, left, right)

        fisher = scipy.misc.derivative(loglik, ret, dx=0.01, n=2, args=(xavg, left, right), order=7)
        error = scipy.stats.norm.isf((1.0-0.5)/2.0, 0, numpy.sqrt(1.0/fisher))
        # print "\td^2 L / dx^2:", fisher
        # print "\t50% confidence interval error bar using observed fisher information:", error
        # print "\tlibsize from %.1f to %.1f" % (xcount / K(max(0.01, ret-error), left, right, False), xcount / K(ret+error, left, right, False))
        # print "d^2 L / dr^2:", scipy.misc.derivative(lambda x,a,b,c: loglik2([x,ret2[1]], a,b,c), ret2[0], dx=0.01, n=2, args=(counts, left, right), order=7)
        # print "d^2 L / dp^2:", scipy.misc.derivative(lambda x,a,b,c: loglik2([ret2[0],x], a,b,c), ret2[1], dx=0.01, n=2, args=(counts, left, right), order=7)


    MLlambda = ret # API changing?
    libsize = xcount / K(MLlambda, left, right, False)
    libsize2 = xcount / K(ret2, left, right, False, cdf=scipy.stats.nbinom.cdf)
    libsize3 = xcount / K(ret3, left, right, False, cdf=scipy.stats.logser.cdf)

    fish_alpha = xsum * (1.0 - ret3[0]) / ret3[0]
    print "LSD alpha estimate (a bit wrong because of right truncation):", fish_alpha
    print "Unique molecules in 1T reads using LSD:", fish_alpha * numpy.log(1.0 + 1e12 / fish_alpha)
    new_beta = beta * 1e12 / xsum
    print "Unique molecules in 1T reads using NB:", libsize2 * (1.0 - scipy.stats.nbinom.pmf(0, ret2[0], 1.0 - (1e12/libsize)/(ret2[0]+(1e12/libsize)))) #, (1.0 - scipy.stats.nbinom.pmf(0, ret2[0], ret2[1]))
    print "Unique molecules in 1T reads using poisson:", libsize * (1.0 - scipy.stats.poisson.pmf(0, 1e12/libsize))

    # print "maximum useful coverage with NB model (assuming 100bp reads onto mappable hg19):\t%.2f" % find_max_coverage(0.1, libsize2, numpy.sqrt(ret2[0]))
    # if True or ret2[0] >= 100.0: print "1x coverage predicted dupfrac:", dupfrac(libsize2, numpy.sqrt(ret2[0]), 1.0)

    # what about truncation???
    pr0_po = scipy.stats.poisson.pmf(0, ret)
    pr0_nb = scipy.stats.nbinom.pmf(0, *ret2)

    # print "alternate estimates and variances:"
    temp = scipy.stats.nbinom.stats(xcount, 1.0 - pr0_po)
    # print xcount+temp[0], numpy.sqrt(temp[1])
    temp = scipy.stats.nbinom.stats(xcount, 1.0 - pr0_nb)
    # print xcount+temp[0], numpy.sqrt(temp[1])
    # print "K_NB:", 1.0*xcount/libsize2

    if out:
        print "ML lambda:\t%.3f\nPoisson observed molecules fraction:\t%.3f\nNB observed molecules fraction:\t%.3f\nPoisson ML library size:\t%.1f\nNB ML library size:\t%.1f" % (MLlambda, K(MLlambda, left, right), 1.0*xcount/libsize2, libsize, libsize2)

    if verbose:
        expected, observed, rangex, expected2 = [], [], [], []
        if out:
            print "Hits\tObs.\tPoi.\tNB\tLSD" #\tFrac1\tFrac2"
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
                # print "cutoff:", index
                cutoff = index
            total += temp2
            if count >= 0 and (True or count >= 10 or index==0):
                # print 1.0*count/libsize, 1.0*count/libsize2
                if out:
                    # print "%d:\t%d\t%.0f\t%.0f\t%.3f\t%.3f\t%f\t%f" % (index, count, temp, temp4, temp2, temp3, scipy.stats.binom.sf(count, libsize, 1.0*temp/libsize), scipy.stats.binom.sf(count, libsize2, 1.0*temp4/libsize2))
                    print "%d:\t%d\t%.0f\t%.0f\t%.1f" % (index, count, temp, temp4, temp5)
                expected.append(temp)
                observed.append(count)
                rangex.append(index)
                expected2.append(temp4)
        expected = numpy.array(expected)
        expected2 = numpy.array(expected2)
        observed = numpy.array(observed)
        start, finish = 1, 21

        # print >>sys.stderr, "vectors we're testing:", observed[start:finish], expected[start:finish], expected2[start:finish]
        if out and False:
            print "chisq obs against expected (high pval good):", scipy.stats.chisquare(observed[start:finish], expected[start:finish]), scipy.stats.chisquare(observed[start:finish], expected2[start:finish])
            print >>sys.stderr, "K-S 2-sample test, obs against expected (high pval good):", scipy.stats.ks_2samp(observed[start:finish], expected[start:finish]), scipy.stats.ks_2samp(observed[start:finish], expected2[start:finish])
            print >>sys.stderr, "Count cutoff at FDR of 0.01 is:\t", cutoff
            print >>sys.stderr, "Expected number of unique reads (molecules):\t%.2f\t%.2f" % (libsize - expected[0], libsize2 - expected2[0])
            print >>sys.stderr, "Expected number of new unique reads (molecules) with..."
        new_unique = None
        for more_reads in [xsum, xsum*2] + range(1000000, 10000001, 1000000):
            break
            L = 1.0*(xsum+more_reads)/libsize
            sig = numpy.sqrt(alpha*beta**2.0)
            new_beta = beta * L / (xsum / libsize) # scale beta to stay proportional to lambda
            # poisson estimate
            temp1 = -libsize * scipy.stats.poisson.pmf(0, L) + libsize * scipy.stats.poisson.pmf(0, MLlambda)
            # correct (??) NB estimate
            temp2 = -libsize2 * scipy.stats.nbinom.pmf(0, ret2[0], 1.0 - new_beta/(1.0+new_beta)) + libsize2 * scipy.stats.nbinom.pmf(0, ret2[0], ret2[1])
            print "\t", scipy.stats.nbinom.stats(ret2[0], 1.0 - new_beta/(1.0+new_beta)), new_beta, L, L*new_beta+L
            # wrong (old) NB estimate
            # temp3 = -libsize2 * scipy.stats.nbinom.pmf(0, ?, ?) + libsize2 * scipy.stats.nbinom.pmf(0, ret2[0], ret2[1])
            if more_reads == xsum:
                new_unique = (temp1, temp2)
            if out:
                print "\t%dM more reads:\t%.2f\t%.2f" % (more_reads/1000000, temp1, temp2)
        if out and False:
            print >>sys.stderr, "Sum of expected:", sum(expected)
            print 1.0*(xsum+xsum)/libsize, 2.0*MLlambda, K(MLlambda, left, right, False)

    return libsize, xcount, xavg, new_unique, ret2

def map_to_histo_vec(temp):
    counts = numpy.zeros((max([-1]+temp.keys())+1))
    for k,v in temp.iteritems():
        counts[k] = v
    return counts

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
    """
    hits = defaultdict(int)
    for index, count in enumerate(histo):
        for c in xrange(count):
            hits[scipy.stats.binom.rvs(index, frac)] += 1
    return map_to_histo_vec(hits)
    """

# doesn't do anything with items that were observed 0 times
def sample_histo(histo, frac=0.5):
    choices = []
    for index, count in enumerate(histo):
        choices.extend(range(len(choices),len(choices)+int(count))*index)
    sampled = []
    for choice in choices:
        if random.random() < frac:
            sampled.append(choice)
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
    if len(sys.argv) in [2,3,4,5]:
        # get counts from command line
        temp = {}
        rev = ("REV" in sys.argv[1:])
        if sys.argv[1] == "-":
            f = sys.stdin
        else:
            f = open(sys.argv[1])
            print sys.argv[1]
        if len(sys.argv) >= 4:
            left, right = map(int, sys.argv[2:4])
        else:
            left, right = 1.0, 10.0
        for line in f:
            # index, count = map(int, line.strip().split())
            if line[0] == "#": continue
            index, count = map(int, line.strip().split()[:2])
            if rev:
                index, count = count, index
            temp[index] = count
        counts = map_to_histo_vec(temp)
    elif len(sys.argv) > 5:
        rev = False
        left, right = map(int, sys.argv[1:3])
        for fname in sys.argv[3:]:
            f = open(fname)
            temp = {}
            for line in f:
                if line[0] == "#": continue
                index, count = map(int, line.strip().split()[:2])
                if rev:
                    index, count = count, index
                temp[index] = count
            counts = map_to_histo_vec(temp)
            print "*** processing %s ***" % fname
            for rep in xrange(1):
                #print "** taking %d/10 reads **" % (rep+1)
                #subsample = sample_histo_fast(counts, (rep+1.0)/10.0)
                subsample = counts
                estimate(subsample, left, right, out=True, verbose=True)
            print "*** done ***\n"
            # sys.exit(0)
        sys.exit(0)
    else:
        print >>sys.stderr, "usage: python %s counts.txt [min_count max_count] [REV]" % sys.argv[0]
        sys.exit(1)

    if True:
        if counts[0] != 0:
            print "true 0 counts:", counts[0]
        counts[0] = 0
        orig1counts = float(counts[1])
        for rep in xrange(1):
            for innerrep in xrange(1):
                # print "rep", (innerrep+1)
                # subsample = sample_histo_fast(counts, (rep+1.0)/10.0)
                # subsample = sample_histo_fast(counts, 0.05)
                # print >>sys.stderr, "downsampling to 5% of reads!"
                subsample = counts
                libsize, unique_old, ignored, predicted_unique, ignored = estimate(subsample, left, right, out=True, verbose=True, prior=1.0)

            """
            libsize, unique_new, ignored = estimate(counts, left, right, out=True, verbose=True)
            print "new unique reads: %d, prediction was:" % (unique_new - unique_old), predicted_unique
            for reps in xrange(10):
                subsample = sample_histo(counts, (rep+1.0)/10.0)
                libsize, unique_old, predicted_unique = estimate(subsample, left, right, out=True, verbose=True)

                # def loglik2(L, x, left, right, dist=scipy.stats.nbinom):
                # ret2 = scipy.optimize.fmin(loglik2, [1.5, 0.8], args=(counts, left, right), maxiter=1000)
                # i want to plot the likelihood function for values of k for fixed (true) lambda.
            """
    else:
        print "true 0 counts:", counts[0]
        counts[0] = 0
        # counts = sample_histo(counts, 0.1)
        ignored, xcount, xavg, ignored, params = estimate(counts, left, right, out=True, verbose=True)
        # estimate_EM(counts, left, right)
        
        import pylab
        pylab.subplot(311)
        krange = numpy.arange(0.01,2.0,0.1)
        L = 0.8
        pylab.plot(krange, [loglik2([1.0/k, 1.0 - L/(1.0/k+L)], counts, left, right) for k in krange])
        L = params[0]*(1.0 - params[1]) / params[1]
        pylab.plot(krange, [loglik2([1.0/k, 1.0 - L/(1.0/k+L)], counts, left, right) for k in krange])

        Lrange = numpy.arange(0.02, 3.0, 0.025)
        pylab.subplot(312)
        for k in krange:
            # pylab.plot(krange, [loglik2([1.0/k, 1.0 - L/(1.0/k+L)], counts, left, right) for L in Lrange])
            def target(L):
                return loglik2([1.0/k, 1.0 - L/(1.0/k+L)], counts, left, right)
            Lhat = scipy.optimize.fminbound(target, min(Lrange), max(Lrange))
            print "\t", k, Lhat, xcount / K([1.0/k, 1.0 - Lhat/(1.0/k+Lhat)], left, right, False, cdf=scipy.stats.nbinom.cdf), -loglik2([1.0/k, 1.0 - Lhat/(1.0/k+Lhat)], counts, left, right), (1.0 - scipy.stats.nbinom.pmf(0, 1.0/k, 1.0 - Lhat/(1.0/k+Lhat))), xcount*(1.0 - scipy.stats.nbinom.cdf(4, 1.0/k, 1.0 - Lhat/(1.0/k+Lhat))), 1.0 - scipy.stats.gengamma.cdf(1e-6, 1.0/k, 1.0, scale=Lhat*k)

            #print "\t", scipy.integrate.quad(lambda x: numpy.exp(-loglik(x, xavg, left, right)), Lhat-0.1, Lhat+0.1)
            # print "\t", scipy.integrate.quad(lambda x: numpy.exp(-loglik(x, xavg, left, right)), 0, 100)
            # print "\t", scipy.integrate.dblquad(lambda k,L: numpy.exp(-loglik2([1.0/k, 1.0 - L/(1.0/k+L)], counts, left, right)), k-0.01, k+0.01, lambda x: L-0.01, lambda x: L+0.01)

            # i need a prior on L that integrates to 1...

        pylab.subplot(313)

        # idiot, P(data|params) summed over params is not a probability distribution... need to read more gelman
        # to figure out reasonable things to do.  MCMC??? ridiculous.  or closed-form intelligence with EM.
        Lrange = numpy.arange(0.0005, 10.0, 0.01)
        xcount = 1.575
        print "\t", xcount
        pylab.plot(Lrange, [numpy.exp(-loglik(L, xavg, left, right)*xcount) for L in Lrange])

        # xavg and xcount should be from the complete data here, not the truncated data...
        pylab.plot(Lrange, [scipy.stats.gengamma.pdf(L, xavg*xcount, 1.0, scale=1.0/xcount) for L in Lrange], ":")

        pylab.ylabel("loglik")
        print "\t", scipy.integrate.quad(lambda x: numpy.exp(-loglik(x, xavg, left, right)*xcount), 0, Lhat)
        print "\t", scipy.integrate.quad(lambda x: numpy.exp(-loglik(x, xavg, left, right)*xcount), Lhat, 30)
        print "\t", scipy.integrate.quad(lambda x: numpy.exp(-loglik(x, xavg, left, right)*xcount), Lhat, Lhat+1.25)
        
        print "mean lambda:", sum([L*numpy.exp(-loglik(L, xavg, left, right)*xcount) for L in Lrange]) / sum([numpy.exp(-loglik(L, xavg, left, right)*xcount) for L in Lrange])

        """
        pylab.clf()
        negloglik = numpy.array([[loglik2([1.0/k, 1.0 - L/(1.0/k+L)], counts, left, right) for L in Lrange] for k in krange])
        pylab.imshow(-negloglik, interpolation="nearest", origin="lower", extent=(0,max(Lrange),0,max(krange)))
        pylab.xlabel("$\lambda$")
        pylab.ylabel("k")
        pylab.colorbar()
        """
        pylab.show()
