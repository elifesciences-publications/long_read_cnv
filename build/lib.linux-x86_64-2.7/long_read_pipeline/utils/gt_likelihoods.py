# Copyright (c) 2017 Boocock James <james.boocock@otago.ac.nz>
# Author: Boocock James <james.boocock@otago.ac.nz>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import math

def _log_choose(n,k):
    """
        Log n choose k 
        
        @author James Boocock   
        @date 8 Feb 2017
    """

    r = 0.0

    if k * 2 > n:
        k = n - k

    for d in xrange(1,k+1):
        r += math.log(n,10)
        r -= math.log(d,10)
        n -= 1

    return r


def genotype_likelihoods(supporting,  not_supporting):
    """
        Get the likelihoods for the reads using counts of alternate and non-alternate bases.

        @author James Boocock
        @date 8 Feb 2017
    """

    prior_s = [0.1,0.9]
    not_supporting = int(not_supporting)
    supporting = int(supporting)
    total_reads = not_supporting + supporting
    log_combo = _log_choose(total_reads, supporting)
    lp_ref = log_combo + supporting * math.log(prior_s[0],10) + not_supporting*math.log(1-prior_s[0], 10)
    lp_alt = log_combo + supporting * math.log(prior_s[1],10) + not_supporting*math.log(1-prior_s[1], 10)
    gt_list = [lp_ref, lp_alt]
    gt_idx = gt_list.index(max(gt_list)) 
    gt_sum = 0
    for gt in gt_list: 
        try:
            gt_sum += 10**gt
        except OverflowError:
            gt_sum += 0 
    if gt_sum > 0:
        gt_sum_log = math.log(gt_sum, 10)
        GQ= abs(-10 * (gt_list[0] - gt_sum_log))
        if GQ > 200:
            GQ = 200
        if lp_alt > lp_ref:
            # Must be alternate
            GT = 1
        else:
            GT = 0
    else:
        if gt_idx == 0:
            GT = 0
        else:
            GT = 1
        GQ = 0 
    return(GQ, GT)
