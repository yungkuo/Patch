"""
Function for burst search in cython.
"""

import numpy as np
cimport numpy as np


def burstsearch_c(np.float64_t[:] signal,
                  np.int32_t m,
                  np.float64_t threshold):
    """Sliding window burst search. Cython version.

    Returns:
        Record array of burst data, one element per burst.
        Each element is a composite data type containing
        burst start, burst stop and burst score.
    """
    cdef char in_burst
    cdef np.int32_t i, ii, start, stop
    cdef np.float64_t score, totalscore

    bursts = []
    in_burst = 0

    score = 0.
    for i in range(m):
        score += signal[i]

    for i in range(signal.size - m + 1):
        if i > 0:
            score += signal[i + m - 1] - signal[i - 1]
        if score > threshold:
            if not in_burst:
                # index of first cycle in burst
                start = i
                in_burst = 1
                totalscore = 0
                for ii in range(m-1):
                    totalscore += signal[i + ii]
            totalscore += signal[i + m - 1]
        elif in_burst:
            # index of last cycle in burst
            stop = i + m - 2
            bursts.append((start, stop, totalscore))
            in_burst = 0

    # Create numpy recarray
    dt = np.dtype([('start','int32'), ('stop','int32'), ('score', 'float64')])
    return np.array(bursts, dtype=dt).view(np.recarray)
