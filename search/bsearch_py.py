"""
Functions for burst search.
"""

import numpy as np


def burstsearch_py(signal, m, threshold):
    """Sliding window burst search. Pure python version.

    Returns:
        Record array of burst data, one element per burst.
        Each element is a composite data type containing
        burst start, burst stop and burst score.
    """
    bursts = []
    in_burst = False
    score = signal[:m].sum()
    deltasignal = signal[m:] - signal[:-m]
    for i, delta in enumerate(deltasignal):
        if score > threshold:
            if not in_burst:
                # index of first cycle in burst
                start = i
                in_burst = True
        elif in_burst:
            # index of last cycle in burst
            stop = i + m - 2
            totalscore = signal[start:stop + 1].sum()
            bursts.append((start, stop, totalscore))
            in_burst = False
        score += delta

    # Create numpy recarray
    dt = np.dtype([('start','int32'), ('stop','int32'), ('score', 'float64')])
    return np.array(bursts, dtype=dt).view(np.recarray)




