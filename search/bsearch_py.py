"""
Functions for burst search.
"""

import numpy as np


def burstsearch_py(signal, m, threshold):
    """Sliding window burst search. Pure python version.
    signal is barray (burst array: F[n+2]-2F[n+1]+F[n])
    Returns:
        Record array of burst data, one element per burst.
        Each element is a composite data type containing
        burst start, burst stop and burst score.
    """
    bursts = []
    in_burst = False
    score = signal[:m].sum()
    deltascore = signal[m:] - signal[:-m]
    for i, delta in enumerate(deltascore):
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

def burstarray(fmean_th, normalize=True):
    Favg = np.mean(fmean_th)
    a = fmean_th[2:] - 2*fmean_th[1:-1] + fmean_th[:-2]
    if normalize:
        a /= Favg
    return a[::2]

def burstmask(nframes, period, bursts):
    burstmask = np.zeros(nframes, dtype=bool)
    for start, stop, score in bursts:
        burstmask[start*period : stop*period] = True
    return burstmask


def burstmask2(bursts, barray):
    burstmask = np.zeros(len(barray), dtype=bool)
    for start, stop, score in bursts:
        burstmask[start : stop] = True
