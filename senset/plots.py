import numpy as np
import seaborn as sns


def pairwise_overlap_heatmap(data, ticklabels=None):
    """Plots size of overlap for every pair of lists in data."""
    N = len(data)
    hm = np.zeros((N, N), dtype=int) - 1
    for i in range(N):
        for j in range(i, N):
            hm[i, j] = np.intersect1d(data[i], data[j]).size

    kwargs = {}
    if ticklabels is not None:
        kwargs['xticklabels'] = ticklabels
        kwargs['yticklabels'] = ticklabels

    ax = sns.heatmap(
        hm,
        annot=True,
        fmt='g',
        mask=hm == -1,
        square=True,
        **kwargs,
    )
    return ax
