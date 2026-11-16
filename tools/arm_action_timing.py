"""Sample the outward and return trace consistently across replay exports."""

import numpy as np


def geometric_replay_indices(trace_length):
    peak = trace_length // 2
    outward = [round(position) for position in np.linspace(0, peak, 22)]
    return [*outward, *[peak] * 7, *outward[-2::-1]]
