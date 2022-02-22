from copy import copy

import numpy as np
from scipy.spatial import cKDTree


def hyper_angles_to_cartesian(r, Phi):
    n = len(Phi) + 1

    X = []
    for i in range(n):

        assert i <= n
        i = i + 1
        if i == 1:
            X.append(r * np.cos(Phi[0]))
        elif i == n:
            X.append(r * np.product([np.sin(Phi[0:n - 1])]))
        else:
            X.append(r * np.product([np.sin(Phi[0:i - 1])]) * np.cos(Phi[i - 1]))
    return X


def bridson(k, BP1, BP2, density_function=lambda x: 50):
    BP1 = np.array(BP1)
    BP2 = np.array(BP2)

    assert len(BP1) == len(BP2)
    N = len(BP1)
    assert N >= 2

    xc = np.array([BP1, BP2]).mean(axis=0)
    x_0 = copy(xc)

    active_list = [x_0]
    X = [x_0]

    while len(active_list) > 0:

        tree = cKDTree(X)

        ai = np.random.randint(0, len(active_list))
        x_i = active_list[ai]
        r = density_function(x_i)
        for i in range(k):

            Phi = list(np.random.uniform(0, 180, (N - 2,))) + list(np.random.uniform(0, 360, (1,)))
            R = np.random.uniform(r, 2 * r)
            Xnew = hyper_angles_to_cartesian(R, Phi)
            x_new = x_i + Xnew

            if len(tree.query_ball_point(x_new, r, p=2)) > 0:
                continue
            elif np.all(x_new > BP1) and np.all(x_new < BP2):
                active_list.append(x_new)
                X.append(x_new)
                break

        if i == k - 1:
            del active_list[ai]
    return np.array(X)
