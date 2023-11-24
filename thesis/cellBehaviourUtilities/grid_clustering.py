import KDEpy as kde
import numpy as np


def get_apc_positions(cell_grid_positions, no_apcs=5):
    """
    helper function to uniformly distribute apcs in cell_grid axis aligned bounding box.


    :param cell_grid_positions: Iterable with shape (N,3). The list of possible position for cells
    :param no_apcs: number of apcs
    :return: Iterable with shape (N,3).  The list of clustering centers. Must not be on cell grid
    """

    f = lambda x: np.random.uniform(min(x), max(x), (no_apcs))

    apcs = np.array([
        f(cell_grid_positions[:, 0]),
        f(cell_grid_positions[:, 1]),
        f(cell_grid_positions[:, 2])
    ]).T

    return apcs


def get_cell_grid_positions(x_dim, y_dim, z_dim, distance=20):
    """
    Helper function to produce cell grid from grid dimensions

    :param x_dim:
    :param y_dim:
    :param z_dim:
    :return: Iterable with shape (N,3). The list of possible position for cells
    """

    x = np.arange(0, x_dim, distance)
    y = np.arange(0, y_dim, distance)
    z = np.arange(0, z_dim, distance)

    xx, yy, zz = np.meshgrid(x, y, z)

    xx = np.ravel(xx)
    yy = np.ravel(yy)
    zz = np.ravel(zz)

    cells = np.array([xx, yy, zz]).T

    return cells


def make_clusters(cell_grid_positions, apcs, fractions_dict, cluster_strengths, seed=0):
    """
    Generates clusters of cell types around an arbitrary number of apc positions.
    Also works in pseudo 2D.

    :param cell_grid_positions: Iterable with shape (N,3). The list of possible position for cells
    :param apcs: Iterable with shape (N,3).  The list of clustering centers. Must not be on cell grid
    :param fractions: Iterable with shape (N_types,). List of fractions for each cell type. Sum must be <= 1.
    :param cluster_strengths: Iterable with shape (N_types). Clustering strength corresponding to each entry in fractions.
        0 -> max. clustering; 1 -> homogenous distribution

    :return: Iterable with shape (N,). Indicates cell type for each grid position. fractions[0] corresponds to index 1,
        because 0 indicates background. Values go from 0 to N_types + 1.
    """

    # np.random.seed(seed)
    fractions = list(fractions_dict.values())
    individual = True
    bws = np.ones(shape=(len(fractions, ))) * np.max(cell_grid_positions)
    cluster_strengths = [1e-100 if cs == 0 else cs for cs in cluster_strengths]
    number_of_cell_types = len(fractions)
    cell_types = np.zeros(shape=(cell_grid_positions.shape[0], number_of_cell_types + 1), dtype=bool)
    sorted_positions = np.zeros((len(cell_grid_positions), len(fractions)), dtype=int)

    for t, N in enumerate(fractions):

        if individual:
            p_list = np.zeros((len(apcs), cell_grid_positions.shape[0]))
            for i in range(len(apcs)):
                density = kde.TreeKDE("tri", bw=bws[t]).fit(np.array([apcs[i, :]]))
                p_list[i] = density.evaluate(cell_grid_positions)
            p = np.max(p_list, axis=0)
        else:
            density = kde.TreeKDE("tri", bw=bws[t]).fit(apcs)
            p = density.evaluate(cell_grid_positions)

        ai = np.flip(np.argsort(p))
        sorted_positions[:, t] = ai

        q = cluster_strengths[t]
        l = int(np.around(len(cell_grid_positions) * N))  # Number of cells for this type
        a = 0.1

        def exp_choice(q, l, s):

            if l == 0:
                return []

            b = (1 / a - 1)
            expected_value = 1 / (a * (q)) - b

            def f(x, p, c):
                return p * np.power((1 - p), (c * x))

            x = np.arange(0, s, 1)
            p = f(x, 1 / (expected_value + 1), 1 / (1 * l * a))
            p = (1 / np.sum(p)) * p  # normalization for truncation error for large mean_distance

            cutoff_index = np.max(np.nonzero(p)) + 1
            x = x[:cutoff_index]
            p = p[:cutoff_index]

            draws = np.random.choice(x, (l,), replace=False, p=p)
            return draws

        draws = exp_choice(q, l, len(cell_grid_positions))

        cell_types[ai[draws], t + 1] = t + 1

    def resolve_conflict(conflict_row, cluster_strengths, fractions_dict_keys=()):

        draw_list = conflict_row
        p = np.array([cluster_strengths[o] for o in conflict_row])
        if np.all(p == 0):
            p = np.ones(p.shape)

        p = (1 / np.sum(p)) * p
        if "Tsec" in np.array(fractions_dict_keys)[conflict_row]:
            ti = np.argwhere(np.array(fractions_dict_keys)[conflict_row] == "Tsec")[0]
        else:
            not_index = np.random.choice(draw_list, 1, replace=False, p=p)[0]
            ti = np.argwhere(draw_list == not_index)[0, 0]

        winner_index = np.argwhere(conflict_row == draw_list[ti])[0, 0]
        remaining_indices = np.delete(conflict_row, winner_index)

        return winner_index, remaining_indices

    # This loop looks for conflicts and resovles them

    positions_to_check = np.where(np.count_nonzero(cell_types[:, 1:], axis=1) > 1)[0]

    counter = 0
    while len(positions_to_check) > 0:

        if counter > 1e4:
            raise StopIteration

        counter += 1
        i = np.random.choice(len(positions_to_check), 1, replace=False)[0]
        ii = positions_to_check[i]
        nz = np.where(cell_types[ii, 1:])[0]  # check where two celltypes would be chosen
        assert len(nz) > 1  # make sure there is a conflict
        winner_index, remaining_indices = resolve_conflict(nz, cluster_strengths, list(
            fractions_dict.keys()))  # winner_index is the index of nz,
        # remaining_indices are the entries of nz, so the indices of cell_types[ii, 1:]

        for o in remaining_indices:
            cell_types[ii, o + 1] = False

            def g(a, b):
                r = np.random.uniform(0, 1, a.shape)
                return (b - a) > r

            candidates = np.argwhere(
                np.all(
                    np.logical_or(
                        g(np.tile(cluster_strengths, (cell_types.shape[0], 1)), cluster_strengths[o]),
                        cell_types[:, 1:][sorted_positions[:, 0]] == False)
                    , axis=1)
            )[:, 0]

            assert (len(positions_to_check) == 1) or (len(candidates > 0))
            if len(candidates) == 0:
                break

            p = exp_choice(cluster_strengths[o], 1, len(candidates))[0]
            cell_types[sorted_positions[candidates[p], o], o + 1] = True

        positions_to_check = np.where(np.count_nonzero(cell_types[:, 1:], axis=1) > 1)[0]

    cell_types[:, 0] = np.invert(np.any(cell_types[:, 1:], axis=1))
    return np.where(cell_types)[1]
