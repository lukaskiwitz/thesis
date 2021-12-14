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


def make_clusters(cell_grid_positions, apcs, fractions, cluster_strengths):
    """
    Generates clusters of cell types around an arbitrary number of apcs positions.
    Works in pseudo 2D.

    :param cell_grid_positions: Iterable with shape (N,3). The list of possible position for cells
    :param apcs: Iterable with shape (N,3).  The list of clustering centers. Must not be on cell grid
    :param fractions: Iterable with shape (N_types,). List of fractions for each cell type. Sum must be <= 1.
    :param cluster_strengths: Iterable with shape (N_types). Clustering strength corresponding to each entry in fractions.
        0 -> max. clustering; 1 -> homogenous distribution

    :return: Iterable with shape (N,). Indicates cell type for each grid position. fractions[0] corresponds to index 1,
        because 0 indicates background. Values go from 0 to N_types + 1.
    """

    np.random.seed(0)

    individual = True
    bws = np.ones(shape=(len(fractions, ))) * np.max(cell_grid_positions)

    number_of_cell_types = len(fractions)
    cell_types = np.zeros(shape=(cell_grid_positions.shape[0], number_of_cell_types), dtype=bool)
    cell_type = np.zeros(shape=(cell_grid_positions.shape[0],), dtype=int)

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
        l = int(len(cell_grid_positions) * N)  # Number of cells for this type

        draw_length = int(q * len(ai))
        draw_length = draw_length if draw_length > l else l

        draws = np.random.choice(draw_length, (l,), replace=False) if q > 0 else np.arange(0, l, dtype=int)
        cell_types[ai[draws], t] = t + 1

    # This loop looks for conflicts and resovles them
    for i in np.random.choice(len(cell_types), len(cell_types), replace=False):

        nz = np.where(cell_types[i])[0]
        if len(nz) == 0:
            continue
        else:

            ti = np.random.choice(nz - 1, 1, replace=False)[0]
            cell_type[i] = nz[ti] + 1
            empty_position_indices = np.where(np.all(np.invert(cell_types), axis=1))[0]
            for o in np.delete(nz, ti):
                for p in sorted_positions[:, o]:
                    if p in empty_position_indices:
                        cell_type[p] = o + 1
                        empty_position_indices = np.delete(empty_position_indices,
                                                           np.where(empty_position_indices == p))
                        cell_types[p, o] = True
                        break

    return cell_type
