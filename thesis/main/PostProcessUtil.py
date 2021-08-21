import fenics as fcs
import numpy as np


def get_concentration_conversion(unit_length_exponent: int):
    assert isinstance(unit_length_exponent, int)
    return float(10 ** (-1 * (unit_length_exponent * 3 + 3)))


def get_gradient_conversion(unit_length_exponent: int, target_exponent = -6):
    #converts gradient to nM/um for target_exponent = -6.
    assert isinstance(unit_length_exponent, int)
    # exp = (-1 * (unit_length_exponent * 3 + 3)) - 1
    # exp -= (6 + unit_length_exponent)  # to nM/um

    return float(10 ** (-3-4 * unit_length_exponent + target_exponent))


def get_mesh_volume(mesh):
    sum = 0
    for cell in fcs.cells(mesh):
        sum += cell.volume()
    return sum


def get_rectangle_plane_mesh(u, res = (200,200)):

    coords = u.function_space().mesh().coordinates()
    coords_t = np.transpose(coords)

    x_limits = (
        np.min(coords_t[0]), np.max(coords_t[0])
    )
    y_limits = (
        np.min(coords_t[1]), np.max(coords_t[1])
    )
    rec_mesh = fcs.RectangleMesh(fcs.Point(x_limits[0], y_limits[0]), fcs.Point(x_limits[1], y_limits[1]), res[0], res[1])

    return rec_mesh, [x_limits,y_limits]