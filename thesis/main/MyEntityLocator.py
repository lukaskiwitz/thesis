from abc import ABC, abstractmethod

import numpy as np
from scipy.spatial import cKDTree, distance_matrix

from thesis.main.Entity import Entity, Cell
from thesis.main.EntityType import EntityType, CellType
from thesis.main.ParameterSet import ParameterSet
from thesis.main.SimComponent import SimComponent


# module_logger = logging.getLogger(__name__)


class MyEntityLocator(ABC, SimComponent):

    def __init__(self):
        super(MyEntityLocator, self).__init__()

    @abstractmethod
    def get_entity_list(self, entity_types: [EntityType], global_p: ParameterSet) -> [Entity]:
        return None


class MyCellListLocator(MyEntityLocator):
    """
    Locator class to arbitrarily place a set of cells in the simulation.

    """

    def __init__(self, cell_pos, cell_types):
        """
        Accepts a (n,3) list of cells positions and a list of matching cell types.
        If len(cell_pos) > len(cell_types), the last element in cell_types is used for the remaining cells.
        Note: Cell types can be changed later (in pre_scan, pre_step, etc), but a dummy is needed to get the
        the cell radius parameter (rho) before meshing.

        :param cell_pos: list of entity positions
        :param cell_types: corresponding list of entity templates.
        """

        assert np.array(cell_pos).shape[1] == 3
        assert len(cell_pos) > 0

        super(MyCellListLocator, self).__init__()

        self.cell_pos = cell_pos
        self.cell_types = cell_types

    def get_entity_list(self, cell_type: CellType, global_p: ParameterSet) -> [Cell]:
        assert issubclass(type(cell_type), CellType)

        cell_list = []
        for i, p in enumerate(self.cell_pos):
            if len(self.cell_types) > i:
                assert isinstance(self.cell_types[i], CellType)
                cell_type = self.cell_types[i]
            else:
                cell_type = self.cell_types[-1]

            r = cell_type.p.get_physical_parameter("rho", "rho").get_in_sim_unit()
            assert r is not None and r > 0
            cell = Cell(p, r, [])
            cell.set_cell_type(cell_type, None, 0)

            cell_list.append(cell)

        return cell_list


class MyCellGridLocator(MyEntityLocator):
    """
    Locator class to place cells in a primitive cubic lattice.
    Dimensions are taken from "geometry" parameter collection.
    """

    def __init__(self):
        super(MyCellGridLocator, self).__init__()

    def get_entity_list(self, cell_type: CellType, global_p: ParameterSet) -> [Cell]:

        assert issubclass(type(cell_type), CellType)

        xx = self.make_grid(global_p, "x_grid")
        yy = self.make_grid(global_p, "y_grid")
        zz = self.make_grid(global_p, "z_grid")

        x, y, z = np.meshgrid(xx, yy, zz)

        points = np.array([np.ravel(x), np.ravel(y), np.ravel(z)])
        cell_list = []
        for i, p in enumerate(points.T):
            r = cell_type.p.get_physical_parameter("rho", "rho").get_in_sim_unit()
            assert r is not None and r > 0

            cell = Cell(p, r, [])
            cell.set_cell_type(cell_type, None, 0)
            cell_list.append(cell)

        return cell_list

    def make_grid(self, p, key: str):

        margin = p.get_misc_parameter("margin", "geometry").get_in_sim_unit(type=float)
        distance = p.get_misc_parameter("distance", "geometry").get_in_sim_unit(type=float)
        grid = p.get_misc_parameter(key, "geometry")

        if grid is not None:
            grid = grid.get_in_sim_unit(type=float)
            return np.round(np.arange(margin, grid, distance), 2)
        else:
            return [0]


class MyRandomCellLocator(MyCellGridLocator):
    """
    Randomizes grid positions into collision free uniform distribution of cells
    """

    def __init__(self):

        pass

    def get_entity_list(self, cell_type: CellType, global_p: ParameterSet) -> [Cell]:
        assert issubclass(type(cell_type), CellType)
        r = cell_type.p.get_physical_parameter("rho", "rho").get_in_sim_unit()

        self.cell_pos = self.get_random_pos(r, 3, global_p)

        cell_list = []
        for i, p in enumerate(self.cell_pos):
            assert r is not None and r > 0
            cell = Cell(p, r, [])
            cell.set_cell_type(cell_type, None, 0)

            cell_list.append(cell)

        return cell_list

    def get_random_pos(self, rho, penalty, global_p):

        xx = self.make_grid(global_p, "x_grid")
        yy = self.make_grid(global_p, "y_grid")
        zz = self.make_grid(global_p, "z_grid")

        length = len(xx) * len(yy) * len(zz)
        x = np.random.uniform(np.min(xx), np.max(xx), (1, length))
        y = np.random.uniform(np.min(yy), np.max(yy), (1, length))
        z = np.random.uniform(np.min(zz), np.max(zz), (1, length))

        positions = np.vstack(list(zip(x.ravel(), y.ravel(), z.ravel())))

        iterations = 1000

        for i in range(iterations):
            len_neighbors = 0
            tree = cKDTree(positions)
            for j in range(len(positions)):
                neighbors = tree.query(positions[j], 2, distance_upper_bound = (rho + penalty) * 2, n_jobs = 1)[1][1:]
                neighbors = neighbors[neighbors != len(positions)]
                length = len(neighbors)
                len_neighbors += length
                if length != 0:
                    cell_1 = np.mean(positions[neighbors], axis=0)
                    cell_2 = positions[j]
                    factor = 1.01 if length > 4 else 0.1
                    c_c_v = ((cell_2 - cell_1) / 2) * factor
                    positions[j] += + c_c_v

            if len_neighbors == 0:
                break


        return positions

