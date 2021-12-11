import logging
from abc import ABC, abstractmethod

import numpy as np

from thesis.main.Entity import Entity, Cell
from thesis.main.EntityType import EntityType, CellType
from thesis.main.ParameterSet import ParameterSet
from thesis.main.SimComponent import SimComponent

module_logger = logging.getLogger(__name__)


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
                assert isinstance(self.cell_types[i],CellType)
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
