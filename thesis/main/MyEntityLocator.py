import numpy as np

from thesis.main.Entity import Entity, Cell
from thesis.main.EntityType import EntityType, CellType
from thesis.main.ParameterSet import ParameterSet

class MyEntityLocator:

    def __init__(self):
        pass

    def get_entity_list(self, entity_types: [EntityType], global_p: ParameterSet) -> [Entity]:
        return None

class MyCellListLocator(MyEntityLocator):

    def __init__(self, cell_pos, cell_types):


        assert np.array(cell_pos).shape[1] == 3
        assert len(cell_pos) == len(cell_types)

        self.cell_pos = cell_pos
        self.cell_types = cell_types


    def get_entity_list(self, cell_type: CellType, global_p: ParameterSet) -> [Cell]:

        assert issubclass(type(cell_type), CellType)

        cell_list = []
        for i, p in enumerate(self.cell_pos):
            r = self.cell_types[i].p.get_physical_parameter("rho", "rho").get_in_sim_unit()
            assert r is not None and r > 0

            cell = Cell(p, r, [])
            cell.set_cell_type(self.cell_types[i], None)
            cell_list.append(cell)

        return cell_list

class MyCellGridLocator(MyEntityLocator):

    def __init__(self):

        pass

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
            cell.set_cell_type(cell_type, None)
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
