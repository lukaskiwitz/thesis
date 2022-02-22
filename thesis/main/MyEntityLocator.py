import json
import os
from abc import ABC, abstractmethod

import numpy as np
from scipy.spatial import cKDTree

from thesis.main.Entity import Entity, Cell
from thesis.main.EntityType import EntityType, CellType
from thesis.main.ParameterSet import ParameterSet
from thesis.main.SimComponent import SimComponent
from thesis.main.my_debug import message


# module_logger = logging.getLogger(__name__)


class MyEntityLocator(ABC, SimComponent):

    def __init__(self):
        super().__init__()

    def get_entity_list(self, entity_types: [EntityType], global_p: ParameterSet, path: str, mesh_folder_path,
                        overwrite_cache=True) -> [
        Entity]:

        if mesh_folder_path is not None:
            fp = os.path.join(mesh_folder_path, "entity_id_to_pos_map.json")
        else:
            fp = os.path.join(path, "entity_id_to_pos_map.json")

        os.makedirs("/".join(fp.split("/")[:-1]), exist_ok=True)

        if os.path.exists(fp) and not overwrite_cache:
            with open(fp, "r") as f:
                id_to_pos_map = {int(k): v for k, v in json.load(f).items()}

            message("loaded id-to-pos map from file for {n} entities".format(n=len(id_to_pos_map)), self.logger)
            c = 0

            entity_list = []
            for id, pos in id_to_pos_map.items():
                r = entity_types.p.get_physical_parameter("rho", "rho").get_in_sim_unit()
                cell = Cell(pos, r, [])
                cell.set_cell_type(entity_types, None, 0)
                cell.id = int(id)
                entity_list.append(cell)

            for e in entity_list:
                if e.id in id_to_pos_map.keys():
                    c = c + 1
                    e.center = id_to_pos_map[e.id]

            message("restored id-to-pos map for {n} entities".format(n=c), self.logger)

        else:
            entity_list = self._get_entity_list(entity_types, global_p, path)
            id_to_pos_map = {e.id: list(e.center) for e in entity_list}
            message("saving id-to-pos map from file for {n} entities".format(n=len(id_to_pos_map)), self.logger)
            with open(fp, "w") as f:
                json.dump(id_to_pos_map, f)

        return entity_list

    @abstractmethod
    def _get_entity_list(self, entity_types: [EntityType], global_p: ParameterSet, path: str) -> [Entity]:
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

    def _get_entity_list(self, cell_type: CellType, global_p: ParameterSet, path: str) -> [Cell]:
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
            cell.id = i
            cell_list.append(cell)

        return cell_list


class MyCellGridLocator(MyEntityLocator):
    """
    Locator class to place cells in a primitive cubic lattice.
    Dimensions are taken from "geometry" parameter collection.
    """

    def __init__(self):
        super(MyCellGridLocator, self).__init__()

    def _get_entity_list(self, cell_type: CellType, global_p: ParameterSet, path: str) -> [Cell]:

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
            cell.id = i
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
        super().__init__()

    def correct_overlaps(self, positions, iterations, rho, penalty):
        for i in range(iterations):
            len_neighbors = 0
            tree = cKDTree(positions)
            for j in range(len(positions)):
                neighbors = tree.query(positions[j], 3, distance_upper_bound=(rho + penalty) * 2, n_jobs=1)[1][1:]
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

    def random_3D_vector(self, length):
        """
        Generates a random 3D unit vector (direction) with a uniform spherical distribution
        Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
        :return:
        """
        phi = np.random.uniform(0, np.pi * 2)
        costheta = np.random.uniform(-1, 1)

        theta = np.arccos(costheta)
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        return np.array([x, y, z]) * length

    def _get_entity_list(self, cell_type: CellType, global_p: ParameterSet, path: str) -> [Cell]:
        assert issubclass(type(cell_type), CellType)
        r = cell_type.p.get_physical_parameter("rho", "rho").get_in_sim_unit()

        self.cell_pos = self.get_random_pos(r, 3, global_p)

        cell_list = []
        for i, p in enumerate(self.cell_pos):
            assert r is not None and r > 0
            cell = Cell(p, r, [])
            cell.set_cell_type(cell_type, None, 0)
            cell.id = i
            cell_list.append(cell)

        return cell_list

    def get_random_pos(self, rho, penalty, global_p):

        xx = self.make_grid(global_p, "x_grid")
        yy = self.make_grid(global_p, "y_grid")
        zz = self.make_grid(global_p, "z_grid")

        length = len(xx) * len(yy) * len(zz)

        x, y, z = np.meshgrid(xx, yy, zz)

        positions = np.array([np.ravel(x), np.ravel(y), np.ravel(z)]).T

        steps = global_p.get_misc_parameter("steps", "geometry").get_in_sim_unit(type=int)
        step_size = global_p.get_misc_parameter("step_size", "geometry").get_in_sim_unit(type=float)

        f = lambda x: x + self.random_3D_vector(step_size)
        for i in range(steps):
            positions = np.array([f(xi) for xi in positions])

        iterations = 1000
        positions = self.correct_overlaps(positions, iterations, rho, penalty)

        return positions


class MyBridsonCellLocator(MyCellGridLocator):
    """
    Implements bridson sampling to generate blue noise distribution
    https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
    """

    def __init__(self):
        super().__init__()

    def _get_entity_list(self, cell_type: CellType, global_p: ParameterSet, path: str) -> [Cell]:

        assert issubclass(type(cell_type), CellType)
        r = cell_type.p.get_physical_parameter("rho", "rho").get_in_sim_unit()

        self.cell_pos = self.get_random_pos(r, 3, global_p)

        cell_list = []
        for i, p in enumerate(self.cell_pos):
            assert r is not None and r > 0
            cell = Cell(p, r, [])
            cell.set_cell_type(cell_type, None, 0)
            cell.id = i
            cell_list.append(cell)

        return cell_list

    def get_random_pos(self, rho, penalty, global_p):

        from thesis.cellBehaviourUtilities.bridson_sampling import bridson
        # steps = global_p.get_misc_parameter("steps", "geometry")
        # steps = steps.get_in_sim_unit(type=int) if steps is not None else 30
        steps = 30
        distance = global_p.get_misc_parameter("distance", "geometry").get_in_sim_unit(type=int)

        xx = self.make_grid(global_p, "x_grid")
        yy = self.make_grid(global_p, "y_grid")
        zz = self.make_grid(global_p, "z_grid")

        BB = [xx, yy, zz]

        BP1 = [np.min(i) for i in BB if len(i) > 1]
        BP2 = [np.max(i) for i in BB if len(i) > 1]

        X = bridson(steps, BP1, BP2, density_function=lambda x: distance)

        positions = np.ndarray(shape=(len(X), 3))

        c = 0
        for i, x in enumerate(BB):
            if len(x) > 1:
                positions[:, i] = X[:, c]
                c = c + 1
            else:
                positions[:, i] = np.ones(shape=(len(X))) * x[0]
        return positions
