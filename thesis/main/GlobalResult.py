from abc import ABC, abstractmethod
import fenics as fcs
import numpy as np
from numbers import Number
import mpi4py.MPI as MPI
import os

class ResultEmptyError(Exception):pass

class GlobalResult(ABC):

    def __init__(self, directory_path: str, field_quantity: str):

        self.path:str = directory_path
        self.field_quantity: str = field_quantity

    @abstractmethod
    def save(self):pass

    @abstractmethod
    def load(self):pass

    @abstractmethod
    def get(self):pass

    @abstractmethod
    def set(self):pass

class ScalarResult(GlobalResult):

    def __init__(self, *args):
        super().__init__(*args)
        self.u: float = None

    def set(self, u: float):

        assert isinstance(u,Number)

        self.u = u

    def get(self):

        if self.u is None: raise ResultEmptyError("The result value was not set!")

        return self.u

    def save(self, time_index: int):

        file_name = "field_{fq}.npy".format(fq=self.field_quantity)
        partial_path = "sol/"+file_name

        file = os.path.join(self.path,partial_path)
        if not os.path.exists(os.path.join(self.path,"sol/")):
            os.makedirs(os.path.join(self.path,"sol/"))

        if os.path.exists(file) and time_index > 1:
            u = np.load(file)
            u = np.insert(u,len(u),self.u)
        else:
            u = np.array([self.u])

        np.save(file,u)

        return None, partial_path, type(self)

    def load(self, time_index: int):

        file_name = "field_{fq}.npy".format(fq=self.field_quantity)
        partial_path = "sol/" + file_name

        file = os.path.join(self.path, partial_path)
        if os.path.exists(file):
            self.u =  np.load(file)[time_index-1]# todo assumes 1 as start index
        else:
            raise FileNotFoundError


class ScalarFieldResult(GlobalResult):


    def __init__(self,  *args):
        super().__init__(*args)
        self.u: fcs.Function = None

    def set(self, u):

        assert isinstance(u,fcs.Function)
        self.u = u

    def get(self):

        if self.u is None: raise ResultEmptyError("The result value was not set!")
        return self.u

    def save(self, time_index:int):

        file_name = "field_{fq}".format(fq=self.field_quantity)
        distplot = "sol/distplot/{fn}_{ti}_distPlot.h5".format(fn=file_name, ti=str(time_index))
        sol = "sol/{fn}_{ti}.xdmf".format(fn=file_name, ti=str(time_index))
        u = self.get()

        u.rename(self.field_quantity, self.field_quantity)
        with fcs.HDF5File(fcs.MPI.comm_world, os.path.join(self.path,distplot), "w") as f:
            f.write(u, self.field_quantity)
        with fcs.XDMFFile(fcs.MPI.comm_world, os.path.join(self.path, sol)) as f:
            f.write(u, time_index)
        return (distplot, sol, type(self))

    def load(self, time_index, mesh_path):

        comm = MPI.COMM_WORLD


        file_name = "field_{fq}".format(fq=self.field_quantity)
        distplot = "sol/distplot/{fn}_{ti}_distPlot.h5".format(fn=file_name, ti=str(time_index))
        sol = "sol/{fn}_{ti}.xdmf".format(fn=file_name, ti=str(time_index))

        mesh = fcs.Mesh()
        with fcs.XDMFFile(mesh_path ) as f:
            f.read(mesh)

        function_space = fcs.FunctionSpace(mesh, "P", 1)

        u: fcs.Function = fcs.Function(function_space)
        with fcs.HDF5File(comm, os.path.join(self.path,distplot), "r") as f:
            f.read(u, "/" + self.field_quantity)

        self.u = u