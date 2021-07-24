import os
from thesis.main.PostProcess import PostProcessor, PostProcessComputation
import fenics as fcs
from thesis.main.ParameterSet import PhysicalParameter, PhysicalParameterTemplate
from parameters import path, boundary
os.environ["LOG_PATH"] = path


class errornorm(PostProcessComputation):
    """Defines a custom computation"""


    def __init__(self):
        """this name will appear in output dataframe"""
        self.name = "error"

    def __call__(self, *args, **kwargs):

        return self.high_density(*args,**kwargs)

    def high_density(self, u, grad, c_conv, grad_conv, mesh_volume, V = None, cell_df = None, path = None, **kwargs):
        from scipy.constants import N_A
        ule = -6

        t = {
        "R": PhysicalParameterTemplate(PhysicalParameter("R", 0, to_sim=N_A ** -1 * 1e9)),
        "k_on": PhysicalParameterTemplate(PhysicalParameter("k_on", 111.6, to_sim=1e15 / 60 ** 2, is_global=True)),
        "q": PhysicalParameterTemplate(PhysicalParameter("q", 0, to_sim=N_A ** -1 * 1e9)),
        "D": PhysicalParameterTemplate(PhysicalParameter("D", 10, to_sim=1, is_global=True)),
        "kd": PhysicalParameterTemplate(PhysicalParameter("kd", 0.1, to_sim=1 / (60 ** 2), is_global=True)),
        "rho": PhysicalParameterTemplate(PhysicalParameter("rho", 0, to_sim=10 ** (-6 - ule))),

        }

        R = t["R"](1e2).get_in_sim_unit()
        Resp = t["R"](1e2).get_in_sim_unit()

        q = t["q"](10).get_in_sim_unit()
        rho = t["rho"](5).get_in_sim_unit()
        L = t["rho"](15).get_in_sim_unit()

        kon = t["k_on"](116).get_in_sim_unit()
        D = t["D"](10).get_in_sim_unit()

        r = "(pow(pow(x[0],2) + pow(x[1],2) + pow(x[2],2),0.5))"

        enum = "4*pi *D*{r}*(L + rho) + kon*(L- {r} + rho) * N*Resp".format(r=r)
        denum = "kon*L*R*N*Resp + 4*pi*D*rho*(L+rho)*(R + N*Resp)"

        v = fcs.Expression("((q*rho/(kon*{r}) * (({enum})/({denum}))))".format(r = r, enum = enum, denum = denum)
                           , degree=1,q=q, rho=rho, kon=kon, R=R, D=D, N =12, Resp = Resp, L = L)
        s = fcs.interpolate(v,V)
        error = fcs.errornorm(u,s,norm_type="l2", degree_rise=1)

        return error * c_conv

class n_nodes(PostProcessComputation):
    """Defines a custom computation"""


    def __init__(self):
        """this name will appear in output dataframe"""
        self.name = "nodes"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, V = None, cell_df = None, **kwargs):

        return  len(u.vector())


pp = PostProcessor(path)
pp.unit_length_exponent = -6
pp.computations.append(errornorm())
pp.computations.append(n_nodes())
pp.run_post_process(2)