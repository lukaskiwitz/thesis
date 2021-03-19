import numpy as np
from scipy.constants import N_A
from thesis.main.ParameterSet import MiscParameter, ParameterCollection, ScannablePhysicalParameter, ParameterSet


class updateState():

    def __init__(self,t,geometry,templates,Tsec_distribution_array,Th_distribution_array,offset=0):
        self.Tsec_distribution_array = Tsec_distribution_array
        self.Th_distribution_array = Th_distribution_array
        self.geometry = geometry
        self.offset = offset
        self.templates = templates
        self.t = t
        self.drawn_fraction = 0

    def get_offset_ids(self,sc):
        no_of_cells = len(sc.entity_list)
        offset = self.offset
        try:
            a = self.geometry["z_grid"]
            dims = 3
            del a
        except KeyError:
            dims = 2
        print("dims = ", dims)
        if dims == 3:
            cube_size = round(np.cbrt(no_of_cells), 0)
        else:
            cube_size = round(np.sqrt(no_of_cells), 0)
        xr = np.arange(0, cube_size, 1)
        yr = np.arange(0, cube_size, 1)

        if dims == 3:
            zr = np.arange(0, cube_size, 1)
            z_offset = np.array([zr[:offset], zr[-offset:]]).flatten()
        else:
            zr = [None]
            z_offset = []

        x_offset = np.array([xr[:offset], xr[-offset:]]).flatten()
        y_offset = np.array([yr[:offset], yr[-offset:]]).flatten()

        anti_draws = []
        counter = 0
        for x in xr:
            for y in yr:
                for z in zr:
                    if x in x_offset or y in y_offset or z in z_offset:
                        anti_draws.append(counter)
                    counter += 1
        return anti_draws


    def get_draws(self,sc,exluded_ids, fraction, at_least_one = True):
        possible_draws = np.setdiff1d(range(len(sc.entity_list)), exluded_ids)
        print("possible draws: ", len(possible_draws))
        fraction = fraction + self.drawn_fraction
        if np.isclose(1.0,fraction):
            self.drawn_fraction = 0
        else:
            self.drawn_fraction = fraction
        print("fraction = ", fraction)
        amount_of_draws = int(round(len(possible_draws) * fraction))
        if at_least_one == True:
            if amount_of_draws == 0:
                amount_of_draws = 1
                print("set to 1")
        try:
            draws = np.random.choice(possible_draws, amount_of_draws, replace=False)
        except ValueError:
            draws = np.random.choice(possible_draws, amount_of_draws - 1, replace=False)
        return draws


    def set_parameters(self,sc,Tsec_draws, q_il2_sum):
        t_R_start = self.templates["R_start"]
        t_pSTAT5 = self.templates["pSTAT5"]
        t_EC50 = self.templates["EC50"]
        for i, e in enumerate(sc.entity_list):
            e.p.add_parameter_with_collection(t_pSTAT5(0, in_sim=False))
            e.p.add_parameter_with_collection(t_EC50(0, in_sim=False))
            if len(Tsec_draws) != 0 and e.change_type == "Tsec":
                e.p.get_physical_parameter("q", "IL-2").set_in_sim_unit(q_il2_sum / len(Tsec_draws))
                e.p.add_parameter_with_collection(t_R_start(1e2, in_sim=False))
            elif e.change_type == "Th":
                e.p.add_parameter_with_collection(t_R_start(1e4, in_sim=False))
            elif e.change_type == "Tnaive":
                e.p.add_parameter_with_collection(t_R_start(1e2, in_sim=False))
            elif e.change_type == "blank":
                e.p.add_parameter_with_collection(t_R_start(0, in_sim=False))


    def set_R_lognorm_parameters(self,sc,Tsec_draws,q_il2_sum):
        var = sc.p.get_physical_parameter("sigma", "IL-2").get_in_sim_unit()
        for i, e in enumerate(sc.entity_list):
            t_R_start = self.templates["R_start"]
            t_pSTAT5 = self.templates["pSTAT5"]
            t_EC50 = self.templates["EC50"]
            e.p.add_parameter_with_collection(t_EC50(0, in_sim=False))
            if len(Tsec_draws) != 0 and e.change_type == "Tsec":
                e.p.get_physical_parameter("q", "IL-2").set_in_sim_unit(q_il2_sum / len(Tsec_draws))
                e.p.add_parameter_with_collection(t_R_start(1e2, in_sim=False))
                E = 0
                # e.p.add_parameter_with_collection(t_R_start(5000, in_sim=False))
            elif e.change_type == "Th":
                e.p.add_parameter_with_collection(t_R_start(1e4, in_sim=False))
                E = 1e4 * N_A ** -1 * 1e9
                # e.p.add_parameter_with_collection(t_R_start(27000, in_sim=False))
            elif e.change_type == "Tnaive":
                e.p.add_parameter_with_collection(t_R_start(1e2, in_sim=False))
                E = 0
                # e.p.add_parameter_with_collection(t_R_start(0, in_sim=False))
            elif e.change_type == "blank":
                e.p.add_parameter_with_collection(t_R_start(0, in_sim=False))
            if E != 0:
                tmp_sigma = np.sqrt(np.log(var ** 2 / E ** 2 + 1))
                mean = np.log(E) - 1 / 2 * tmp_sigma ** 2
                R_draw = np.random.lognormal(mean, tmp_sigma)

                e.p.get_physical_parameter("R", "IL-2").set_in_sim_unit(R_draw)
                e.p.add_parameter_with_collection(t_R_start(R_draw, in_sim=True))


    def set_cell_types(self,sc):
        offset_ids = self.get_offset_ids(sc)
        Tsec_fraction = sc.entity_list[0].p.get_physical_parameter("Tsec_fraction", "IL-2").get_in_post_unit()
        Tsec_draws = self.get_draws(sc, exluded_ids=offset_ids, fraction=Tsec_fraction)

        if len(self.Tsec_distribution_array) != len(Tsec_draws):
            Tsec_distribution_array = Tsec_draws
        else:
            Tsec_draws = self.Tsec_distribution_array

        Th_fraction = sc.entity_list[0].p.get_physical_parameter("Th_fraction", "IL-2").get_in_sim_unit()
        Th_draws = self.get_draws(sc, exluded_ids=list(offset_ids) + list(Tsec_draws), fraction=Th_fraction)

        if len(self.Th_distribution_array) != len(Th_draws):
            Th_distribution_array = Th_draws
        else:
            Th_draws = self.Th_distribution_array

        no_of_Th = 0
        q_il2_sum = len(sc.entity_list) * 0.25 * 10 * N_A ** -1 * 1e9

        for i, e in enumerate(sc.entity_list):
            e.change_type = "Tnaive"
            if i in Tsec_draws:
                e.change_type = "Tsec"
                # e.p.get_physical_parameter("R", "IL-2").set_in_post_unit(5000)
                # print(e.p.get_physical_parameter("R", "IL-2").get_in_post_unit())
            if i in offset_ids:
                e.change_type = "Th"
                no_of_Th += 1
            if i in Th_draws:
                e.change_type = "Th"
                no_of_Th += 1
            e.p.add_parameter_with_collection(MiscParameter("id", int(i)))

        print("Number of secreting cells: ", len(Tsec_draws))
        print("Number of Ths: ", no_of_Th)
        if len(Tsec_draws) != 0:
            print("Cells q:", q_il2_sum / len(Tsec_draws) / (N_A ** -1 * 1e9))
        else:
            print("cells q = 0")
        return Tsec_draws, Th_draws, q_il2_sum


    def step(self,sc):
        Tsec_draws, Th_draws, q_il2_sum = self.set_cell_types(sc)
        self.set_parameters(sc,Tsec_draws, q_il2_sum)
        # print("drawn fraction:", self.drawn_fraction)


    def step_R_lognorm(self,sc):
        Tsec_draws, Th_draws, q_il2_sum = self.set_cell_types(sc)
        self.set_R_lognorm_parameters(sc,Tsec_draws, q_il2_sum)
        # print("drawn fraction:", self.drawn_fraction)