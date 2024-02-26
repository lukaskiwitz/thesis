p = {
	"N_cells": 1000, # amount of cells
	"L": 10, # µm, distance between cells surface
	"f_Tsec": 0.05, # cells
	# "N_Tresp": 0.9, # cells
	"rho": 5, # µm
	"d": 20, # µm, cell-to-cell distance
	"m": 20, # µm, margin around volume
	"D":10, # µm^2/s
	"q":10, # molecules/s/cell
	"q_sys": 30*0.1, # systemic q, molecules/s/all cells
	"R":1e2, # molecules/cell
	"k_on": 3.1e7, # 1/(M s)
	"k_off": 23.06e-5, # 1/s
	"k_endo": 0.00046, # 1/s
	"K_D": 45343.81237655793, # molecules. k_off/k_on #care for correct units, k_on in molecules!
	"kd": 0.1/3600, # 1/s
	"R_Tsec": 1e2, # molecules/cell
	"R_Th": 1e4, # molecules/cell
	"R_Treg": 5e3, # molecules/cell
	# feedback
	"feedback": 0, # binary
	"gamma": 100, # unit
	"eta": 1/72000, # 1/s receptor decay
	"nu": 1e-3, # 1/s receptor production
	"c_0": 79000, # molecules
	"N_R": 3, # unit
	"k": 0, # nM, to be set
	"K_m": 79000, #molecules
	"k_min": 1e-3,  # receptors/s
	"k_max": 10, # receptors/s
	"R_mean": 5e3, # molecules/cell
	# saturation
	"a_max": 1, # molecules/(cell s)
	"K_c": 0.010, # nM
	"K_R": 1e2, # molecules/cell
	"K_EC50": 5e4, # molecules/cell
	"EC50_max": 0.1, # nM
	"EC50_min": 0.001, # nM
	"N_EC50": 0.8, # unit
	"N_pSTAT5": 1, # unit
	"sigma": 1, # unit
}