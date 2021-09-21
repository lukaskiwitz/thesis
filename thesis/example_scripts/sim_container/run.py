from parameters import cytokines, cell_types_dict, geometry, numeric, boundary
from thesis.main.SimContainer import SimContainer
from thesis.scenarios.box_grid import setup
from thesis.main.assign_fractions import assign_fractions

path = "/extra/kiwitz/simcontainer_example/test_1/"
sc = setup(
    cytokines,
    cell_types_dict,
    boundary,
    geometry,
    numeric,
    path,
    "")

assign_fractions(sc,0)
sc.run([0,1])
sc.save_fields(1)
sc.save_markers(1)