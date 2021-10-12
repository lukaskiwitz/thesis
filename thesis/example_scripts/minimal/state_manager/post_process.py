import os

from thesis.main.PostProcess import PostProcessor, ParaviewRender

os.environ["PARAVIEW_PATH"] = "/home/kiwitz/ParaView-5.9.0-osmesa"
path = "/extra/kiwitz/statemanager_example/test_1/"
pp = PostProcessor(path)
pp.unit_length_exponent = -6

pp.computations.append(ParaviewRender)
pp.run_post_process(4)
