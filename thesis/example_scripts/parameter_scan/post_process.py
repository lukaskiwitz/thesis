from thesis.main.PostProcess import PostProcessor

path = "/extra/kiwitz/statemanager_example/test_1/"
pp = PostProcessor(path)
pp.unit_length_exponent = -6
pp.run_post_process(4)
