Installation Instructions (with anaconda)
===========================================

Requirements
------------

1. All python dependencies (provided as an anaconda yml-file)
2. correctly set LD_LIBRARY_PATH environment variable
3. installed thesisproject-package


Installation
------------

1. Install anaconda: :code:`https://www.anaconda.com/products/individual`

2. Run :code:`conda env create -f fenics.yml -n "enviroment_name"`
to create a new anaconda environment with the necessary packages.
This could take a moment.

3.Run :code:`conda activate "enviroment_name"`
to activate the newly created environment.

4.Install the package using pip with :code:`pip install wheel-file`
if you got the wheel file or :code:`pip install -e "repository_dir"`
if you cloned this repository.

5.Run :code:`conda env config vars set LD_LIBRARY_PATH=~/anaconda3/envs/enviroment_name/lib`
and if you want to use paraview rendering :code:`conda env config vars set PARAVIEW_PATH=~path_to_paraview_dir`
to set up necessary environment variables.

6.Reactivate the conda environment to apply changes :code:`conda activate enviroment_name`

Paraview Integration
--------------------

Only works for simulation data generated with the current version

1. Download an omesa build of paraview (necessary to run without X-server)

2. Set the PARAVIEW_PATH env. variable to your paraview directory. For instance in you conda env. or in the post_process.py

4. pass the keyword argument "render_paraview=True" to run_post_process()

5. various settings can be passed. "example_min" contains a commented example


The cell type markers are floats and will be set according to "SimContainer.marker_lookup".
SimContainer now has the "markers" attribute, which lists entity
properties (parameters or python attributes) for which to export a MeshFunction.
Assigment of non numerical values is controlled by marker_lookup.