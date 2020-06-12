<h1>
    Requirements
</h1>

1. All python dependencies (provided as an anaconda yml-file)
2. correctly set LD_LIBRARY_PATH environment variable
3. installed thesisproject-package (provided as wheel file for pip)

<h1>
Installation Instructions (with anaconda)
</h1>

1. Install anaconda:
`https://www.anaconda.com/products/individual`


2. Run

`
conda env create -f fenics_min.yml -n enviroment_name
`

to create a new anaconda environment with the necessary packages.
This could take a moment.

3.Run

`conda activate enviroment_name`

to activate the newly created environment.

4.Install the package using pip

`
pip install wheel-file
`

5.Run

`conda env config vars set LD_LIBRARY_PATH=~/anaconda3/envs/enviroment_name/lib
`

to set up necessary environment variables.

6.Reactivate the conda environment to apply changes

`conda activate enviroment_name`


<h1>
    Usage
</h1>
1.Copy the boxed_static folder to a location of your choice.

2.Open the parameters.py file and edit the path variable at the bottom to point to the desired top level directoy for your simulation results.
Note: If your are on the itb computers this will work well without a change.

3.Run

`python run.py; python post_process.py number_of_threads; python plot.py`

to start the simulation and plot the results.


 

 

