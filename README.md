# thesis project



1. Run

`
cond env create -f fenics_min.yml -n enviroment_name
 `
 
 to create a new anaconda environment with the necessary packages.
 This could take a moment.
 
 2. Run
  
`conda activate enviroment_name`
 
 to activate the newly created environment.
 
 3. Install the package using pip
 
 `
 pip install thesisproject-0.0.1-py3-none-any.whl
 `
 
 4. copy the boxed_static folder to a location of your choice.
 5. open the parameters.py file and edit the path variable a the bottom to change to top level directoy for your simulation result.
 Note: If your are on the itb computers this will work well without a change.
 
6. Run

`python run.py; python post_process.py number_of_threads; python plot.py`

to start the simulation.
 
 
 

