### Installing the right packages
`conda env create -f environment.yml`

### Activate the environment
`conda activate py2`

### Running biopy
Go to biopy directory and run the following:\

```cd biopy```\
`python make_consensus.py -i inp_file -o output_file`\
inp_file is the `clusters_bipart_splits.txt' with the path 

### Setting up biopy2
This is a standard python package. To install run 'python setup.py install'.
To install locally run 'python setup.py install --home=DIR', and add 'DIR/lib/python'
to the PYTHONPATH environment variable.

For example:
  python setup.py install  --home=./`python --version 2>&1 | sed -s "s/Python /version-/"`
