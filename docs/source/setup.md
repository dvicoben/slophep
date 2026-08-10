# Set-up

SLOP is python-based. There is currently no `pyPI` installation, but SLOP can be easily installed after cloning the repository.


## Requirements
Requirements are listed in [`requirements.txt`](https://github.com/dvicoben/slophep/blob/master/requirements.txt)
- Predictions use `flavio` to go from FFs and WCs to amplitudes and observables
- Internally also uses common libraries (`numpy`, `matplotlib`, `scipy`), and `pypmc` for MCMC sampling

## Quick
Ensure you are in a python environment with all the requirements in [`requirements.txt`](https://github.com/dvicoben/slophep/blob/master/requirements.txt)
```
git clone https://github.com/dvicoben/slophep.git
cd slophep
```
Then append the `src` directory to the `PYTHONPATH`
```
export PYTHONPATH="${PYTHONPATH}:/path/to/slophep/src"
```
such that contents therein will be found when running scripts.

The script `setup.sh` performs this assuming it is run from the root directory of the repository. You will need to `source ./setup.sh` whenever you start a new terminal session.

## Using pip
There is no `PyPI` release yet, but you can use `pip` to install the package from its source. In the python environment of your choice, 
```
git clone https://github.com/dvicoben/slophep.git
cd slophep
pip install -e .
```
which should install the package (`slophep`) and the required dependencies.