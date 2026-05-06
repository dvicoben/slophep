# Set-up

SLOP is python-based. There is currently no `pyPI` installation, but SLOP can be easily installed after cloning the repository.

## Using pip
In the python environment of your choice, 

```
git clone https://github.com/dvicoben/slophep.git
cd slophep
pip install -e .
```

## Quick
Ensure you are in a python environment with all the requirements in [`requirements.txt`](https://github.com/dvicoben/slophep/blob/master/requirements.txt)
```
git clone https://github.com/dvicoben/slophep.git
cd slophep
source ./setup.sh
```
The script setup.sh simply appends `src/` to the `PYTHONPATH` so that contents therein will be found when running scripts. You will need to `source ./setup.sh` whenever you start a new terminal session.