#/bin/sh

# needs https://pypi.org/project/coverage/

python3 -m coverage run -m unittest runall.py

python3 -m coverage html  

