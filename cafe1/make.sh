#! /bin/tcsh
# make all the figures 

echo " set environmnent"
conda activate myenv
conda env list

nb2py.sh bias.ipynb
python < bias.py
exit

nb2py.sh coh.ipynb
python < coh.py

nb2py.sh ENSO.ipynb
python < ENSO.py

nb2py.sh PDO.ipynb
python < PDO.py


