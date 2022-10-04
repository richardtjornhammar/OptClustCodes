import pandas as pd

if __name__ == '__main__' :
    print ( 'HERE' )
    import os
    os.system('wget http://ftp.theochem.ruhr-uni-bochum.de/outgoing/praktikum/cpmd-300k/coord.pdb')
    os.system('wget http://www.theochem.ruhr-uni-bochum.de/~legacy.akohlmey/files/32spce-h3op-1ns.xyz')
    os.system('mv 32spce-h3op-1ns.xyz ../data/.')

