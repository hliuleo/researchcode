# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 15:03:59 2016

@author: hliu
"""

import pandas as pd
import os
from StringIO import StringIO

os.chdir('/home/hliu/Desktop')


def readMol(mol_pth):

    with open(mol_pth, 'r') as f:
        mol_str = f.readlines()
    mol_str = ''.join([a.replace('X', ' ') for a in mol_str])
    
    mol = pd.read_csv(StringIO(mol_str), delim_whitespace=True, header=None, skiprows=1, na_filter=False)
    mol.index = list(range(len(mol.index)))
    mol = mol.dropna(axis=0, how='all')
    mol = mol.dropna(axis=1, how='all')
    mol = mol.drop(mol.index[-1])
    
    def group_by_HexOrTen(index, cutoff):
        if index <= cutoff:
            return 'ten'
        else:
            return 'hex'        
    
    atomNum = len(mol.index)
    resnum = mol.iloc[:, 4]
    if atomNum < 10000:
        mol[4] = mol[4].apply(int)
        return mol
    else:
        cutoff = resnum[resnum=='9999'].index[-1]
        grouped = resnum.groupby(lambda x: group_by_HexOrTen(x, cutoff))
        new_resnum = pd.concat([grouped.get_group('ten').apply(int),
                                grouped.get_group('hex').apply(lambda x: int('0x'+x, 16))])
        mol.iloc[:, 4] = new_resnum
    
        return mol


def fetchSeq(mol):
    resIdx = mol.groupby(4).apply(lambda x: x.index[0])
    seq = mol.iloc[resIdx, [3,4]]
    return seq
