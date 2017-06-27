# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 15:03:59 2016

@author: hliu
"""

import pandas as pd
import os
from StringIO import StringIO

# os.chdir('/home/hliu/Desktop')


def readMolStr(mol_pth):
    with open(mol_pth, 'r') as f:
        mol_str = f.readlines()
    atom_identifier = ['ATOM']
    def isAtomInfo(line):
        results = []
        for i in atom_identifier:
            results.append(line.startswith(i))
        return any(results)
    mol_str = ''.join([a for a in mol_str if isAtomInfo(a)])
    return mol_str


def molStr2DF(mol_str):

    cols = ['Identifier',
            'AtomID',
            'AtomName',
            'ResName',
            'SegName',
            'ResID',
            'X',
            'Y',
            'Z',
            'BK1',
            'BK2',
            'BK3']
    
    mol = pd.read_csv(StringIO(mol_str), delim_whitespace=True, header=None, na_filter=False)
    mol.index = list(range(len(mol.index)))
    mol = mol.dropna(axis=0, how='all')
    mol = mol.dropna(axis=1, how='all')
    mol_cols = cols[:len(mol.columns)]
    mol.columns = mol_cols
    
    def group_by_HexOrTen(index, cutoff):
        if index <= cutoff:
            return 'ten'
        else:
            return 'hex'        
    
    atomNum = len(mol.index)
    resnum = mol.loc[:, 'ResID']
    if atomNum < 10000:
        mol['ResID'] = mol['ResID'].apply(int)
        return mol
    else:
        cutoff = resnum[resnum=='9999'].index[-1]
        grouped = resnum.groupby(lambda x: group_by_HexOrTen(x, cutoff))
        new_resnum = pd.concat([grouped.get_group('ten').apply(int),
                                grouped.get_group('hex').apply(lambda x: int('0x'+x, 16))])
        mol.iloc[:, 4] = new_resnum
    
        return mol

def fetchSeq(mol):
    idx = 1
    resi = mol.ResID.ix[0]
    info = [resi, idx]
    def label_residue(i):
        if i == info[0]:
            return info[1]
        else:
            idx = info[1]+1
            info[0] = i
            info[1] = idx
            return idx
    resi_label = [label_residue(i) for i in mol.ResID]
    resIdx = mol.groupby(resi_label).apply(lambda x: x.index[0])
    seq = mol.loc[resIdx, ['ResID', 'ResName']]
    return seq
