import pandas as pd
from tqdm.auto import tqdm

allele_sum = lambda x, base: int(x.a1 == base) + int(x.a2 == base)

def geno_sum(x, base):
    x = list(x.geno)
    if len(x) == 2:
        if (x[0] == base) and (x[1] == base):
            return 2
        elif (x[0] == base) or (x[1] == base):
            return 1
        else:
            return 0
    elif len(x) == 1:
        if x[0] == base:
            return 1
        else:
            return 0
    else:
        return 0

col_map = lambda col, abbrv: f'{col}_{abbrv}'

bases = ['A', 'C', 'G', 'T']

def load_ancestry(file, abbrv='A', verbose=False, **params):
    loaded = pd.read_csv(file, header=0, names=['rsid', 'chrome', 'position', 'a1', 'a2'], comment='#', sep='\t', **params)
    loaded['name'] = abbrv
    loaded.chrome = loaded.chrome.astype(str)
    for base in tqdm(bases, disable=(not verbose)):
        loaded[base] = loaded.apply(lambda x: allele_sum(x, base), axis=1)
    loaded = loaded.drop(['a1', 'a2'], axis=1)
    return loaded

def load_ancestry2(file, abbrv='A', verbose=False, **params):
    loaded = pd.read_csv(file, header=0, names=['rsid', 'chrome', 'position', 'a1', 'a2'], comment='#', sep='\t', **params)
    loaded['name'] = abbrv
    loaded.chrome = loaded.chrome.astype(str)
    loaded['genotype'] = loaded.a1 + loaded.a2
    loaded = loaded.drop(['a1', 'a2'], axis=1)
    return loaded



def load_common(file, abbrv='C', verbose=False, **params):
    loaded = pd.read_csv(file, header=None, names=['rsid', 'chrome', 'position', 'geno'], comment='#', sep='\t', **params)
    loaded['name'] = abbrv
    loaded.chrome = loaded.chrome.astype(str)
    for base in tqdm(bases, disable=(not verbose)):
        loaded[base] = loaded.apply(lambda x: geno_sum(x, base), axis=1)
    loaded = loaded.drop(['geno'], axis=1)
    return loaded

ancestry = load_ancestry('./AncestryDNA.txt', "Ancestry", True)
twenty3andme = load_common('./23andme.txt', '23&me', True)
living = load_common('./LivingDNA.txt', 'LivingDNA', True)
combined = pd.concat([ancestry, twenty3andme, living])