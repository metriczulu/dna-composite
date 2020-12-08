import pandas as pd

class DNA:

    def __init__(self, dna_file, gender='male'):
        if type(dna_file) == str:
            self.dna_file = [dna_file]
        else:
            self.dna_file = dna_file
        self.total_types = list()
        self.gender = gender
        self.valid = {'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT',
                      'A', 'C', 'G', 'T'}
        self.dna = self._load_dna(self.dna_file[0])
        if len(self.dna_file) > 1:
            for additional_dna in self.dna_file[1:]:
                self.dna = self._stack(self.dna, self._load_dna(additional_dna))

    def _find_first(self, series):
        for row in series:
            if row in self.valid:
                return row
        return "--"

    def _load_dna(self, dna_file):
        file_type = self._detect_type(dna_file)
        if file_type == 'Ancestry':
            return self._load_ancestry(dna_file)
        else:
            return self._load_common(dna_file, abbrv=file_type)

    def _detect_type(self, dna_file):
        with open(dna_file, 'r') as files:
            first_line = files.readline()
            if 'AncestryDNA' in first_line:
                self.total_types.append('Ancestry')
                return 'Ancestry'
            elif '23andMe' in first_line:
                self.total_types.append('23andMe')
                return '23andMe'
            elif 'Living DNA' in first_line:
                self.total_types.append('Living DNA')
                return 'LivingDNA'
            else:
                self.total_types.append('Common')
                return 'Common'

    def _remove_missing(self, x):
        if x in self.valid:
            return x
        else:
            return '--'

    def _load_ancestry(self, file, abbrv='Ancestry', **params):
        df = pd.read_csv(file, header=0, names=['rsid', 'chromosome', 'position', 'a1', 'a2'], comment='#', sep='\t', **params)
        df['name'] = abbrv
        df.chromosome = df.chromosome.astype(str).replace({'23': 'X', '25': 'X', '24': 'Y', '26': 'MT'})
        if self.gender == 'male':
            df.a2[df.chromosome == 'X'] = ""
            df.a2[df.chromosome == 'Y'] = ""
        df.a2[df.chromosome == 'MT'] = ""
        df['genotype'] = df.a1 + df.a2
        df.genotype = df.genotype.apply(self._remove_missing)
        df = df.drop(['a1', 'a2'], axis=1)
        return df

    def _load_common(self, file, abbrv='Common', **params):
        df = pd.read_csv(file, header=None, names=['rsid', 'chromosome', 'position', 'genotype'], comment='#', sep='\t', **params)
        df['name'] = abbrv
        df.chromosome = df.chromosome.astype(str)
        df.genotype = df.genotype.apply(self._remove_missing)
        return df

    def _stack(self, dna, new_dna):
        combined = pd.concat([dna, new_dna], sort=False)
        grouped = combined.groupby('rsid').agg({'rsid': 'first', 'chromosome': 'first', 'position': 'first',
                                             'name': "/".join, 'genotype': self._find_first})
        return grouped.sort_values(['chromosome', 'position'])

    def to_csv(self, file_path="/DNA.txt", **params):
        with open(file_path, 'w') as f:
            f.write(f"# Combined DNA files - {len(self.total_types)}; 23andMe format\n")
            f.write("#\n")
            f.write("# Raw DNA files combined using program written by Shane Stephenson - www.github.com/metriczulu")
            f.write("#\n")
            f.write("#\n")
        self.dna.drop(['name'], axis=1).to_csv(file_path, index=False, sep='\t', mode="a", **params)