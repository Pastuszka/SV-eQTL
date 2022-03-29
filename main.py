import vcf
import pandas as pd
from scipy.stats import linregress
import numpy as np
import matplotlib.pyplot as plt


class RegressionResult:
    def __init__(self, p_val, variant, gene_symbol, coord):
        self.p_val = p_val
        self.chromosome = variant.CHROM
        self.variant = variant.ID
        self.variant_location = variant.POS
        self.gene_location = coord
        self.gene = gene_symbol

    def __str__(self):
        row = f'{self.chromosome}\t{self.gene}\t{self.variant}\t{self.gene_location}\t{self.variant_location}'
        row += f'\t{self.p_val}'
        return row


def get_n_mutations(gt):
    a = int(gt[0])
    b = int(gt[2])
    if a == 2 or b == 2:
        return None
    else:
        return a + b


def read_variant(variant):
    hets = variant.get_hets()
    variant_dict = {}
    for het in hets:
        gt = het.data.GT
        gt_val = get_n_mutations(gt)
        if gt_val:
            variant_dict[het.sample] = gt_val
    return variant_dict


def make_model(variant, row, specimen):
    variant_vals = read_variant(variant)
    x = np.empty(len(specimen))
    y = np.empty(len(specimen))
    keys = variant_vals.keys()
    for i, sample in enumerate(specimen):
        if sample in keys:
            y[i] = variant_vals[sample]
        else:
            y[i] = 0
        x[i] = row['specimen']
    return linregress(x, y)


def eqtl_regression(row, sv_reader, p_val, specimen, max_distance):
    chromosome = row['Chr']
    gene_symbol = row['Gene_Symbol']
    coord = row['Coord']
    variants = sv_reader.fetch(chromosome, coord-max_distance, coord+max_distance)
    regression_results = []
    for variant in variants:
        model = make_model(variant, row, specimen)
        if model.pvalue < p_val:
            result = RegressionResult(model.pvalue, variant, gene_symbol, coord)
            regression_results.append(result)
    return regression_results


def eqtl_analysis(sv_path, expression_path, max_distance, p_val):
    sv_reader = vcf.Reader(open(sv_path, 'r'))
    expressions = pd.read_table(expression_path)
    specimen = list(expressions.columns[4:])
    found = []
    for index, row in expressions.iterrows():
        regression_results = eqtl_regression(row, sv_reader, p_val, specimen, max_distance)
        found.extend(regression_results)
    return found


def main():
    sv_path = '../data/go_files/eqtl_anal/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf'
    expression_path = '../data/go_files/eqtl_anal/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt'
    max_distance = 1000
    p_val = 0.05
    results = eqtl_analysis(sv_path, expression_path, max_distance, p_val)


if __name__ == '__main__':
    main()

