import pandas as pd
import cell2cell as c2c
import numpy as np
import gseapy
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns
from collections import defaultdict
import os
from tabulate import tabulate

import_path = '/Users/Troks27/Desktop/TrokaChat for Publication - R Scripts/TrokaChat Publication/TrokaChat/TrokaChat Output/'
output_path = '/Users/Troks27/Desktop/TrokaChat for Publication - R Scripts/TrokaChat Publication/TrokaChat/GSEA/'
files_list = ['output_H_allpathways', 'output_NL_allpathways', 'output_LS_allpathways']
sample = ["H","NL", "LS"]
Source = [0,1,2,3,4,5,6,7,8,9,10,11]
Target = [0,1,2,3,4,5,6,7,8,9,10,11]
output_folder = '/Users/Troks27/Desktop/TrokaChat for Publication - R Scripts/TrokaChat Publication/TrokaChat/'
species = "human"
database = '/Users/Troks27/Desktop/9. Human Skin TrokaChat Analysis/Python TrokaChat7/c5.go.v2022.1.Hs.symbols.gmt.txt'



df = pd.DataFrame(columns=['Ligand', 'Receptor', 'Source','Target'])
for i in range(len(files_list)):
    locals()["files" + str(i)] = pd.read_csv(import_path + files_list[i] + '.csv')
    locals()["files" + str(i)] = locals()["files" + str(i)].loc[locals()["files" + str(i)]['Source'].isin(Source)]
    locals()["files" + str(i)] = locals()["files" + str(i)].loc[locals()["files" + str(i)]['Target'].isin(Target)]
    locals()["files" + str(i)] = locals()["files" + str(i)].drop('Uniqueness Score', axis=1).drop('Pathway Label', axis=1)
    locals()["files" + str(i)].rename(columns={"jank_mass_action": sample[i]}, inplace=True)
    locals()["files" + str(i)] = locals()["files" + str(i)].reset_index(drop=True)
    print(locals()["files" + str(i)])
    df = df.merge(locals()["files" + str(i)], on=['Ligand', 'Receptor', 'Source','Target'], how='outer')
df['Receptor'] = df['Receptor'].str.replace('_','&').str.upper()
df['Ligand'] = df['Ligand'].str.upper()
df['df_list'] = df[['Ligand', 'Receptor']].apply(lambda x: '^'.join(x[x.notnull()]), axis = 1)
df = df.drop('Ligand', axis=1).drop('Receptor', axis=1)
for i in range(len(files_list)):
    df[sample[i]] = df.groupby('df_list')[sample[i]].transform('mean')
df = df.drop('Source', axis=1).drop('Target', axis=1)
df = df.drop_duplicates(subset=["df_list"], keep='first').reset_index(drop=True)
first_column = df.pop('df_list')
df.insert(0, 'df_list', first_column)
df = df.drop(df.filter(like='Sample Ident').columns, axis=1)
print("TAKE TBUS SJBAFO")
print(df)

data_folder = './'
directory = os.fsencode(data_folder)


if not os.path.isdir(output_folder):
    os.mkdir(output_folder)



pathway_per_gene = defaultdict(set)
with open(database, 'r') as f:
    print(f)
    for i, line in enumerate(f):
        l = line.split('\t')
        l[-1] = l[-1].replace('\n', '')
        l = [pw for pw in l if ('http' not in pw)] # Remove website info
        for gene in l[1:]:
            pathway_per_gene[gene] = pathway_per_gene[gene].union(set([l[0]]))

print(pathway_per_gene['IL6'])


if species == "human":
    lr_pairs = pd.read_excel('/Users/Troks27/Desktop/8. Leticia All/Python TrokaChat/ALL PATHWAYS.xlsx', sheet_name=0)
    lr_pairs = lr_pairs.astype(str)
elif species == "mouse":
    lr_pairs = pd.read_excel('/Users/Troks27/Desktop/8. Leticia All/Python TrokaChat/human_ALL_PATHWAYS_UPDATED.xlsx', sheet_name=0)
    lr_pairs = lr_pairs.astype(str)


# If the LR db include protein complexes.
# This is the character separating members
complex_sep = '&'

# Dictionary to save the LR interaction (key) and the annotated pathways (values).
pathway_sets = defaultdict(set)

# Iterate through the interactions in the LR DB.
for idx, row in lr_pairs.iterrows():
    lr_label = row['interaction_symbol']
    lr = lr_label.split('^')

    # Gene members of the ligand and the receptor in the LR pair
    if complex_sep is None:
        ligands = [lr[0]]
        receptors = [lr[1]]
    else:
        ligands = lr[0].split(complex_sep)
        receptors = lr[1].split(complex_sep)

    # Find pathways associated with all members of the ligand
    for i, ligand in enumerate(ligands):
        if i == 0:
            ligand_pathways = pathway_per_gene[ligand]
        else:
            ligand_pathways = ligand_pathways.intersection(pathway_per_gene[ligand])

    # Find pathways associated with all members of the receptor
    for i, receptor in enumerate(receptors):
        if i == 0:
            receptor_pathways = pathway_per_gene[receptor]
        else:
            receptor_pathways = receptor_pathways.intersection(pathway_per_gene[receptor])

    # Keep only pathways that are in both ligand and receptor.
    lr_pathways = ligand_pathways.intersection(receptor_pathways)
    for p in lr_pathways:
        pathway_sets[p] = pathway_sets[p].union([lr_label])

K = 15

lr_set = defaultdict(set)

for k, v in pathway_sets.items():
    if len(v) >= K:
        lr_set[k] = v

print(len(lr_set))
#print(lr_set)

weight = 1
min_size = 15
permutations = 2000
significance_threshold = 0.05

print(tabulate(df.head(), headers="keys"))
for factor in df.columns[1:]:
    print("READ ME FACTOR" + factor)
    # Rank LR pairs of each factor by their respective loadings
    test = df[['df_list', factor]]
    test.columns = [0, 1]
    test = test.sort_values(by=1, ascending=False)
    test.reset_index(drop=True, inplace=True)
    test = test.dropna()
    print(test)
    #print(lr_set)


    # RUN GSEA
    pre_res = gseapy.prerank(rnk=test,
                   gene_sets=lr_set,
                   min_size=min_size,
                   weighted_score_type=weight,
                   processes=6,
                   permutation_num=permutations, # reduce number to speed up testing
                   outdir=output_folder + 'GSEA/' + factor, format='png', seed=6)

    print(pre_res.res2d.sort_index())

pvals = []
terms = []
factors = []
nes = []
for factor in df.columns[1:]:
    p_report = pd.read_csv(output_folder + 'GSEA/' + factor + '/gseapy.gene_set.prerank.report.csv')
    pval = p_report['NOM p-val'].values.tolist()
    pvals.extend(pval)
    terms.extend(p_report.Term.values.tolist())
    factors.extend([factor] * len(pval))
    nes.extend(p_report['NES'].values.tolist())
pval_df = pd.DataFrame(np.asarray([factors, terms, nes, pvals]).T, columns=['Factor', 'Term', 'NES', 'P-value'])
pval_df = pval_df.loc[pval_df['P-value'] != 'nan']
pval_df['P-value'] = pd.to_numeric(pval_df['P-value'])
pval_df['P-value'] = pval_df['P-value'].replace(0., 1. / (permutations + 1))
pval_df['NES'] = pd.to_numeric(pval_df['NES'])

pval_df['Adj. P-value'] = fdrcorrection(pval_df['P-value'].values,
                                        alpha=significance_threshold)[1]

print(pval_df.loc[(pval_df['Adj. P-value'] < significance_threshold) & (pval_df['NES'] > 0.)])
print(pval_df.loc[(pval_df['Adj. P-value'] < significance_threshold) & (pval_df['NES'] < 0.)])

pval_df.to_excel(output_folder + 'GSEA-Adj-Pvals.xlsx')

pval_pivot = pval_df.pivot(index="Term", columns="Factor", values="Adj. P-value").fillna(1.)
scores = pval_df.pivot(index="Term", columns="Factor", values="NES").fillna(0)

with sns.axes_style("darkgrid"):
    dotplot = c2c.plotting.pval_plot.generate_dot_plot(pval_df=pval_pivot,
                                                       score_df=scores,
                                                       significance=significance_threshold,
                                                       xlabel='',
                                                       ylabel='KEGG Pathways',
                                                       cbar_title='NES',
                                                       cmap='PuOr',
                                                       figsize=(9,30),
                                                       label_size=20,
                                                       title_size=20,
                                                       tick_size=14,
                                                       filename=output_folder + 'GSEA-Dotplot.png'
                                                      )