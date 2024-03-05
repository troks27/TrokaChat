import pandas as pd
from tabulate import tabulate
import statistics
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from labellines import labelLine, labelLines
from plotapi import Chord
import plotapi
plotapi.api_key("e04a1729-af91-4b04-8df0-bca9f42a2133")
import statistics
import itertools
import logging
from collections import Counter
from functools import reduce
import matplotlib.cm as cm

mpl.use('macosx')
conditions_list = ['NG','DIAB']
colors_list = ["#00CD6C", "#009ADE"]
conditions_to_test = 2
counts_filepath = "/Volumes/LaCie/Bushra_study_2.28.24/TrokaChat/DEG Output/"
import_filepath = '/Volumes/LaCie/Bushra_study_2.28.24/TrokaChat/TrokaChat Compare/'
nulldist_filepath = '/Volumes/LaCie/Bushra_study_2.28.24/TrokaChat/TrokaChat Output/'
output_filepath = '/Volumes/LaCie/Bushra_study_2.28.24/TrokaChat/TrokaChat Charts/'
twoway_filenames = ["DIAB_NG_output_int"]
threeway_filename = []
pval_cutoff = 0.05
condition_of_interest = "NG"
signaling_group = ["Secreted Signaling"]
ligand = []
receptor = []
source = []
target = []
fig_labels = "No"
top = 100000
names = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
#Either "LR_name" or "Source to Target"
label_source = "LR_name"
min_num_cells = 30









def movecol(df, cols_to_move=[], ref_col='', place='After'):
    cols = df.columns.tolist()
    if place == 'After':
        seg1 = cols[:list(cols).index(ref_col) + 1]
        seg2 = cols_to_move
    if place == 'Before':
        seg1 = cols[:list(cols).index(ref_col)]
        seg2 = cols_to_move + [ref_col]
    seg1 = [i for i in seg1 if i not in seg2]
    seg3 = [i for i in cols if i not in seg1 + seg2]
    return (df[seg1 + seg2 + seg3])



conditions_list_clean = conditions_list
def polymul(p, q):
    """
    Multiply two polynomials, represented as lists of coefficients.
    """
    r = [0]*(len(p) + len(q) - 1)
    for i, c in enumerate(p):
        for j, d in enumerate(q):
            r[i+j] += c*d
    return r
def ncombinations(it, k):
    """
    Number of combinations of length *k* of the elements of *it*.
    """
    counts = Counter(it).values()
    prod = reduce(polymul, [[1]*(count+1) for count in counts], [1])
    return prod[k] if k < len(prod) else 0
number_of_condition_combinations = ncombinations(conditions_list, 2)
print(number_of_condition_combinations)



if len(conditions_list) == 3:
    num_3_way = 1
else:
    num_3_way = 0



d = {}
for i in range(num_3_way):
    data = pd.read_csv(import_filepath + threeway_filename[i] + ".csv")
    data.dropna(how='all', axis=1, inplace=True)
    data.drop(data.columns[[16, 17, 18]], axis=1, inplace=True)
    name1 = data.iloc[0,7]
    name2 = data.iloc[0,11]
    name3 = data.iloc[0, 15]
    data.rename(columns={data.columns[4]: 'Communication Score_' + name1, data.columns[5]: 'Uniqueness Score_' + name1, data.columns[6]: 'Pathway Label_' + name1, data.columns[7]: 'Sample Ident_'+ name1}, inplace=True)
    data.rename(columns={data.columns[8]: 'Communication Score_' + name2, data.columns[9]: 'Uniqueness Score_' + name2, data.columns[10]: 'Pathway Label_' + name2, data.columns[11]: 'Sample Ident_'+ name2}, inplace=True)
    data.rename(columns={data.columns[12]: 'Communication Score_' + name3, data.columns[13]: 'Uniqueness Score_' + name3, data.columns[14]: 'Pathway Label_' + name3, data.columns[15]: 'Sample Ident_'+ name3}, inplace=True)
    name4 = threeway_filename[i].replace('_output_int', '')
    d[name4] = data



for i in range(number_of_condition_combinations):
    data = pd.read_csv(import_filepath + twoway_filenames[i] + ".csv")
    data.dropna(how='all', axis=1, inplace=True)
    data.drop(data.columns[[12, 13]], axis=1, inplace=True)
    name1 = data.iloc[0,7]
    name2 = data.iloc[0,11]
    data.rename(columns={data.columns[4]: 'Communication Score_' + name1, data.columns[5]: 'Uniqueness Score_' + name1, data.columns[6]: 'Pathway Label_' + name1, data.columns[7]: 'Sample Ident_'+ name1}, inplace=True)
    data.rename(columns={data.columns[8]: 'Communication Score_' + name2, data.columns[9]: 'Uniqueness Score_' + name2, data.columns[10]: 'Pathway Label_' + name2, data.columns[11]: 'Sample Ident_'+ name2}, inplace=True)
    name3 = twoway_filenames[i].replace('_output_int', '')
    d[name3] = data



for i in range(len(conditions_list)):
    data = pd.read_csv(import_filepath + conditions_list[i] + "_outputfinal.csv")
    data.dropna(how='all', axis=1, inplace=True)
    data.rename(columns={data.columns[4]: 'Communication Score_' + conditions_list[i], data.columns[5]: 'Uniqueness Score_' + conditions_list[i], data.columns[6]: 'Pathway Label_' + conditions_list[i], data.columns[7]: 'Sample Ident_'+ conditions_list[i]}, inplace=True)
    d[conditions_list[i]] = data
combined = pd.concat(d.values(), ignore_index=True)
print(tabulate(combined.head(), headers="keys"))



name_dict = {}
for i in range(len(conditions_list)):
    name_dict["name"+str(i)] = combined.iloc[0, 7+(4*i)]
    combined.iloc[:, 7+(4*i)] = name_dict["name"+str(i)]
print(combined.iloc[:, 7].head())
print(combined.iloc[:, 11].head())



for i in range(len(conditions_list)):
    for a in range(len(conditions_list)):
        combined["new"] = combined.iloc[:, 6+(4*i)].fillna(combined.iloc[:, 6+(4*i)])
    combined.iloc[:,6+(4*i)] = combined["new"]
combined.drop('new', axis=1, inplace=True)
combined.fillna(0, inplace=True)



namedict = {}
for i in range(len(conditions_list)):
    name0 = combined.iloc[0, 7+(4*i)]
    namedict[i] = name0
print(namedict)



for i in range(len(conditions_list)):
    data = pd.read_excel(nulldist_filepath + "NULL_DIST_" + namedict[i] + " DEGs_allpathways.xlsx")
    stdev = statistics.stdev(data.iloc[:,4])
    mean = statistics.mean(data.iloc[:,4])
    combined["zscore_"+namedict[i]] = combined["Communication Score_"+namedict[i]].sub(mean).div(stdev)
    #combined = combined.loc[combined["Communication Score_" + namedict[i]] > 0]
    combined["pval_"+namedict[i]] = scipy.stats.norm.sf(combined["zscore_"+namedict[i]])



# Creates all of the combinations of conditions based on number. For example, if you have three conditions, there will be 3 combinations of numbers.
l = [i for i in range(len(conditions_list))]
conditions_numbers_combinations = list(itertools.combinations(l, 2))
print(conditions_numbers_combinations)



# Creates all of the combinations of columns containing the Communication Score for each condtion. For example, if you have three conditions, there will be 3 combinations of communication scores.
l2 = [(i*4) + 4 for i in range(len(conditions_list))]
columns_numbers_combinations = list(itertools.combinations(l2, 2))
print(columns_numbers_combinations)



# Does the condition-pairwise divisions of Communication Score.
for i in range(number_of_condition_combinations):
    combined[namedict[conditions_numbers_combinations[i][0]] + "/" + namedict[conditions_numbers_combinations[i][1]]] = combined.iloc[:, columns_numbers_combinations[i][0]].div(combined.iloc[:, columns_numbers_combinations[i][1]])



newdf = pd.DataFrame()
newdf['avg'] = combined.iloc[:, l2].mean(axis=1)
newdf["stdev"] = combined.iloc[:, l2].std(axis=1)
combined["rel_stdev"] = ((newdf["stdev"].div(newdf['avg']))*100).round(2)
combined = combined[(combined[['pval_'+name0,'pval_'+name1,'pval_'+name2]] < pval_cutoff).any(axis=1)]
combined = combined[combined['Pathway Label_'+name0].isin(signaling_group)]



if len(receptor):
    combined = combined[combined["Receptor"].isin(receptor)]
else:
    logging.debug("null")

if len(ligand):
    combined = combined[combined["Ligand"].isin(ligand)]
else:
    logging.debug("null")

if len(target):
    combined = combined[combined["Target"].isin(target)]
else:
    logging.debug("null")

if len(source):
    combined = combined[combined["Source"].isin(source)]
else:
    logging.debug("null")



cols = ["Source", "Target"]
combined["Source to Target"] = combined[cols].apply(lambda row: '-->'.join(row.values.astype(str)), axis=1)
combined.to_csv(output_filepath + 'THISISIT.csv', index=False)



if conditions_to_test == 1:
    counts = pd.read_csv('/Users/Troks27/Desktop/9. Human Skin TrokaChat Analysis/Python TrokaChat7/DEG Output/counts' + '.csv')
    counts = counts.fillna(0)
    counts = counts.loc[counts['sample'] == condition_of_interest].reset_index()
    counts = counts.drop('sample', axis=1)
    counts.apply(pd.to_numeric)
    row = counts.apply(lambda row: row[row > 25], axis=1).reset_index(drop=True)
    cols = row.columns
    cols2 = cols.tolist()
    cols3 = [int(i) for i in cols2]
    row2 = list(set(names) - set(cols3))
    logging.debug(row2)
    combined.loc[combined['Source'].isin(row2), 'Communication Score_' + condition_of_interest] = 0
    combined.loc[combined['Target'].isin(row2), 'Communication Score_' + condition_of_interest] = 0



    combined = combined.loc[(combined["pval_" + condition_of_interest] < pval_cutoff)]
    combined = combined.loc[(combined['Communication Score_' + condition_of_interest] > 0)]



    logging.debug(tabulate(combined.head(), headers="keys"))
    combined["LR_name"] = combined.Ligand.str.cat(combined.Receptor, sep='^')
    cols = ["Source", "Target"]
    combined["Source to Target"] = combined[cols].apply(lambda row: '-->'.join(row.values.astype(str)), axis=1)
    combined = combined.head(top)
    combined.to_csv(output_filepath + 'TESTnew.csv', index=False)



    #Creating the dot plot by condition.
    combined_new = combined.loc[(combined['Communication Score_' + condition_of_interest] > 0)]
    labels = np.array(combined_new["LR_name"])
    labels1 = np.array(combined_new["Source to Target"])
    labels_dict = {"LR_name": labels, "Source to Target": labels1}



    plt.figure(figsize=(8, 6))
    plt.scatter(combined_new.iloc[:, 12],
                combined_new.iloc[:, 13],
                s=30,
                c="green",
                label=condition_of_interest)
    plt.legend(loc="upper left")
    plt.yscale("log")
    plt.xlabel("Communication Score", size=20)
    plt.ylabel("Uniqueness Score", size=20)
    plt.title("Ligand-Receptor Interactions in Condition "+condition_of_interest, size=24, fontname="Helvetica Neue", fontweight="bold")
    plt.savefig(output_filepath + condition_of_interest + '.pdf', bbox_inches="tight", pad_inches=0.5, transparent=True, dpi=300)
    plt.show()



    #First bar plot
    newdict = {}
    for i in range(len(condition_of_interest)):
        combined2 = combined.loc[combined["pval_" + condition_of_interest] < pval_cutoff]
        combined2 = combined2.drop(combined2.iloc[:, 2:28], axis=1)
        combined2 = combined2.reset_index()
        newdict[i] = combined2.groupby(['Ligand', 'Receptor']).count()
    df = newdict[0].rename(columns={'index': condition_of_interest})
    df = df.fillna(0)
    df = df.reset_index()
    df["LR_name"] = df.Ligand.str.cat(df.Receptor, sep='^')
    df['avg'] = df[[condition_of_interest]].mean(axis=1)
    df["stdev"] = df[[condition_of_interest]].std(axis=1)
    df["rel_stdev"] = ((df["stdev"].div(df['avg'])) * 100).round(2)
    df["rel_stdev and rel_stdev_Rank"] = '(σ=' + df['rel_stdev'].astype(str) + ')'
    df = df.sort_values('avg')
    logging.debug(tabulate(df, headers="keys"))



    fig = plt.subplots(figsize=(18, 9))
    barWidth = 0.25
    H = df[namedict[2]]
    br1 = np.arange(len(H))
    br2 = [x + barWidth for x in br1]
    br3 = [x + barWidth for x in br2]



    # Horizontal Bar Plot
    plt.barh(br1, H, label=condition_of_interest, height=barWidth, color="green")
    #plt.xlabel('Number of LR Interactions (Including All Source --> Target Interactions)', fontweight='bold', fontsize=15)
    #plt.ylabel('LR Interactions', fontweight='bold', fontsize=15)
    plt.yticks([r + barWidth for r in range(len(H))], df['LR_name'])
    plt.legend()
    plt.savefig(output_filepath + condition_of_interest+ ' first chart.pdf',bbox_inches="tight", pad_inches=0.5, transparent=True, dpi=300)
    plt.show()


    #This is the flow chart
    newdict2 = {}
    # Tabulating the communication network
    for i in range(len(condition_of_interest)):
        combined2 = combined.loc[combined["pval_" + condition_of_interest] < pval_cutoff]
        combined2 = combined2.iloc[:, 2:4]
        combined2 = combined2.reset_index()
        newdict2[i] = combined2.groupby(['Source', 'Target']).count()
    df2 = newdict2[0].rename(columns={'index': condition_of_interest})
    df2 = df2.fillna(0)
    df2 = df2.reset_index()
    source = df2['Source'].tolist()
    target = df2['Target'].tolist()
    names = list(set(target + source))
    main_list_source = list(set(names) - set(source))
    main_list_target = list(set(names) - set(target))
    xtra = {'Source': main_list_source}
    xtra2 = {'Target': main_list_target}
    df2 = df2.append(pd.DataFrame(xtra))
    df2 = df2.append(pd.DataFrame(xtra2))
    df2 = df2.fillna(0)
    logging.debug(tabulate(df2, headers="keys"))

    combined2 = movecol(combined,
                   cols_to_move=['Source', 'Target'],
                   ref_col='Ligand',
                   place='Before').round(decimals=3)
    combined2 = combined[['Source', 'Target','Ligand','Receptor', 'Communication Score_' + condition_of_interest,
                          'Uniqueness Score_' + condition_of_interest, 'Pathway Label_' + condition_of_interest,
                          'pval_' + condition_of_interest, 'rel_stdev', 'Source to Target']]



    for i in range(len(condition_of_interest)):
        pivot = pd.pivot_table(df2, values=condition_of_interest,
                               index='Source',
                               columns='Target')
        logging.debug(pivot)
        pivot = pivot.fillna(0)
        pivot = pivot.to_numpy()
        pivot = np.round(pivot).astype(int)
        pivot = pivot.tolist()
        logging.debug(pivot)
        names = names
        Chord(pivot, names, data_table=combined2.to_csv(index=False), directed=True,
              title="Source to Target Communication in: " + condition_of_interest).to_pdf(
            output_filepath + "plot_" + condition_of_interest)



    #Next Bar Chart
    newdict3 = {}
    for i in range(len(condition_of_interest)):
        combined2 = combined.loc[combined["pval_" + condition_of_interest] < pval_cutoff]
        combined2 = combined2.iloc[:, [0, 1, (12)]]
        logging.debug(combined2)
        combined2 = combined2.reset_index()
        newdict3[i] = combined2.groupby(['Ligand', 'Receptor']).mean()
    df = newdict3[0].rename(columns={'index': namedict[0]})
    df = df.fillna(0)
    df = df.reset_index()
    df["LR_name"] = df.Ligand.str.cat(df.Receptor, sep='^')
    df['avg'] = df[['Communication Score_' + condition_of_interest]].mean(axis=1)
    df = df.sort_values('avg')
    logging.debug(tabulate(df, headers="keys"))



    fig = plt.subplots(figsize=(18, 9))
    barWidth = 0.25
    H = df['Communication Score_' + condition_of_interest]
    br1 = np.arange(len(H))



    # Horizontal Bar Plot
    plt.barh(br1, H, label=condition_of_interest, height=barWidth, color ="green")
   # plt.xlabel('Average Communication Score (Including All Source --> Target Interactions)', fontweight='bold', fontsize=15)
   # plt.ylabel('LR Interactions', fontweight='bold', fontsize=15)
    plt.yticks([r + barWidth for r in range(len(H))], df['LR_name'])
    plt.legend()
    plt.savefig(output_filepath + condition_of_interest + ' second chart.pdf',
        bbox_inches="tight", pad_inches=0.5, transparent=True, dpi=300)
    plt.show()



# Takes the list of conditions and reverses the order.
conditions_list.reverse()
list_combinations = list(itertools.combinations(conditions_list, 2))
print(list_combinations)



# Takes the the reversed list and finds all combinations of the items without duplicates.
joined_combinations = [' vs. '.join(str_list) for str_list in list_combinations]
print(joined_combinations)



# Takes the list of colors and makes all of the combinations without duplicates.
colors_list.reverse()
color_combinations = list(itertools.combinations(colors_list, 2))
print(color_combinations)



l.reverse()
conditions_numbers_combinations = list(itertools.combinations(l, 2))
print(conditions_numbers_combinations)



combined_dict={}
combined_dict2={}
conditions_list = conditions_list_clean
for i in range(len(conditions_list)):
    #combined.to_csv(
    #    '/Users/Troks27/Desktop/9. Human Skin TrokaChat Analysis/Python TrokaChat7/TrokaChat Charts/For Figures/save_3way/combined'+str(i) + '.csv',
    #    index=False)
    #combined
    counts = pd.read_csv(counts_filepath + 'counts'+'.csv')
    counts = counts.fillna(0)
    counts = counts.loc[counts['sample'] == conditions_list[i]].reset_index()
    counts = counts.drop('sample', axis=1)
    counts.apply(pd.to_numeric)
    counts_new = counts
    row = counts_new.apply(lambda row: row[row > min_num_cells], axis=1).reset_index(drop=True)
    cols = row.columns
    cols2 = cols.tolist()
    cols3 = [int(s) for s in cols2]
    row2 = list(set(names) - set(cols3))
    logging.debug(row2)
    logging.debug(conditions_list[i])
    logging.debug(combined[['Communication Score_'+conditions_list[i], "Target"]])
    combined.loc[combined['Source'].isin(row2), ['Communication Score_'+conditions_list[i]]] = 0
    combined.loc[combined['Target'].isin(row2), ['Communication Score_' + conditions_list[i]]] = 0

    combined.loc[combined['Source'].isin(row2), ['Uniqueness Score_'+conditions_list[i]]] = 0
    combined.loc[combined['Target'].isin(row2), ['Uniqueness Score_' + conditions_list[i]]] = 0

    combined.loc[combined['Source'].isin(row2), ['zscore_'+conditions_list[i]]] = 0
    combined.loc[combined['Target'].isin(row2), ['zscore_' + conditions_list[i]]] = 0

    combined.loc[combined['Source'].isin(row2), ['pval_'+conditions_list[i]]] = 0.99999
    combined.loc[combined['Target'].isin(row2), ['pval_' + conditions_list[i]]] = 0.99999

for i in range(len(conditions_list)):
    combined_dict[i] = combined
    combined_dict2[i] = combined
    combined_dict[i].to_csv(
        output_filepath + 'combined_dict'+str(i) + 'overall_filtered.csv',
        index=False)



logging.debug(combined_dict)


combined_dict3 = combined_dict
combined_dict4 = combined_dict

only_condition1 = {}
only_condition2 = {}
if conditions_to_test == 2:
    for i in range(number_of_condition_combinations):
        trust = combined_dict[i].copy(deep=True)
        trust = trust.loc[(trust['Communication Score_' + list_combinations[i][0]] > 0) | (trust['Communication Score_' + list_combinations[i][1]] > 0)]
        trust = trust.loc[(trust["pval_" + list_combinations[i][0]] < pval_cutoff) | (trust["pval_" + list_combinations[i][1]] < pval_cutoff)]
        trust = trust.loc[(trust["zscore_" + list_combinations[i][0]] > 0) | (trust["zscore_" + list_combinations[i][1]] > 0)]
        trust.sort_values("rel_stdev", ascending=False, inplace=True)

        trust["LR_name"] = trust.Ligand.str.cat(trust.Receptor, sep='^')
        cols = ["Source", "Target"]
        trust["Source to Target"] = trust[cols].apply(lambda row: '-->'.join(row.values.astype(str)), axis=1)
        trust = trust.head(top)
        logging.debug(trust)
        trust.to_csv(output_filepath + 'LOOK AT THIS_'+joined_combinations[i]+'.csv', index=False)





        #This is the creation of the dot plots with those pathways shared in both condtions
        combined_new = trust.loc[(trust['Communication Score_'+ list_combinations[i][0]] > 0) & (trust['Communication Score_'+ list_combinations[i][1]] > 0)]
        X_coords = np.array([combined_new['Communication Score_' + list_combinations[i][0]], combined_new['Communication Score_' + list_combinations[i][1]]])
        Y_coords = np.array([combined_new['Uniqueness Score_' + list_combinations[i][0]], combined_new['Uniqueness Score_' + list_combinations[i][1]]])
        labels = np.array(combined_new["LR_name"])
        labels1 = np.array(combined_new["Source to Target"])
        labels_dict = {"LR_name": labels, "Source to Target": labels1}

        plt.figure(figsize=(8, 6))
        plt.scatter(combined_new['Communication Score_' + list_combinations[i][0]],
                    combined_new['Uniqueness Score_' + list_combinations[i][0]],
                    s=30,
                    c=color_combinations[i][0],
                    label=list_combinations[i][0])

        plt.scatter(combined_new['Communication Score_' + list_combinations[i][1]],
                    combined_new['Uniqueness Score_' + list_combinations[i][1]],
                    s=30,
                    c=color_combinations[i][1],
                    label=list_combinations[i][1])
        plt.legend(loc="upper right")
        plt.plot(X_coords,
                 Y_coords,
                 color='gray',
                 linewidth=.25,
                 label=labels_dict[label_source])

        if fig_labels == 'Yes':
            labelLines(plt.gca().get_lines(), zorder=2.5, align=True, fontsize=5)
        else:
            logging.debug("Bruh")



        plt.xlim([0, 20])
        plt.ylim([.05, 100000])
        plt.yscale("log")

        #plt.xlabel("Communication Score", size=14, fontweight = "bold")
        #plt.ylabel("Uniqueness Score", size=14, fontweight = "bold")
        #plt.title(joined_combinations[i], size=18, fontweight = "bold")
        plt.savefig(output_filepath + joined_combinations[i] + '.pdf', bbox_inches="tight", pad_inches=0.5, transparent=True, dpi=300)
        combined_new.to_csv(output_filepath+joined_combinations[i]+'.csv', index=False)

        X_coords_2 = np.array([[combined_new['Communication Score_' + list_combinations[i][0]]],[combined_new['Uniqueness Score_' + list_combinations[i][0]]]])
        Y_coords_2 = np.array([[combined_new['Communication Score_' + list_combinations[i][1]]], [combined_new['Uniqueness Score_' + list_combinations[i][1]]]])

        distances = np.square((Y_coords_2 - X_coords_2))
        distances = np.sqrt(distances)
        avgs = np.average(distances)
        logging.debug(joined_combinations[i]+" Comparrison Similarity "+ str(avgs))








        #This is the creation of the dot plots with those pathways one in one or the other condition
        only_condition1[i] = combined_dict2[i]
        only_condition2[i] = combined_dict2[i]
        only_condition1[i] = only_condition1[i].loc[(only_condition1[i]['Communication Score_' + list_combinations[i][0]] > 0) & (only_condition1[i]['Communication Score_'+ list_combinations[i][1]] == 0) & (only_condition1[i]['pval_'+ list_combinations[i][0]] < pval_cutoff)]
        only_condition2[i] = only_condition2[i].loc[(only_condition2[i]['Communication Score_' + list_combinations[i][0]] == 0) & (only_condition2[i]['Communication Score_'+ list_combinations[i][1]] > 0) & (only_condition2[i]['pval_' + list_combinations[i][1]] < pval_cutoff)]


        labels = np.array(combined_new["LR_name"])
        labels1 = np.array(combined_new["Source to Target"])
        labels_dict = {"LR_name": labels, "Source to Target": labels1}

        plt.figure(figsize=(8, 6))
        plt.scatter(only_condition1[i]['Communication Score_' + list_combinations[i][0]],
                    only_condition1[i]['Uniqueness Score_' + list_combinations[i][0]],
                    s=30,
                    c=color_combinations[i][0],
                    label=list_combinations[i][0])

        plt.scatter(only_condition2[i]['Communication Score_' + list_combinations[i][1]],
                    only_condition2[i]['Uniqueness Score_' + list_combinations[i][1]],
                    s=30,
                    c=color_combinations[i][1],
                    label=list_combinations[i][1])
        plt.legend(loc="upper right")

        plt.xlim([0, 20])
        plt.ylim([.05, 10000])
        plt.yscale("log")
        #plt.xlabel("Communication Score", size=24)
        #plt.ylabel("Uniqueness Score", size=24)
        #plt.title(joined_combinations[i], size=24, fontname="Helvetica Neue", fontweight="bold")
        plt.savefig(
            output_filepath+ 'only_in_one_or_other_'+joined_combinations[i]+'.pdf',
            bbox_inches="tight", pad_inches=0.5, transparent=True, dpi=300)
        only_condition1[i].to_csv(output_filepath + 'only_in_'+list_combinations[i][0]+'_in_'+joined_combinations[i]+'.csv', index=False)
        only_condition2 [i].to_csv(output_filepath + 'only_in_'+list_combinations[i][1]+'_in_'+joined_combinations[i]+'.csv', index=False)




        #This is the sum of all sigificant LR pairs within each condition
        newdict ={}
        for a in range(2):
            combined_dict2[conditions_numbers_combinations[i][a]] = combined_dict2[
                conditions_numbers_combinations[i][a]].copy(deep=True)
            combined_dict2[conditions_numbers_combinations[i][a]] = combined_dict2[conditions_numbers_combinations[i][a]].loc[
                (combined_dict2[conditions_numbers_combinations[i][a]]['Communication Score_' + list_combinations[i][a]] > 0)]
            combined2 = combined_dict2[conditions_numbers_combinations[i][a]].loc[combined_dict2[conditions_numbers_combinations[i][a]]["pval_"+list_combinations[i][a]] < pval_cutoff]
            combined3 = combined2.drop(combined2.iloc[:, 2:28], axis=1)
            combined3 = combined3.reset_index()
            newdict[a] = combined3.groupby(['Ligand', 'Receptor']).count()
            combined2.to_csv(
                output_filepath + 'save/' +
                joined_combinations[i] +'_in_condition_'+list_combinations[i][a]+ '.csv', index=False)
        logging.debug(newdict[0])
        newdict[0] = newdict[0].rename(columns={'index': list_combinations[i][0]})
        newdict[1] = newdict[1].rename(columns={'index': list_combinations[i][1]})
        df = pd.concat([newdict[0], newdict[1]], axis=1)
        df =df.fillna(0)
        df = df.reset_index()
        df["LR_name"] = df.Ligand.str.cat(df.Receptor, sep='^')
        df['avg'] = df[[list_combinations[i][0], list_combinations[i][1]]].mean(axis=1)
        df["stdev"] = df[[list_combinations[i][0], list_combinations[i][1]]].std(axis=1)
        df["rel_stdev"] = ((df["stdev"].div(df['avg']))*100).round(2)
        df["rel_stdev and rel_stdev_Rank"] = '(σ=' + df['rel_stdev'].astype(str) + ')'
        df = df.sort_values('avg')
        logging.debug(tabulate(df, headers="keys"))

        fig = plt.subplots(figsize =(16, 14))
        barWidth = 0.25

        LS = df[list_combinations[i][0]]
        NL = df[list_combinations[i][1]]


        br1 = np.arange(len(LS))
        br2 = [x + barWidth for x in br1]

        # Horizontal Bar Plot
        plt.barh(br1, LS, label = list_combinations[i][0], height=barWidth, color=color_combinations[i][0])
        plt.barh(br2, NL, label = list_combinations[i][1], height=barWidth, color=color_combinations[i][1])

        plt.xlabel('Number of LR Interactions', fontweight ='bold', fontsize = 24)
        plt.ylabel('LR Interactions', fontweight ='bold', fontsize = 24)
        plt.yticks([r + barWidth for r in range(len(LS))], df['LR_name'] + " " + df["rel_stdev and rel_stdev_Rank"], fontsize=20)
        plt.xticks(fontsize=20)
        plt.xlim([0, 32])

        plt.legend(fontsize=20, loc = 'lower right')
        plt.savefig(
            output_filepath + 'save/communication_sum_' +
            joined_combinations[i] + '.pdf',
            bbox_inches="tight", pad_inches=0.5, transparent=True, dpi=300)
        df.to_csv(output_filepath + 'save/'+joined_combinations[i]+'.csv', index=False)




        # This is the average strength of all significant LR pairs within each condition
        newdict ={}
        for a in range(2):
            combined_dict2[conditions_numbers_combinations[i][a]] = combined_dict2[conditions_numbers_combinations[i][a]].copy(deep=True)
            combined_dict2[conditions_numbers_combinations[i][a]] = combined_dict2[conditions_numbers_combinations[i][a]].loc[
                (combined_dict2[conditions_numbers_combinations[i][a]]['Communication Score_' + list_combinations[i][a]] > 0)]
            combined2 = combined_dict2[conditions_numbers_combinations[i][a]].loc[combined_dict2[conditions_numbers_combinations[i][a]]["pval_"+list_combinations[i][a]] < pval_cutoff]
            combined3 = combined2[['Ligand','Receptor','Communication Score_'+ list_combinations[i][a]]]
            combined3 = combined3.reset_index()
            newdict[a] = combined3.groupby(['Ligand', 'Receptor']).mean()
            combined2.to_csv(output_filepath + 'save2/' + joined_combinations[i] +'_in_condition_'+list_combinations[i][a]+ '.csv', index=False)
        logging.debug(newdict[0])
        newdict[0] = newdict[0].rename(columns={'index': list_combinations[i][0]})
        newdict[1] = newdict[1].rename(columns={'index': list_combinations[i][1]})
        df = pd.concat([newdict[0], newdict[1]], axis=1)
        df =df.fillna(0)
        df = df.reset_index()
        df["LR_name"] = df.Ligand.str.cat(df.Receptor, sep='^')
        df['avg'] = df[['Communication Score_' + list_combinations[i][0], 'Communication Score_' + list_combinations[i][1]]].mean(axis=1)
        df["stdev"] = df[['Communication Score_' + list_combinations[i][0], 'Communication Score_' + list_combinations[i][1]]].std(axis=1)
        df["rel_stdev"] = ((df["stdev"].div(df['avg']))*100).round(2)
        df["rel_stdev and rel_stdev_Rank"] = '(σ=' + df['rel_stdev'].astype(str) + ')'
        df = df.sort_values('avg')
        logging.debug(tabulate(df, headers="keys"))

        fig = plt.subplots(figsize =(18, 14))
        barWidth = 0.25

        LS = df['Communication Score_' + list_combinations[i][0]]
        NL = df['Communication Score_' + list_combinations[i][1]]

        br1 = np.arange(len(LS))
        br2 = [x + barWidth for x in br1]

        # Horizontal Bar Plot
        plt.barh(br1, LS, label = list_combinations[i][0], height=barWidth, color=color_combinations[i][0])
        plt.barh(br2, NL, label = list_combinations[i][1], height=barWidth, color=color_combinations[i][1])

        #plt.xlabel('Average Communication Score', fontweight ='bold', fontsize = 15)
        #plt.ylabel('LR Interactions', fontweight ='bold', fontsize = 15)
        plt.yticks([r + barWidth for r in range(len(LS))], df['LR_name']+ " " + df["rel_stdev and rel_stdev_Rank"], fontsize=20)
        plt.xticks(fontsize=20)
        plt.xlim([0, 10])

        plt.legend()
        plt.savefig(
            output_filepath + 'save2/communication_avg_strength_' +
            joined_combinations[i] + '.pdf',
            bbox_inches="tight", pad_inches=0.5, transparent=True, dpi=300)
        df.to_csv(output_filepath + 'save2/'+joined_combinations[i]+'.csv', index=False)





        #Making the chord diagrams for each condition within the two-way comparrison
        newdict2 ={}
        dis_dict = {}
        for a in range(2):
            combined_dict3[conditions_numbers_combinations[i][a]] = combined_dict3[conditions_numbers_combinations[i][a]].copy(deep=True)
            combined_dict3[conditions_numbers_combinations[i][a]] = combined_dict3[conditions_numbers_combinations[i][a]].loc[
                (combined_dict3[conditions_numbers_combinations[i][a]]['Communication Score_' + list_combinations[i][a]] > 0)]

            if a == 0:
                combined_dict3[conditions_numbers_combinations[i][a]] = combined_dict3[conditions_numbers_combinations[i][a]].loc[(combined_dict3[conditions_numbers_combinations[i][a]]['Communication Score_' + list_combinations[i][0]] > combined_dict3[conditions_numbers_combinations[i][a]]['Communication Score_' + list_combinations[i][1]])]
            elif a == 1:
                combined_dict3[conditions_numbers_combinations[i][a]] = combined_dict3[conditions_numbers_combinations[i][a]].loc[(combined_dict3[conditions_numbers_combinations[i][a]]['Communication Score_' + list_combinations[i][0]] < combined_dict3[conditions_numbers_combinations[i][a]]['Communication Score_' + list_combinations[i][1]])]

            combined2_3 = combined_dict3[conditions_numbers_combinations[i][a]].loc[combined_dict3[conditions_numbers_combinations[i][a]]["pval_"+list_combinations[i][a]] < pval_cutoff]
            dis_dict[a] = combined2_3
            combined2_4 = combined2_3.iloc[:,2:4]
            combined2_4 = combined2_4.reset_index()
            newdict2[a] = combined2_4.groupby(['Source', 'Target']).count()
        df2 = newdict2[0].rename(columns={'index':list_combinations[i][0]}).join(newdict2[1].rename(columns={'index':list_combinations[i][1]}),how='left')
        logging.debug(df2)
        df2 = df2.fillna(0)
        df2 = df2.reset_index()
        source = df2['Source'].tolist()
        target = df2['Target'].tolist()
        names = list(set(target + source))
        print("This is source", source)
        print("This is target", target)
        print("This is names", names)
        main_list_source = list(set(names) - set(source))
        main_list_target = list(set(names) - set(target))
        xtra = {'Source': main_list_source}
        xtra2 = {'Target': main_list_target}
        print("This is main_list_source", main_list_source)
        print("This is main_list_target", main_list_target)
        print("This is extra", xtra)
        print("This is extra", xtra2)
        print("This is df2",df2)
        df2 = pd.concat([df2, pd.DataFrame(xtra)], ignore_index=True)
        df2 = pd.concat([df2, pd.DataFrame(xtra2)], ignore_index=True)
        df2 = df2.fillna(0)
        logging.debug(tabulate(df2, headers="keys"))



        for a in range(2):
            logging.debug(tabulate(dis_dict[a], headers="keys"))
            this = movecol(dis_dict[a],
                         cols_to_move=['Source', 'Target'],
                         ref_col='Ligand',
                         place='Before').round(decimals=3)
            #this = dis_dict[a].iloc[:, [2,3,0,1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]].round(decimals=3)
            this = this[['Source', 'Target','Ligand','Receptor','Communication Score_'+list_combinations[i][a],'Uniqueness Score_'+list_combinations[i][a],"pval_"+list_combinations[i][a],'Pathway Label_'+list_combinations[i][a], 'rel_stdev','Source to Target']]
            this.to_csv(output_filepath + 'save3/' +"Source to Target Communication in "+joined_combinations[i]+' Significantly in '+list_combinations[i][a]+ '.csv', index=False)

            pivot = pd.pivot_table(df2, values=list_combinations[i][a],
                                        index='Source',
                                        columns='Target')
            logging.debug(pivot)
            pivot = pivot.fillna(0)
            pivot = pivot.to_numpy()
            pivot = np.round(pivot).astype(int)
            pivot = pivot.tolist()
            logging.debug(pivot)
            names = names
            #Chord(pivot, names, data_table=this.to_csv(index=False), directed=True, title = "Source --> Target Communication in "+joined_combinations[i]+' - Significant in '+list_combinations[i][a]).to_svg('/Users/Troks27/Desktop/9. Human Skin TrokaChat Analysis/Python TrokaChat7/TrokaChat Charts/For Figures/save3/'+"Source --> Target Communication in "+joined_combinations[i]+' Significantly in '+list_combinations[i][a]+".svg")
            Chord(pivot, names, directed=True).to_pdf(output_filepath + 'save3/'+"Source to Target Communication in "+joined_combinations[i]+' Significantly in '+list_combinations[i][a]+'.pdf')
            Chord(pivot, names, data_table=this.to_csv(index=False), directed=True).to_html(output_filepath + 'save3/'+"Source to Target Communication in "+joined_combinations[i]+' Significantly in '+list_combinations[i][a]+'.html')


conditions_list = conditions_list_clean
if conditions_to_test == 3:
    for i in range(3):
        combined_dict[i] = combined_dict[i].copy(deep=True)
        combined_dict[i] = combined_dict[i].loc[(combined_dict[i]['Communication Score_' + "LS"] > 0) & (combined_dict[i]['Communication Score_' + "NL"] > 0) & (combined_dict[i]['Communication Score_' + "H"] > 0)]
        combined_dict[i] = combined_dict[i].loc[(combined_dict[i]["pval_" + "LS"] < pval_cutoff) | (combined_dict[i]["pval_" + "NL"] < pval_cutoff) | (combined_dict[i]["pval_" + "H"] < pval_cutoff)]
        combined_dict[i] = combined_dict[i].loc[(combined_dict[i]["zscore_" + "LS"] > 0) | (combined_dict[i]["zscore_" + "NL"] > 0) | (combined_dict[i]["zscore_" + "H"] > 0)]
        combined_dict[i].sort_values("rel_stdev", ascending=False, inplace=True)

        combined_dict[i]["LR_name"] = combined_dict[i].Ligand.str.cat(combined_dict[i].Receptor, sep='^')
        cols = ["Source", "Target"]
        combined_dict[i]["Source to Target"] = combined_dict[i][cols].apply(lambda row: '-->'.join(row.values.astype(str)), axis=1)
        combined_dict[i] = combined_dict[i].head(top)
        logging.debug(combined_dict[i])



    #This is the creation of the dot plots with those pathways shared in both condtions
    combined_new = combined_dict[0].loc[(combined_dict[i]['Communication Score_' + "LS"] > 0) | (combined_dict[0]['Communication Score_' + "NL"] > 0)| (combined_dict[0]['Communication Score_' + "H"] > 0)]

    plt.figure(figsize=(8, 6))
    plt.scatter(combined_new['Communication Score_LS'],
                combined_new['Uniqueness Score_LS'],
                s=30,
                c="#FF1F5B",
                label="LS")

    plt.scatter(combined_new['Communication Score_NL'],
                combined_new['Uniqueness Score_NL'],
                s=30,
                c="#009ADE",
                label="NL")

    plt.scatter(combined_new['Communication Score_H'],
                combined_new['Uniqueness Score_H'],
                s=30,
                c="#00CD6C",
                label="H")

    plt.legend(loc="upper right")


    if fig_labels == 'Yes':
        labelLines(plt.gca().get_lines(), zorder=2.5, align=True, fontsize=5)
    else:
        logging.debug("Bruh")


    plt.yscale("log")
    plt.xlim([1,50])

    #plt.xlabel("Communication Score", size=14, fontweight = "bold")
    #plt.ylabel("Uniqueness Score", size=14, fontweight = "bold")
    #plt.title(joined_combinations[i], size=18, fontweight = "bold")
    plt.savefig(output_filepath + '3 way dot plot.pdf',
        bbox_inches="tight", pad_inches=0.5, transparent=True, dpi=300)
    combined_new.to_csv(output_filepath + '3 way dot plot.csv', index=False)







    # This is the average strength of all significant LR pairs within each condition
    conditions_list = conditions_list_clean
    newdict ={}
    for a in range(3):
        logging.debug(a)
        yes = combined_dict4[a].copy(deep=True)
        this = yes.loc[(yes['Communication Score_' + conditions_list[a]] > 0)]
        combined2 = this.loc[this["pval_"+conditions_list[a]] < pval_cutoff]
        combined3 = combined2[['Ligand','Receptor','Communication Score_'+ conditions_list[a]]]
        combined3 = combined3.reset_index()
        newdict[a] = combined3.groupby(['Ligand', 'Receptor']).mean()
        logging.debug(newdict[a])
    logging.debug(newdict[0].index)
    newdict[0] = newdict[0].rename(columns={'index':conditions_list[0]})
    newdict[1] = newdict[1].rename(columns={'index':conditions_list[1]})
    newdict[2] = newdict[2].rename(columns={'index':conditions_list[2]})
    df = pd.concat([newdict[0], newdict[1],newdict[2]], axis=1)
    logging.debug(tabulate(df, headers="keys"))
    df =df.fillna(0)
    df = df.reset_index()
    df["LR_name"] = df.Ligand.str.cat(df.Receptor, sep='^')
    df['avg'] = df[['Communication Score_' + conditions_list[0], 'Communication Score_' + conditions_list[1], 'Communication Score_' + conditions_list[2]]].mean(axis=1)
    df["stdev"] = df[['Communication Score_' + conditions_list[0], 'Communication Score_' + conditions_list[1], 'Communication Score_' + conditions_list[2]]].std(axis=1)
    df["rel_stdev"] = ((df["stdev"].div(df['avg']))*100).round(2)
    df["rel_stdev and rel_stdev_Rank"] = '(σ=' + df['rel_stdev'].astype(str) + ')'
    df = df.sort_values('avg')


    fig = plt.subplots(figsize =(18, 9))
    barWidth = 0.25

    H = df['Communication Score_' + conditions_list[0]]
    NL = df['Communication Score_' + conditions_list[1]]
    LS = df['Communication Score_' + conditions_list[2]]

    br1 = np.arange(len(H))
    br2 = [x + barWidth for x in br1]
    br3 = [x + barWidth for x in br2]

    # Horizontal Bar Plot
    plt.barh(br1, H, label = conditions_list[0], height=barWidth, color="#00CD6C")
    plt.barh(br2, NL, label = conditions_list[1], height=barWidth, color="#009ADE")
    plt.barh(br3, LS, label= conditions_list[2], height=barWidth, color="#FF1F5B")

    #plt.xlabel('Average Communication Score', fontweight ='bold', fontsize = 15)
    #plt.ylabel('LR Interactions', fontweight ='bold', fontsize = 15)
    plt.yticks([r + barWidth for r in range(len(LS))], df['LR_name']+ " " + df["rel_stdev and rel_stdev_Rank"])


    plt.legend()
    plt.savefig(output_filepath + 'save_3way/communication_avg_strength_'+ '.pdf', bbox_inches="tight", pad_inches=0.5, transparent=True, dpi=300)
    df.to_csv(output_filepath + 'save_3way/final'+'.csv', index=False)








