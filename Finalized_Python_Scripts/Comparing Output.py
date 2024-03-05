import pandas as pd
from tabulate import tabulate

conditions_list = ['NG','DIAB']
import_filepath = '/Volumes/LaCie/Bushra_study_2.28.24/TrokaChat/TrokaChat Output/'
export_filepath = '/Volumes/LaCie/Bushra_study_2.28.24/TrokaChat/TrokaChat Compare/'

#conditions_list = ['H','NL','LS']
#import_filepath = '/Users/Troks27/Desktop/TrokaChat for Publication - R Scripts/TrokaChat Publication/TrokaChat/CvsT/'
#export_filepath = '/Users/Troks27/Desktop/TrokaChat for Publication - R Scripts/TrokaChat Publication/TrokaChat/CvsT/'



d = {}
if len(conditions_list) == 3:
       letters = ['x','y','z']


       for i in range(len(conditions_list)):
              d[letters[i]] = pd.read_excel(import_filepath  + conditions_list[i] + "_allpathways.xlsx")
              d[letters[i]]['Sample Ident'] = conditions_list[i]


       int_y_z_x = d[letters[1]].merge(d[letters[0]], on=['Ligand', 'Receptor', 'Source', 'Target'], how='inner').merge(d[letters[2]], on=['Ligand', 'Receptor', 'Source', 'Target'], how='inner').reset_index(drop=True)
       print(tabulate(int_y_z_x, headers="keys"))
       y_x = d[letters[1]].merge(d[letters[0]],on=['Ligand', 'Receptor', 'Source', 'Target'], how='inner').reset_index(drop=True)
       y_x = int_y_z_x.merge(y_x,on=['Ligand', 'Receptor', 'Source', 'Target'] , how='outer', indicator=True)
       y_x = y_x.loc[y_x._merge == 'right_only'].drop(['_merge'], axis=1)

       y_z = d[letters[1]].merge(d[letters[2]],on=['Ligand', 'Receptor', 'Source', 'Target'], how='inner').reset_index(drop=True)
       y_z = int_y_z_x.merge(y_z,on=['Ligand', 'Receptor', 'Source', 'Target'] , how='outer', indicator=True)
       y_z = y_z.loc[y_z._merge == 'right_only'].drop(['_merge', 'jank_mass_action','Uniqueness Score', 'Sample Ident'], axis=1)

       x_z = d[letters[2]].merge(d[letters[0]],on=['Ligand', 'Receptor', 'Source', 'Target'], how='inner').reset_index(drop=True)
       x_z = int_y_z_x.merge(x_z,on=['Ligand', 'Receptor', 'Source', 'Target'] , how='outer', indicator=True)
       x_z = x_z.loc[x_z._merge == 'right_only'].drop(['_merge', 'jank_mass_action','Uniqueness Score', 'Sample Ident'], axis=1).reset_index(drop=True)

       int_y_x = y_x
       int_y_z = y_z
       int_x_z = x_z

       y_final = d[letters[1]].merge(int_y_x,on=['Ligand', 'Receptor', 'Source', 'Target'], how='outer', indicator=True)
       y_final = y_final.loc[y_final._merge == 'left_only'].drop(['_merge'], axis=1)
       y_final = y_final.merge(int_y_z,on=['Ligand', 'Receptor', 'Source', 'Target'] , how='outer', indicator=True)
       y_final = y_final.loc[y_final._merge == 'left_only'].drop(['_merge'], axis=1)
       y_final = y_final.merge(int_y_z_x,on=['Ligand', 'Receptor', 'Source', 'Target'], how='outer', indicator=True)
       y_final = y_final.loc[y_final._merge == 'left_only'].drop(['_merge'], axis=1).reset_index(drop=True)


       x_final = d[letters[0]].merge(int_y_x,on=['Ligand', 'Receptor', 'Source', 'Target'] , how='outer', indicator=True)
       x_final = x_final.loc[x_final._merge == 'left_only'].drop(['_merge'], axis=1)
       x_final = x_final.merge(int_x_z,  on=['Ligand', 'Receptor', 'Source', 'Target'], how='outer', indicator=True)
       x_final = x_final.loc[x_final._merge == 'left_only'].drop(['_merge'], axis=1)
       x_final = x_final.merge(int_y_z_x,on=['Ligand', 'Receptor', 'Source', 'Target'] , how='outer', indicator=True)
       x_final = x_final.loc[x_final._merge == 'left_only'].drop(['_merge'], axis=1).reset_index(drop=True)


       z_final = d[letters[2]].merge(int_y_z,on=['Ligand', 'Receptor', 'Source', 'Target'] , how='outer', indicator=True)
       z_final = z_final.loc[z_final._merge == 'left_only'].drop(['_merge'], axis=1)
       z_final = z_final.merge(int_x_z,on=['Ligand', 'Receptor', 'Source', 'Target'] , how='outer', indicator=True)
       z_final = z_final.loc[z_final._merge == 'left_only'].drop(['_merge'], axis=1)
       z_final = z_final.merge(int_y_z_x,on=['Ligand', 'Receptor', 'Source', 'Target'] , how='outer', indicator=True)
       z_final = z_final.loc[z_final._merge == 'left_only'].drop(['_merge'], axis=1).reset_index(drop=True)

       int_y_x = int_y_x.drop(['jank_mass_action_x_x','Uniqueness Score_x_x', 'Sample Ident_x_x', 'jank_mass_action_y_x', 'Uniqueness Score_y_x', 'Sample Ident_y_x', 'jank_mass_action','Uniqueness Score', 'Sample Ident', 'Pathway Label_x_x','Pathway Label_y_x','Pathway Label'], axis=1).reset_index(drop=True)
       print(tabulate(int_y_x, headers="keys"))
       int_x_z = int_x_z.drop(['jank_mass_action_x_x','Uniqueness Score_x_x', 'Sample Ident_x_x', 'jank_mass_action_y_x', 'Uniqueness Score_y_x', 'Sample Ident_y_x','Pathway Label_x_x','Pathway Label_y_x','Pathway Label'], axis=1).reset_index(drop=True)
       print(tabulate(int_x_z, headers="keys"))
       int_y_z = int_y_z.drop(['jank_mass_action_x_x','Uniqueness Score_x_x', 'Sample Ident_x_x', 'jank_mass_action_y_x', 'Uniqueness Score_y_x', 'Sample Ident_y_x', 'Pathway Label_x_x','Pathway Label_y_x', 'Pathway Label'], axis=1).reset_index(drop=True)


       col1 = [str(int_y_x.iat[0, 7]) + "/" + str(int_y_x.iat[0, 11])]
       col2 = [str(int_y_x.iat[0, 11]) + "/" + str(int_y_x.iat[0, 7])]
       int_y_x[col1[0]] = int_y_x.iloc[:, 4] / int_y_x.iloc[:, 8]
       int_y_x[col2[0]] = int_y_x.iloc[:, 8] / int_y_x.iloc[:, 4]



       col3 = [str(int_x_z.iat[0, 7]) + "/" + str(int_x_z.iat[0, 11])]
       col4 = [str(int_x_z.iat[0, 11]) + "/" + str(int_x_z.iat[0, 7])]
       int_x_z[col3[0]] = int_x_z.iloc[:, 4] / int_x_z.iloc[:, 8]
       int_x_z[col4[0]] = int_x_z.iloc[:, 8] / int_x_z.iloc[:, 4]



       col5 = [str(int_y_z.iat[0, 7]) + "/" + str(int_y_z.iat[0, 11])]
       col6 = [str(int_y_z.iat[0, 11]) + "/" + str(int_y_z.iat[0, 7])]
       int_y_z[col6[0]] = int_y_z.iloc[:, 4] / int_y_z.iloc[:, 8]
       int_y_z[col5[0]] = int_y_z.iloc[:, 8] / int_y_z.iloc[:, 4]



       col7 = [str(int_y_z_x.iat[0, 7]) + "/" + str(int_y_z_x.iat[0, 11])]
       col8 = [str(int_y_z_x.iat[0, 7]) + "/" + str(int_y_z_x.iat[0, 15])]
       col9 = [str(int_y_z_x.iat[0, 11]) + "/" + str(int_y_z_x.iat[0, 15])]
       int_y_z_x[col7[0]] = int_y_z_x.iloc[:, 4] / int_y_z_x.iloc[:, 8]
       int_y_z_x[col8[0]] = int_y_z_x.iloc[:, 4] / int_y_z_x.iloc[:, 12]
       int_y_z_x[col9[0]] = int_y_z_x.iloc[:, 8] / int_y_z_x.iloc[:, 12]




       y_final = y_final.drop(['jank_mass_action_x_x_x',
              'Uniqueness Score_x_x_x', 'Sample Ident_x_x_x',
              'jank_mass_action_y_x_x', 'Uniqueness Score_y_x_x',
              'Sample Ident_y_x_x', 'jank_mass_action_y_x', 'Uniqueness Score_y_x',
              'Sample Ident_y_x', 'jank_mass_action_x_y_x', 'Uniqueness Score_x_y_x',
              'Sample Ident_x_y_x', 'jank_mass_action_y_y_x',
              'Uniqueness Score_y_y_x', 'Sample Ident_y_y_x',
              'jank_mass_action_x_x_y', 'Uniqueness Score_x_x_y',
              'Sample Ident_x_x_y', 'jank_mass_action_y_x_y',
              'Uniqueness Score_y_x_y', 'Sample Ident_y_x_y',
              'jank_mass_action_x_y_y', 'Uniqueness Score_x_y_y',
              'Sample Ident_x_y_y', 'jank_mass_action_y_y_y',
              'Uniqueness Score_y_y_y', 'Sample Ident_y_y_y', 'jank_mass_action_x_y',
              'Uniqueness Score_x_y', 'Sample Ident_x_y', 'jank_mass_action_y_y',
              'Uniqueness Score_y_y', 'Sample Ident_y_y', 'jank_mass_action',
              'Uniqueness Score', 'Sample Ident'], axis = 1)

       x_final = x_final.drop(['jank_mass_action_x_x_x',
              'Uniqueness Score_x_x_x', 'Sample Ident_x_x_x',
              'jank_mass_action_y_x_x', 'Uniqueness Score_y_x_x',
              'Sample Ident_y_x_x', 'jank_mass_action_y_x', 'Uniqueness Score_y_x',
              'Sample Ident_y_x', 'jank_mass_action_x_y_x', 'Uniqueness Score_x_y_x',
              'Sample Ident_x_y_x', 'jank_mass_action_y_y_x',
              'Uniqueness Score_y_y_x', 'Sample Ident_y_y_x',
              'jank_mass_action_x_x_y', 'Uniqueness Score_x_x_y',
              'Sample Ident_x_x_y', 'jank_mass_action_y_x_y',
              'Uniqueness Score_y_x_y', 'Sample Ident_y_x_y',
              'jank_mass_action_x_y_y', 'Uniqueness Score_x_y_y',
              'Sample Ident_x_y_y', 'jank_mass_action_y_y_y',
              'Uniqueness Score_y_y_y', 'Sample Ident_y_y_y', 'jank_mass_action_x_y',
              'Uniqueness Score_x_y', 'Sample Ident_x_y', 'jank_mass_action_y_y',
              'Uniqueness Score_y_y', 'Sample Ident_y_y', 'jank_mass_action',
              'Uniqueness Score', 'Sample Ident'], axis = 1)

       z_final = z_final.drop(['jank_mass_action_x_x_x',
              'Uniqueness Score_x_x_x', 'Sample Ident_x_x_x',
              'jank_mass_action_y_x_x', 'Uniqueness Score_y_x_x',
              'Sample Ident_y_x_x', 'jank_mass_action_x_y_x',
              'Uniqueness Score_x_y_x', 'Sample Ident_x_y_x',
              'jank_mass_action_y_y_x', 'Uniqueness Score_y_y_x',
              'Sample Ident_y_y_x', 'jank_mass_action_x_x_y',
              'Uniqueness Score_x_x_y', 'Sample Ident_x_x_y',
              'jank_mass_action_y_x_y', 'Uniqueness Score_y_x_y',
              'Sample Ident_y_x_y', 'jank_mass_action_x_y_y',
              'Uniqueness Score_x_y_y', 'Sample Ident_x_y_y',
              'jank_mass_action_y_y_y', 'Uniqueness Score_y_y_y',
              'Sample Ident_y_y_y', 'jank_mass_action_y', 'Uniqueness Score_y',
              'Sample Ident_y', 'jank_mass_action_y', 'Uniqueness Score_y',
              'Sample Ident_y'], axis = 1)


       int_int_x_y = int(len(int_y_x))
       int_int_x_z = int(len(int_x_z))
       int_int_y_z = int(len(int_y_z))
       int_int_y_z_x = int(len(int_y_z_x))
       int_unique_x = int(len(x_final))
       int_unique_y = int(len(y_final))
       int_unique_z = int(len(z_final))

       print([int_int_x_y,int_int_x_z,int_int_y_z,int_int_y_z_x,int_unique_x,int_unique_y,int_unique_z])

       #Import libraries
       from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
       from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
       from matplotlib import pyplot as plt
       import matplotlib.pyplot as plt
       import time


       venn3_unweighted(subsets = (int_unique_x, int_unique_y, int_int_x_y, int_unique_z, int_int_x_z, int_int_y_z, int_int_y_z_x), set_labels = conditions_list, alpha = 0.5);

       plt.title("Signaling Group Comparisons")
       plt.savefig(export_filepath+'Signaling Group Comparisons.png', bbox_inches ='tight', dpi = 300, transparent = True)
       plt.clf()

       print(int_unique_z + int_int_y_z_x + int_int_y_z + int_int_x_z)
       print(int_unique_x + int_int_y_z_x + int_int_x_y + int_int_x_z)
       print(int_unique_y + int_int_y_z_x + int_int_x_y + int_int_y_z)


       int_y_x.to_csv(export_filepath+conditions_list[1]+"_"+conditions_list[0]+'_output_int.csv', index=False)
       int_x_z.to_csv(export_filepath+conditions_list[0]+"_"+conditions_list[2]+'_output_int.csv', index=False)
       int_y_z.to_csv(export_filepath+conditions_list[1]+"_"+conditions_list[2]+'_output_int.csv', index=False)
       int_y_z_x.to_csv(export_filepath+conditions_list[1]+"_"+conditions_list[2]+"_"+conditions_list[0]+'_output_int.csv', index=False)
       x_final.to_csv(export_filepath+conditions_list[0]+'_outputfinal.csv', index=False)
       y_final.to_csv(export_filepath+conditions_list[1]+'_outputfinal.csv', index=False)
       z_final.to_csv(export_filepath+conditions_list[2]+'_outputfinal.csv', index=False)
else:
       letters = ['x', 'y']

       for i in range(len(conditions_list)):
              d[letters[i]] = pd.read_excel(import_filepath + conditions_list[i] + "_allpathways.xlsx")
              d[letters[i]]['Sample Ident'] = conditions_list[i]

       int_y_x = d[letters[1]].merge(d[letters[0]], on=['Ligand', 'Receptor', 'Source', 'Target'], how='inner').reset_index(drop=True)
       print(tabulate(int_y_x, headers="keys"))

       y_final = d[letters[1]].merge(int_y_x, on=['Ligand', 'Receptor', 'Source', 'Target'], how='outer',
                                     indicator=True)
       y_final = y_final.loc[y_final._merge == 'left_only'].drop(['_merge'], axis=1).reset_index(drop=True)
       print(tabulate(y_final, headers="keys"))
       x_final = d[letters[0]].merge(int_y_x, on=['Ligand', 'Receptor', 'Source', 'Target'], how='outer',
                                     indicator=True)
       x_final = x_final.loc[x_final._merge == 'left_only'].drop(['_merge'], axis=1).reset_index(drop=True)
       print(tabulate(x_final, headers="keys"))

       col1 = [str(int_y_x.iat[0, 7]) + "/" + str(int_y_x.iat[0, 11])]
       col2 = [str(int_y_x.iat[0, 11]) + "/" + str(int_y_x.iat[0, 7])]
       int_y_x[col1[0]] = int_y_x.iloc[:, 4] / int_y_x.iloc[:, 8]
       int_y_x[col2[0]] = int_y_x.iloc[:, 8] / int_y_x.iloc[:, 4]

       int_int_x_y = int(len(int_y_x))
       int_unique_x = int(len(x_final))
       int_unique_y = int(len(y_final))


       from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
       from matplotlib import pyplot as plt
       import matplotlib.pyplot as plt

       venn2_unweighted(
              subsets=(int_unique_y, int_unique_x, int_int_x_y),
              set_labels= reversed(conditions_list), alpha=0.5);

       plt.title("Signaling Group Comparisons")
       plt.savefig(export_filepath + 'Signaling Group Comparisons.png', bbox_inches='tight', dpi=300, transparent=True)
       plt.clf()

       print(y_final + int_y_x)
       print(x_final + int_y_x)

       int_y_x.to_csv(export_filepath + conditions_list[1] + "_" + conditions_list[0] + '_output_int.csv', index=False)
       x_final.to_csv(export_filepath + conditions_list[0] + '_outputfinal.csv', index=False)
       y_final.to_csv(export_filepath + conditions_list[1] + '_outputfinal.csv', index=False)
