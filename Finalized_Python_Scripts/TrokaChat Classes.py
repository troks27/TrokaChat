# Import required packages
import scipy
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt

#Class Definition
class TrokaChat_import:
    def __init__(self, input_path: str, file_prefixes: list, output_path: str, pct_min: float):
        """
        Initializes the TrokaChat instance.

        Args:
            input_folder (str): The folder where the input CSV files are located.
            file_prefixes (list): The list of prefixes for the CSV files to be loaded.
            path (str): The path to the input and output folders.
            output_folder (str): The folder where the output Excel files will be saved.
            pct_min (float): The minimum percentage for filtering the data.
        """
        self.input_path = input_path
        self.file_prefixes = file_prefixes
        self.output_path = output_path
        self.pct_min = pct_min
        self.frame = None

    def load_data(self):
        """
        Loads data from CSV files into a pandas DataFrame.

        For each file prefix, it loads the corresponding CSV file and adds a 'sample' column with the condition derived from the prefix.
        """
        frames = []
        for prefix in self.file_prefixes:
            df = pd.read_csv(f"{self.input_path}/{prefix}.csv")
            # Extract the condition from the file_prefix
            if prefix.startswith('nulldistavgs_'):
                condition = 'nulldistavgs_' + prefix.split('_')[1]
            else:
                condition = prefix.split('_')[0]
            df['sample'] = condition
            frames.append(df)
        self.frame = pd.concat(frames)

    def filter_transform_data(self):
        """
        Filters and transforms the data.

        It filters the data by 'pct.1' column and performs transformations on the 'avg_log2FC' column.
        """
        self.frame = self.frame[self.frame["pct.1"] > self.pct_min]
        self.frame['gene2'] = self.frame['gene']
        self.frame['avg_log2FC'] = ((self.frame['avg_log2FC'] - self.frame['avg_log2FC'].min()) *
                                    (self.frame['avg_log2FC'].max() + abs(self.frame['avg_log2FC'].min()) + 1 -
                                    (self.frame['avg_log2FC'].min() + abs(self.frame['avg_log2FC'].min()) + 1)) /
                                    (self.frame['avg_log2FC'].max() - self.frame['avg_log2FC'].min()) +
                                    (self.frame['avg_log2FC'].min() + abs(self.frame['avg_log2FC'].min()) + 1))
        self.frame.drop(columns = ['gene'], inplace=True)
        self.frame.rename(columns={'gene2': 'gene'}, inplace=True)

    def plot_data(self):
        """
        Plots the histogram of the 'avg_log2FC' column.
        """
        plt.hist(self.frame["avg_log2FC"], bins=25, density=True, alpha=0.6, color='g')
        plt.show()

    def describe_data(self):
        """
        Prints the summary statistics of the DataFrame.
        """
        print(self.frame.describe())

    def export_data(self):
        """
        Exports the data for each condition to a separate Excel file.
        """
        for condition, df in self.frame.groupby('sample'):
            df.to_excel(f"{self.output_path}/{condition} DEGs.xlsx", index=False)

    def import_data(self):
        """
        Performs the complete data processing workflow.

        It includes loading, filtering, transforming, plotting, describing, and exporting the data.
        """
        self.load_data()
        self.filter_transform_data()
        self.plot_data()
        self.describe_data()
        self.export_data()





#Class Definition
class TrokaChat_Computation:
    def __init__(self, import_folder, general_filepath, output_path, files_list, species, sample_name_from_R, counts_file, min_clus_percent):
        self.import_folder = import_folder
        self.general_filepath = general_filepath
        self.output_path = output_path
        self.files_list = files_list
        self.species = species
        self.sample_name_from_R = sample_name_from_R
        self.counts_file = counts_file
        self.min_clus_percent = min_clus_percent

    def load_reference(self, species):
        """
        This function loads the reference Excel file based on the species provided.
        Currently supports 'human' and 'mouse'.

        Parameters:
        species (str): The species for which to load the reference. Either 'human' or 'mouse'.

        Returns:
        pandas.DataFrame: DataFrame containing ligand-receptor data from the appropriate reference Excel file.
        """

        # Choose the reference list based on the species
        if species == 'human':
            reference_list = 'human_ALL_PATHWAYS_UPDATED.xlsx'
        elif species == 'mouse':
            reference_list = 'All Signaling Pathways.xlsx'
        else:
            raise ValueError("Invalid species. Expected 'human' or 'mouse'.")

        # Load the reference data
        lr_db_df = pd.read_excel(reference_list, sheet_name=0)

        return lr_db_df

    def preprocess_counts_file(self, min_clus_percent):
        """
        Preprocess the counts file and return a list of clusters.

        This function reads the counts data from a CSV file, replaces missing values with 0, and calculates the
        percentage of counts per sample. The columns of the DataFrame are ordered numerically, and then the function
        filters out the clusters that have less than 'min_clus_percent' percent of counts in all samples.

        The function returns a list of clusters that have more than 'min_clus_percent' percent of counts
        in at least one sample.

        Parameters
        ----------
        min_clus_percent : float
            The minimum percentage of counts for a cluster to be included in the output.

        Returns
        -------
        list
            A list of clusters that have more than 'min_clus_percent' percent of counts in at least one sample.
        """

        # Read counts data from the CSV file, replace missing values with 0, and set the index to sample names
        counts = pd.read_csv(f'{self.general_filepath}{self.import_folder}/{self.counts_file}.csv').fillna(0).set_index(self.sample_name_from_R)

        # Calculate percentage of counts per sample
        counts_percent = counts.div(counts.sum(axis=1), axis=0).multiply(100)

        # Order columns numerically
        ordered_cols = sorted(int(x) for x in counts_percent.columns)
        counts_percent = counts_percent[[str(x) for x in ordered_cols]]

        # Filter out the clusters that have less than 'min_clus_percent' percent of counts in all samples
        new_counts = counts_percent.loc[:, (counts_percent > min_clus_percent).any(axis=0)]

        # Convert column names (clusters) to a list of integers
        new_counts_list = [int(i) for i in new_counts.columns.values]

        return new_counts_list

    def create_lr_dict(self, lr_db_df):
        """
        Create a dictionary representing the ligand-receptor database.

        The dictionary is structured as follows:
        {'Ligand': {'Receptor': ['Subunit_1', 'Subunit_2', ..., 'Subunit_n']}, ...}

        Parameters
        ----------
        lr_db_df : pandas.DataFrame
            DataFrame containing the ligand-receptor database. The dataframe should contain 'ligand', 'receptor',
            and 'subunit_i' (where i ranges from 1 to n) columns.

        Returns
        -------
        dict
            A dictionary representation of the ligand-receptor database.
        """

        # Initialize the dictionary
        lr_dict = {}

        # Iterate over rows in the DataFrame
        for _, row in lr_db_df.iterrows():
            # Extract ligand and receptor from the row
            ligand = row['ligand']
            receptor = row['receptor']

            # Extract subunits; note that if a subunit value is NaN, it is not included
            subunits = [row[f'subunit_{i}'] for i in range(1, 5) if pd.notnull(row[f'subunit_{i}'])]

            # If the ligand is not already a key in the dictionary, add it and initialize its value as an empty dict
            if ligand not in lr_dict:
                lr_dict[ligand] = {}

            # Add the receptor and its corresponding subunits to the ligand's dict
            lr_dict[ligand][receptor] = subunits

        # Return the completed ligand-receptor dictionary
        return lr_dict

    def calculate_communication_score(self, row, deg_data):
        """
        Calculate the communication score for a ligand-receptor pair.

        This function calculates the communication score based on the expression changes (avg_log2FC) and the percentage
        of cells expressing the gene (pct.1) for both the ligand and its receptor subunits. The communication score is
        defined as the product of the ligand's score and the average score of the receptor subunits.

        Parameters
        ----------
        row : pandas.Series
            A row of a DataFrame, which represents a ligand-receptor pair. The row should have the columns 'Ligand',
            'Subunits', 'Source' and 'Target', representing the ligand name, a list of subunits names, the source cluster
            and the target cluster respectively.
        deg_data : pandas.DataFrame
            A DataFrame with differentially expressed genes data. The DataFrame should have the columns 'gene', 'cluster',
            'avg_log2FC' and 'pct.1', representing the gene name, the cluster number, the average log2 fold change, and
            the percentage of cells expressing the gene respectively.

        Returns
        -------
        float
            The communication score for the ligand-receptor pair.
        """

        # Get ligand name and subunits from the row
        ligand = row['Ligand']
        subunits = row['Subunits']

        # Filter the differentially expressed genes data for the ligand and its source cluster
        ligand_data = deg_data[(deg_data['gene'] == ligand) & (deg_data['cluster'] == row['Source'])]

        # Calculate the ligand score as the average log2 fold change raised to the power of the percentage of cells expressing the ligand
        ligand_score = (ligand_data['avg_log2FC'].values[0]) ** (ligand_data['pct.1'].values[0])

        # Initialize an empty list to store the scores of the subunits
        subunit_scores = []

        # Iterate over each subunit
        for subunit in subunits:
            # Filter the differentially expressed genes data for the subunit and its target cluster
            subunit_data = deg_data[(deg_data['gene'] == subunit) & (deg_data['cluster'] == row['Target'])]

            # Calculate the subunit score as the average log2 fold change raised to the power of the percentage of cells expressing the subunit
            subunit_score = (subunit_data['avg_log2FC'].values[0]) ** (subunit_data['pct.1'].values[0])

            # Append the subunit score to the list
            subunit_scores.append(subunit_score)

        # Calculate the average score of the subunits
        average_subunit_score = sum(subunit_scores) / len(subunit_scores)

        # Calculate the communication score as the product of the ligand's score and the average score of the subunits
        communication_score = ligand_score * average_subunit_score

        return communication_score

    def calculate_uniqueness_score(self, row, deg_data):
        """
        Calculate the uniqueness score for a ligand-receptor pair.

        This function calculates the uniqueness score based on the specificity of expression (pct.1/pct.2) for both
        the ligand and its receptor subunits. The uniqueness score is defined as the product of the ligand's score
        and the average score of the receptor subunits.

        Parameters
        ----------
        row : pandas.Series
            A row of a DataFrame, which represents a ligand-receptor pair. The row should have the columns 'Ligand',
            'Subunits', 'Source' and 'Target', representing the ligand name, a list of subunits names, the source cluster
            and the target cluster respectively.
        deg_data : pandas.DataFrame
            A DataFrame with differentially expressed genes data. The DataFrame should have the columns 'gene', 'cluster',
            'pct.1' and 'pct.2', representing the gene name, the cluster number, the percentage of cells expressing the gene
            in the cluster, and the percentage of cells expressing the gene in all clusters respectively.

        Returns
        -------
        float
            The uniqueness score for the ligand-receptor pair.
        """

        # Get ligand name and subunits from the row
        ligand = row['Ligand']
        subunits = row['Subunits']

        # Filter the differentially expressed genes data for the ligand and its source cluster
        ligand_data = deg_data[(deg_data['gene'] == ligand) & (deg_data['cluster'] == row['Source'])]

        # Calculate the ligand score as the ratio of the percentage of cells expressing the ligand in the cluster to
        # the percentage of cells expressing the ligand in all clusters
        ligand_score = ligand_data['pct.1'].values[0] / ligand_data['pct.2'].values[0]

        # Initialize an empty list to store the scores of the subunits
        subunit_scores = []

        # Iterate over each subunit
        for subunit in subunits:
            # Filter the differentially expressed genes data for the subunit and its target cluster
            subunit_data = deg_data[(deg_data['gene'] == subunit) & (deg_data['cluster'] == row['Target'])]

            # Calculate the subunit score as the ratio of the percentage of cells expressing the subunit in the cluster to
            # the percentage of cells expressing the subunit in all clusters
            subunit_score = subunit_data['pct.1'].values[0] / subunit_data['pct.2'].values[0]

            # Append the subunit score to the list
            subunit_scores.append(subunit_score)

        # Calculate the average score of the subunits
        average_subunit_score = sum(subunit_scores) / len(subunit_scores)

        # Calculate the uniqueness score as the product of the ligand's score and the average score of the subunits
        uniqueness_score = ligand_score * average_subunit_score

        return uniqueness_score

    def find_matching_lr_pairs(self, deg_data, lr_dict):
        '''
        Determines which ligand-receptor interactions are present in input data (deg_data)
        based on an input dictionary of known ligand-receptor interactions. It considers all possible
        combinations of ligand-receptor pairs, but only adds pairs where all subunits of the receptor
        come from the same cluster.

        Parameters:
        ----------
        deg_data : pandas.DataFrame
            Input dataframe of DEG data processed from a Seurat object in R.
        lr_dict : dict
            Dictionary of ligand-receptor interactions from the CellChat database.

        Returns:
        -------
        matching_pairs_df : pandas.DataFrame
            A dataframe of ligand-receptor pairs where all receptor subunits come from the same cluster.
        '''

        # Initialize a DataFrame to hold matching pairs
        matching_pairs_df = pd.DataFrame(
            columns=['Ligand', 'Receptor', 'Subunits', 'Source', 'Target', 'Communication Score', 'Uniqueness Score'])

        # Initialize an empty list to hold row data
        rows_list = []

        # Iterate over ligands and their associated receptors in the ligand-receptor dictionary
        for ligand, receptors in lr_dict.items():
            if ligand in deg_data['gene'].values: # Only process ligands present in the DEG data
                ligand_clusters = deg_data.loc[deg_data['gene'] == ligand, 'cluster'].values

                # Iterate over receptors and their associated subunits for a given ligand
                for receptor, subunits in receptors.items():
                    receptor_subunits_data = deg_data[deg_data['gene'].isin(subunits)]
                    receptor_clusters = receptor_subunits_data['cluster'].unique()

                    # Iterate over unique receptor clusters
                    for receptor_cluster in receptor_clusters:
                        # Check if all subunits of the receptor come from the same cluster
                        if set(subunits).issubset(
                                receptor_subunits_data[receptor_subunits_data['cluster'] == receptor_cluster][
                                    'gene'].values):
                            # Iterate over ligand clusters
                            for ligand_cluster in ligand_clusters:
                                # Create a new row for the matching pair
                                new_row = {
                                    'Ligand': ligand,
                                    'Receptor': receptor,
                                    'Subunits': subunits,
                                    'Source': ligand_cluster,
                                    'Target': receptor_cluster
                                }
                                # Append the new row data to the list
                                rows_list.append(new_row)

        # Convert the list of rows into a DataFrame
        matching_pairs_df = pd.DataFrame(rows_list)

        # Calculate the scores for each row and add them to the DataFrame
        for i, row in matching_pairs_df.iterrows():
            communication_score = self.calculate_communication_score(row, deg_data)
            uniqueness_score = self.calculate_uniqueness_score(row, deg_data)

            matching_pairs_df.loc[i, 'Communication Score'] = communication_score
            matching_pairs_df.loc[i, 'Uniqueness Score'] = uniqueness_score

        return matching_pairs_df

    def process_files(self):
        """
        Processes a list of files containing differentially expressed genes data.

        For each file, this function reads the data, preprocesses it, creates a ligand-receptor dictionary from
        a ligand-receptor database, finds the matching ligand-receptor pairs, calculates their respective
        communication and uniqueness scores, and writes the results to an Excel file.

        Parameters
        ----------
        files_list : list
            A list of file names (strings) to be processed. Each file should be an Excel file
            containing differentially expressed genes data.
        species : str
            The species for which to load the reference. Either 'human' or 'mouse'.

        Returns
        -------
        None. The function writes the results to Excel files.
        """

        # Load the ligand-receptor reference data for the specified species
        lr_db_df = self.load_reference(self.species)

        # Iterate over all files in the list
        for file_name in tqdm(self.files_list, desc="Processing files"):
            # Read the differentially expressed genes data from the Excel file
            deg_data = pd.read_excel(self.general_filepath + self.import_folder + '/' + file_name + '.xlsx', sheet_name=0)

            # Preprocess the differentially expressed genes data
            deg_data = deg_data.loc[deg_data['avg_log2FC'] > 0]
            deg_data['pct.2'] = deg_data['pct.2'].replace(0, 0.001)
            deg_data = deg_data.reset_index(drop=True)

            # Preprocess the counts file
            new_counts_list = self.preprocess_counts_file(self.min_clus_percent)

            # Create the ligand-receptor dictionary from the ligand-receptor database
            lr_dict = self.create_lr_dict(lr_db_df)

            # Find the matching ligand-receptor pairs and calculate the scores
            matching_pairs_df = self.find_matching_lr_pairs(deg_data, lr_dict)

            # Filter the matching pairs based on the counts and uniqueness score
            matching_pairs_df = matching_pairs_df.set_index(['Target'])
            matching_pairs_df = matching_pairs_df.loc[matching_pairs_df.index.isin(new_counts_list)].reset_index()
            matching_pairs_df = matching_pairs_df.set_index(['Source'])
            matching_pairs_df = matching_pairs_df.loc[matching_pairs_df.index.isin(new_counts_list)].reset_index()
            matching_pairs_df = matching_pairs_df[
                ['Ligand', 'Receptor', 'Source', 'Target', 'Communication Score', 'Uniqueness Score']]
            matching_pairs_df = matching_pairs_df.loc[matching_pairs_df["Uniqueness Score"] > 0]

            matching_pairs_df = pd.merge(matching_pairs_df, lr_db_df[['ligand', 'receptor', 'annotation']],
                                         left_on=['Ligand', 'Receptor'],
                                         right_on=['ligand', 'receptor'],
                                         how='left')

            # drop the additional columns 'ligand' and 'receptor' from lr_db_df
            matching_pairs_df.drop(['ligand', 'receptor'], axis=1, inplace=True)
            matching_pairs_df = matching_pairs_df.rename(columns={'annotation': 'Pathway Label'})

            # Generate the output file name
            if file_name.startswith('nulldistavgs_'):
                matching_pairs_df['Sample Ident'] = f"null_{file_name.split('_')[1]}"
                output_file_name = f"NULL_DIST_{file_name.split('_')[1]}_allpathways.xlsx"
            else:
                matching_pairs_df['Sample Ident'] = f"{file_name.split(' ')[0]}"
                output_file_name = f"{file_name.split(' ')[0]}_allpathways.xlsx"

            # Write the matching pairs data to an Excel file
            matching_pairs_df.to_excel(f'{self.output_path}{output_file_name}', index=False)



filename = TrokaChat_import(
    input_path = '/Users/troks27/Desktop/TrokaChat All for Publication/TrokaChat/TrokaChat/csv1',
    file_prefixes = ['H_vs_H','NL_vs_H','LS_vs_H', 'nulldistavgs_H_vs_H','nulldistavgs_NL_vs_H','nulldistavgs_LS_vs_H'],
    output_path = '/Users/troks27/Desktop/TrokaChat All for Publication/TrokaChat/TrokaChat',
    pct_min = 0.1)

#filename.import_data()

TrokaObject = TrokaChat_Computation(
    import_folder = 'exp',
    general_filepath = '/Users/troks27/Desktop/TrokaChat All for Publication/TrokaChat/TrokaChat/',
    output_path = '/Users/troks27/Desktop/TrokaChat All for Publication/TrokaChat/TrokaChat/',
    files_list = ['H DEGs', 'NL DEGs', 'LS DEGs','nulldistavgs_H DEGs','nulldistavgs_NL DEGs','nulldistavgs_LS DEGs'],
    species = 'human',
    sample_name_from_R = 'sample',
    counts_file = 'counts',
    min_clus_percent = 0
)
#TrokaObject.process_files()


