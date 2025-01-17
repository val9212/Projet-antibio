import csv
from GData import *
import pandas as pd
from tqdm import tqdm

class Inout():
    @classmethod
    def read_file(cls, file):
        """
        Lit les données à partir d'un fichier CSV et crée des objets GenomeData à partir de ces données.

        :param file: Le chemin vers le fichier CSV à lire.
        :return: Un dictionnaire d'objets GenomeData et une liste d'antibiotiques uniques.
        """
        genome_id_set = set()
        antibiotics_set = set()

        with open(file, mode='r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            raw_data = [line for line in reader]

        genome_objects = {}

        for record in tqdm(raw_data):

            genome_id_set.add(record['genome_id'])
            antibiotics_set.add(record['antibiotic'])


            gid = record['genome_id']
            if gid not in genome_objects:
                genome_objects[gid] = GenomeData(
                    genome_id=gid,
                    genome_name=record['genome_name'],
                    taxon_id=record['taxon_id']
                )


            # Ajoutez les informations spécifiques à l'antibiotique à l'objet.
            genome_objects[gid].add_antibiotic_data(
                antibiotic=record['antibiotic'],
                measurement_value=record['measurement_value'],
                measurement_unit=record['measurement_unit'],
                resistant_phenotype=record['resistant_phenotype'] if 'resistant_phenotype' in record else None,
                laboratory_typing_method=record['laboratory_typing_method']
            )
        return genome_objects, list(antibiotics_set)

    @classmethod
    def to_table(cls, genome_objects, antibiotics, path = "./results/table2.csv"):
        """
        Convertit les objets GenomeData en un tableau CSV contenant les données d'antibiotiques.

        :param genome_objects: Un dictionnaire d'objets GenomeData.
        :param antibiotics: Une liste d'antibiotiques uniques.
        :param path: Le chemin du fichier CSV de sortie.
        """
        headers = ['genome_id', 'genome_name', 'taxon_id'] + antibiotics

        data = {gid: {antibiotic: 'NA' for antibiotic in antibiotics} for gid in genome_objects}

        for gid in genome_objects:
            data[gid]['genome_name'] = genome_objects[gid].get_name()
            data[gid]['taxon_id'] = genome_objects[gid].get_taxon_id()
            info_antibio = genome_objects[gid].get_antibiotics_data()
            for antibiotic_data in info_antibio:
                antibiotic = antibiotic_data['antibiotic']
                if antibiotic in antibiotics:
                    if data[gid][antibiotic] == "NA" or data[gid][antibiotic] == "DATA_ERROR":
                        if antibiotic_data['resistant_phenotype'] == "":
                            data[gid][antibiotic] = "DATA_ERROR"
                        else:
                            data[gid][antibiotic] = antibiotic_data['resistant_phenotype']

        with open(path, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=headers)
            writer.writeheader()
            for gid in genome_objects:
                if gid in data:
                    writer.writerow({'genome_id': gid, 'genome_name': data[gid]['genome_name'], **data[gid]})


    @classmethod
    def write_list_to_file(cls,data_list, file_name='test.txt'):
        """
        Writes the elements of data_list to file_name, separated by commas.

        :param data_list: List of items to write to the file.
        :param file_name: The name of the file to create and write to.
        """
        try:
            # Open the file in write mode ('w') to create it or overwrite if it already exists
            with open(file_name, 'w') as file:
                # Join the list items into a single string separated by commas
                file.write(','.join(str(item) for item in data_list))
            print(f"Data successfully written to {file_name}")
        except IOError as e:
            # Handle the error and print an error message
            print(f"An error occurred: {e.strerror}")

    @classmethod
    def create_csv_from_dicts(cls, dict1, dict2, dict3, dict4, output_file):
        """
        Crée un fichier CSV à partir de quatre dictionnaires en fusionnant les données.

        :param dict1: Premier dictionnaire.
        :param dict2: Deuxième dictionnaire.
        :param dict3: Troisième dictionnaire.
        :param dict4: Quatrième dictionnaire.
        :param output_file: Le nom du fichier CSV de sortie.
        """
        df1 = pd.DataFrame(list(dict1.items()), columns=['Antibiotique', 'Argannot'])
        df2 = pd.DataFrame(list(dict2.items()), columns=['Antibiotique', 'card'])
        df3 = pd.DataFrame(list(dict3.items()), columns=['Antibiotique', 'ncbi'])
        df4 = pd.DataFrame(list(dict4.items()), columns=['Antibiotique', 'resfinder'])

        df = pd.merge(df1, df2, on='Antibiotique', how='outer')
        df = pd.merge(df, df3, on='Antibiotique', how='outer')
        df = pd.merge(df, df4, on='Antibiotique', how='outer')

        df.fillna(0, inplace=True)
        df.to_csv(output_file, index=False)

        print(f"Fichier CSV '{output_file}' créé avec succès !")


