import csv
from Gene import *

genome_id_set = set()
antibiotics_set = set()
with open('C:/Users/val92/PycharmProjects/Projet-antibio/data/genome_AMR.csv', mode='r') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    raw_data = [line for line in reader]

genome_objects = {}

for record in raw_data:

    genome_id_set.add(record['genome_id'])
    antibiotics_set.add(record['antibiotic'])


    gid = record['genome_id']
    if gid not in genome_objects:
        genome_objects[gid] = GenomeData(
            genome_id=gid,
            genome_name=record['genome_name']
        )


    # Ajoutez les informations spécifiques à l'antibiotique à l'objet.
    genome_objects[gid].add_antibiotic_data(
        antibiotic=record['antibiotic'],
        measurement_value=record['measurement_value'],
        measurement_unit=record['measurement_unit'],
        resistant_phenotype=record['resistant_phenotype'] if 'resistant_phenotype' in record else None,
        laboratory_typing_method=record['laboratory_typing_method']
    )

def to_table(genome_objects, antibiotics):
    # Les en-têtes seront 'genome_id', 'genome_name', suivis des antibiotiques
    headers = ['genome_id', 'genome_name'] + antibiotics

    # Initialiser un dictionnaire pour stocker les données
    data = {gid: {antibiotic: 'NA' for antibiotic in antibiotics} for gid in genome_objects}

    # Remplir le dictionnaire avec les valeurs de résistance
    for gid in genome_objects:
        data[gid]['genome_name'] = genome_objects[gid].get_name()
        info_antibio = genome_objects[gid].get_antibiotics_data()
        for antibiotic_data in info_antibio:
            antibiotic = antibiotic_data['antibiotic']
            if antibiotic in antibiotics:
                    # Utiliser le phenotype de résistance
                if antibiotic_data['resistant_phenotype'] == "":
                    data[gid][antibiotic] = "DATA_ERROR"
                else:
                    data[gid][antibiotic] = antibiotic_data['resistant_phenotype']

    # Écrire les données dans un fichier CSV
    with open("C:/Users/val92/PycharmProjects/Projet-antibio/table.csv", mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=headers)
        writer.writeheader()
        for gid in genome_objects:
            if gid in data:
                writer.writerow({'genome_id': gid, 'genome_name': data[gid]['genome_name'], **data[gid]})


print(genome_objects['32002.4'].get_id())
print(genome_objects['32002.4'].get_name())
x = genome_objects['32002.4'].get_antibiotics_data()
for y in x:
    print(y)

to_table(genome_objects, list(antibiotics_set))


