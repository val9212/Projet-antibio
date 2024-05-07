from InOut import Inout
from all import *

if __name__ == '__main__':
    """
    Scripte permettant la lecture du fichier de données Patric GENOME_AMR.txt
    """
    genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")
    listes = []
    for i in genome_objects:
        integer_part = i.split('.')[0]
        listes.append(int(integer_part))
    y = set(listes)
    ly = list(y)
    dico_list = methods.devide_requesterV2(ly)

    for i in genome_objects:
        if genome_objects[i].get_taxon_id() in dico_list:
            genome_objects[i].add_genus(dico_list[genome_objects[i].get_taxon_id()])
    methods.get_more_tested_genome(genome_objects, True)
    Inout.to_table(genome_objects, antibiotics, "table.csv")