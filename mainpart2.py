from InOut import Inout
from all import *

if __name__ == '__main__':
    """
    Scripte permettant de recuperer la liste des genus.
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
    liste_genus = []
    for i in genome_objects:
        genus = genome_objects[i].get_genus()
        if genus is not None:
            liste_genus.append(genus)

    unique_genus = set(liste_genus)
    print(sorted(list(unique_genus)))
    print(len(list(unique_genus)))