from InOut import Inout
from all import *

if __name__ == '__main__':
    """
    scripte retournant la liste des antibiotiques testés, par base de données et par classes.
    """
    argannot = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/argannot/sequences"#argannot ok
    card = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/card/sequences" #card ok
    ncbi = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/ncbi/sequences" #ecoh ok
    resfinder = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/resfinder/sequences" #card ok

    liste_antibio = []
    header_argannot = methods.get_header(argannot)
    x = methods.get_anti(header_argannot, "argannot")
    y = sorted(x)
    liste_antibio.append(y)
    header_card = methods.get_header(card)
    x = methods.get_anti(header_card, "card")
    y = sorted(x)
    liste_antibio.append(y)
    header_ncbi = methods.get_header(ncbi)
    x = methods.get_anti(header_ncbi, "ecoh")
    y = sorted(x)
    liste_antibio.append(y)
    header_resfinder = methods.get_header(resfinder)
    x = methods.get_anti(header_resfinder, "card")
    liste_antibio.append(x)

    name = ["argannot", "card", "ncbi", "resfinder"]

    liste_dico = []
    liste_dico_class = []
    e = 0
    for i in liste_antibio:

        z = methods.normalize_antibiotic_names(i)
        f = methods.get_anti_class(z)
        g = methods.present_elements(z)
        h = methods.present_elements(f)
        n = methods.count_element(z)

        liste_dico.append(g)
        liste_dico_class.append(h)

    dico_argannot = liste_dico[0]
    dico_card = liste_dico[1]
    dico_ncbi = liste_dico[2]
    dico_resfinder = liste_dico[3]

    Inout.create_csv_from_dicts(dico_argannot, dico_card, dico_ncbi, dico_resfinder, "antibiotiques_db.csv")

    dico_argannot = liste_dico_class[0]
    dico_card = liste_dico_class[1]
    dico_ncbi = liste_dico_class[2]
    dico_resfinder = liste_dico_class[3]

    Inout.create_csv_from_dicts(dico_argannot, dico_card, dico_ncbi, dico_resfinder, "class_db.csv")

    """
    Récuperer les genomes les plus testé par genus.
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