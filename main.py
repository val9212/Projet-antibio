
from InOut import *
from all import *
from tqdm import tqdm
import numpy as np


# programme principal
if __name__ == '__main__':

    #lire les DB abricate
    argannot = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/argannot/sequences"#argannot ok
    card = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/card/sequences" #card ok
    ecoh = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/ecoh/sequences" #echo ok
    ecoli_vf = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/ecoli_vf/sequences" #?
    megares = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/megares/sequences" #ecoh ok
    ncbi = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/ncbi/sequences" #ecoh ok
    plasmidfinder = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/plasmidfinder/sequences" #?
    resfinder = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/resfinder/sequences" #card ok
    vfdb = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/vfdb/sequences" #echo no


    liste_antibio = []
    header_argannot = methods.get_header(argannot)
    x = methods.get_anti(header_argannot, "argannot")
    y =sorted(x)
    liste_antibio.append(y)
    header_card = methods.get_header(card)
    x = methods.get_anti(header_card, "card")
    y =sorted(x)
    liste_antibio.append(y)
    header_ncbi = methods.get_header(ncbi)
    x = methods.get_anti(header_ncbi, "ecoh")
    y = sorted(x)
    liste_antibio.append(y)
    header_resfinder = methods.get_header(resfinder)
    x = methods.get_anti(header_resfinder, "card")
    liste_antibio.append(x)

    name = ["argannot", "card", "ncbi", "resfinder"]
    liste_antibio_2 = []
    liste_dico_class = []
    e = 0

    for antibiotique in liste_antibio:
        noms_normalises = methods.normalize_antibiotic_names(antibiotique)
        for nom in noms_normalises:
            if nom not in liste_antibio_2:
                liste_antibio_2.append(nom)

    x = methods.parse_argannot_results("./data/ARGANNOT.txt")
    z = methods.info_from_results(x)
    genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")
    genus, ids = methods.get_selected_genome("genomelist.txt")
    count = 0
    total_count = []
    for i in range(len(ids)):
        liste_anti = []
        simple_list = []
        gid = ids[i]
        count_m = 0
        with open(f"./data/genome3/{genus[i]}.txt", 'w') as file:
            file.write(f"{genus[i]}\tPATRIC\tABRICATE\n")
        antibiotic_data = genome_objects[gid].get_antibiotics_data()
        for antibiotics in antibiotic_data:
            dico_anti = {
                "antibiotic": '',
                "resistance": '',
                "abricate": 'NOT'
            }
            antibiotics["antibiotic"] = antibiotics["antibiotic"].replace("/", "+")
            antibiotic = antibiotics["antibiotic"]
            resistance = antibiotics["resistant_phenotype"]
            for dico in z:
                if dico["file"] == genus[i]:
                    list_tested = dico["product_resistance"]
                    if antibiotic in list_tested:
                        dico_anti["abricate"] = "R"
                    if antibiotic not in liste_antibio_2:
                        if dico_anti["abricate"] != "ND":
                            dico_anti["abricate"] = "ND"
                            count +=1
                            count_m +=1
            if antibiotic not in simple_list:
                simple_list.append(antibiotic)
                dico_anti["antibiotic"] = antibiotic
                if resistance == 'Resistant' or resistance == "Nonsusceptible":
                    dico_anti["resistance"] = 'R'
                elif resistance == 'Susceptible' or resistance == "IS" or resistance == "Susceptible-dose dependent":
                    dico_anti["resistance"] = "S"
                elif resistance == 'Intermediate' or resistance == "Reduced Susceptibility":
                    dico_anti["resistance"] = "I"
                elif resistance == '':
                    dico_anti["resistance"] = "ND"
                liste_anti.append(dico_anti)
            else:
                for dico in liste_anti:
                    if dico["antibiotic"] == antibiotic:
                        if resistance == 'Resistant' or resistance == "Nonsusceptible":
                            resistance = 'R'
                        elif resistance == 'Susceptible' or resistance == "IS" or resistance == "Susceptible-dose dependent":
                            resistance = 'S'
                        elif resistance == 'Intermediate' or resistance == "Reduced Susceptibility":
                            resistance = 'I'
                        elif resistance == '':
                            resistance = "ND"
                        dico_resi = dico["resistance"]

                        if resistance != dico_resi:
                            dico["resistance"] = "ND"
        print(genus[i], count_m)
        total_count.append(count_m)

        with open(f"./data/genome3/{genus[i]}.txt", 'a') as file:
            for x in liste_anti:
                file.write(f'{x["antibiotic"]}\t{x["resistance"]}\t{x["abricate"]}\n')
        print(count)
        print(sum(total_count)/len(total_count))
        print(len(total_count))
        print(np.median(total_count))
        print(max(total_count))


        # lire la DB de PATRICK
    """
    genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")
    listes = []
    for i in genome_objects:
        integer_part = i.split('.')[0]
        listes.append(int(integer_part))
    y = set(listes)
    ly = list(y)
    dico_list = methods.devide_requesterV2(ly)

    for i in tqdm(genome_objects):
        if genome_objects[i].get_taxon_id() in dico_list:
            genome_objects[i].add_genus(dico_list[genome_objects[i].get_taxon_id()])
    methods.get_more_tested_genome(genome_objects, True)
    Inout.to_table(genome_objects, antibiotics)"""

    """listes = [] 
    # ARBRE phylogenetic (transformation des taxons ID en GENUS)
    for i in genome_objects:
        integer_part = i.split('.')[0]
        listes.append(int(integer_part))
    y= set(listes)
    ly = list(y)

    genus_liste = methods.devide_requester(ly, 10)
    result = methods.unique_genus(genus_liste)
    print(len(result), result)"""

    """
    #lire les DB abricate
    argannot = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/argannot/sequences"#argannot ok
    card = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/card/sequences" #card ok
    ecoh = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/ecoh/sequences" #echo ok
    ecoli_vf = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/ecoli_vf/sequences" #?
    megares = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/megares/sequences" #ecoh ok
    ncbi = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/ncbi/sequences" #ecoh ok
    plasmidfinder = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/plasmidfinder/sequences" #?
    resfinder = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/resfinder/sequences" #card ok
    vfdb = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/vfdb/sequences" #echo no


    liste_antibio = []
    header_argannot = methods.get_header(argannot)
    x = methods.get_anti(header_argannot, "argannot")
    y =sorted(x)
    liste_antibio.append(y)
    header_card = methods.get_header(card)
    x = methods.get_anti(header_card, "card")
    y =sorted(x)
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
        print(f"{name[e]}\n{n}'")

    dico_argannot = liste_dico[0]
    dico_card = liste_dico[1]
    dico_ncbi = liste_dico[2]
    dico_resfinder = liste_dico[3]

    Inout.create_csv_from_dicts(dico_argannot, dico_card, dico_ncbi, dico_resfinder, "antibiotiquesV2_db.csv")

    dico_argannot = liste_dico_class[0]
    dico_card = liste_dico_class[1]
    dico_ncbi = liste_dico_class[2]
    dico_resfinder = liste_dico_class[3]

    Inout.create_csv_from_dicts(dico_argannot, dico_card, dico_ncbi, dico_resfinder, "classV2_db.csv")

    #faire un tableau qui compte le nombre de genome testé et son resultat
    'methods.count_table("./results/table2.csv")'"""



