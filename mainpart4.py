from InOut import Inout
from all import *

if __name__ == '__main__':
    """
    parser resultats Abricate et PATRIC, donne un dossier de resultats globaux.
    """
    #lire les DB abricate
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

    e = 0

    for antibiotique in liste_antibio:
        liste_antibio_2 = []
        noms_normalises = methods.normalize_antibiotic_names(antibiotique)
        for nom in noms_normalises:
            if nom not in liste_antibio_2:
                liste_antibio_2.append(nom)
        print(name[e])
        print(liste_antibio_2)

        m = methods.parse_argannot_results("./data/abricate.txt")
        z = methods.info_from_results(m)
        genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")
        genus, ids = methods.get_selected_genome("genomelist.txt")
        count = 0
        total_count = []
        for i in range(len(ids)):
            liste_anti = []
            simple_list = []
            gid = ids[i]
            count_m = 0
            with open(f"./data/genome/{name[e]}/{genus[i]}.txt", 'w') as file:
                file.write(f"{genus[i]}\tPATRIC\tABRICATE\trelax\tstricte\n")
            antibiotic_data = genome_objects[gid].get_antibiotics_data()
            for antibiotics in antibiotic_data:
                dico_anti = {
                    "antibiotic": '',
                    "resistance": '',
                    "abricate": 'NOT',
                    "relax": "NC",
                    "stricte": "NC"
                }
                antibiotics["antibiotic"] = antibiotics["antibiotic"].replace("/", "+")
                antibiotic = antibiotics["antibiotic"]
                resistance = antibiotics["resistant_phenotype"]
                for dico in z:
                    if dico["file"] == genus[i]:
                        if dico["database"] == name[e]:
                            list_tested = dico["product_resistance"]
                            if antibiotic in list_tested:
                                dico_anti["abricate"] = "R"
                            if antibiotic not in liste_antibio_2:
                                if dico_anti["abricate"] != "ND":
                                    dico_anti["abricate"] = "ND"
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
                                if resistance != "ND" and dico_resi == "ND":
                                    dico["resistance"] = resistance
                                elif resistance == "ND" and dico_resi != "ND":
                                    dico["resistance"] = dico_resi
                                else:
                                    dico["resistance"] = "ND"
            for dico_anti in liste_anti:
                        if dico_anti["resistance"] == "R" and dico_anti["abricate"] == "R":
                            dico_anti["relax"] = "TP"
                            dico_anti["stricte"] = "TP"
                        if dico_anti["resistance"] == "R" and dico_anti["abricate"] == "NOT":
                            dico_anti["relax"] = "FN"
                            dico_anti["stricte"] = "FN"
                        if dico_anti["resistance"] == "R" and dico_anti["abricate"] == "ND":
                            dico_anti["relax"] = "NC"
                            dico_anti["stricte"] = "FN"
                        if dico_anti["resistance"] == "S" and dico_anti["abricate"] == "R":
                            dico_anti["relax"] = "FP"
                            dico_anti["stricte"] = "FP"
                        if dico_anti["resistance"] == "S" and dico_anti["abricate"] == "NOT":
                            dico_anti["relax"] = "TN"
                            dico_anti["stricte"] = "TN"
                        if dico_anti["resistance"] == "S" and dico_anti["abricate"] == "ND":
                            dico_anti["relax"] = "NC"
                            dico_anti["stricte"] = "TN"


            with open(f"./data/genome/{name[e]}/{genus[i]}.txt", 'a') as file:
                for x in liste_anti:
                    file.write(f'{x["antibiotic"]}\t{x["resistance"]}\t{x["abricate"]}\t{x["relax"]}\t{x["stricte"]}\n')
        e+=1
