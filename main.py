
from InOut import *
from all import *


# programme principal
if __name__ == '__main__':


    """genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")

    listes = []
    for i in genome_objects:
        integer_part = i.split('.')[0]
        listes.append(int(integer_part))
    y= set(listes)
    ly = list(y)

    genus_liste = methods.devide_requester(ly, 10)
    result = methods.unique_genus(genus_liste)
    print(len(result), result)"""

    argannot = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/argannot/sequences"#argannot ok
    card = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/card/sequences" #card ok
    ecoh = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/ecoh/sequences" #echo ok
    ecoli_vf = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/ecoli_vf/sequences" #?
    megares = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/megares/sequences" #ecoh ok
    ncbi = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/ncbi/sequences" #ecoh ok
    plasmidfinder = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/plasmidfinder/sequences" #?
    resfinder = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/resfinder/sequences" #card ok
    vfdb = "C:/Users/Valentin/PycharmProjects/Projet-antibio/data/db/vfdb/sequences" #echo no


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
    liste = []
    for i in x:
        y = methods.trouver_mot_cle(i)
        if y is None:
            print(i)
        liste.append(y)
    y = sorted(liste)
    liste_antibio.append(y)
    header_resfinder = methods.get_header(resfinder)
    x = methods.get_anti(header_resfinder, "card")
    liste_antibio.append(x)

    liste_dico = []
    liste_dico_class = []
    for i in liste_antibio:

        z = methods.normalize_antibiotic_names(i)
        f = methods.get_anti_class(z)
        g = methods.present_elements(z)
        h = methods.present_elements(f)
        n = methods.count_element(z)
        count = 0
        for j in n:
            count += n[j]

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



