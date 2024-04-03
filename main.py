
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

    argannot = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/argannot/sequences"
    card = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/card/sequences"
    ecoh = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/ecoh/sequences"
    ecoli_vf = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/ecoli_vf/sequences"
    megares = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/megares/sequences"
    ncbi = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/ncbi/sequences"
    plasmidfinder = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/pasmidfinder/sequences"
    resfinder = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/resfinder/sequences"
    vfdb = "C:/Users/val92/PycharmProjects/Projet-antibio/data/db/vfdb/sequences"

    header_argannot = methods.get_header(argannot)
    print(header_argannot)
