
from InOut import *
from all import *


# programme principal
if __name__ == '__main__':
    genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")

    listes = []
    for i in genome_objects:
        integer_part = i.split('.')[0]
        listes.append(int(integer_part))
    y= set(listes)
    ly = list(y)

    genus_liste = methods.devide_requester(ly, 10)
    result = methods.unique_genus(genus_liste)
    print(len(result), result)