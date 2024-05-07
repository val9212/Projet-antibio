from InOut import Inout
from all import *

if __name__ == '__main__':
    """
    Scripte permettant la lecture du fichier de données Patric GENOME_AMR.txt
    """
    genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")
    Inout.to_table(genome_objects, antibiotics, "table.csv")