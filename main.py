from InOut import *

genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")


print(genome_objects['32002.4'].get_id())
print(genome_objects['32002.4'].get_name())
x = genome_objects['32002.4'].get_antibiotics_data()
for y in x:
    print(y)

Inout.to_table(genome_objects, antibiotics)