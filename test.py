#genome
x = methods.parse_argannot_results("./data/ARGANNOT.txt")
z = methods.info_from_results(x)
genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")
genus, ids = methods.get_selected_genome("genomelist.txt")
for i in range(len(ids)):
    liste_anti = []
    simple_list = []
    gid = ids[i]
    with open(f"./data/genome/{genus[i]}.txt", 'w') as file:
        file.write(f"{genus[i]}\tPATRIC\tABRICATE\n")
    antibiotic_data = genome_objects[gid].get_antibiotics_data()
    print(antibiotic_data)
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
            print(dico)
            if dico["file"] == genus[i]:
                list_tested = dico["product_resistance"]
                if antibiotic in list_tested:
                    dico_anti["abricate"] = "R"
                if antibiotic not in liste_antibio_2:
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
                        resistance = 4
                    elif resistance == 'Susceptible' or resistance == "IS" or resistance == "Susceptible-dose dependent":
                        resistance = 2
                    elif resistance == 'Intermediate' or resistance == "Reduced Susceptibility":
                        resistance = 3
                    elif resistance == '':
                        resistance = 1
                    dico_resi = dico["resistance"]
                    if dico_resi == "R":
                        dico_resi = 4
                    elif dico_resi == "S":
                        dico_resi = 2
                    elif dico_resi == "I":
                        dico_resi = 3
                    else:
                        dico_resi = 1

                    if resistance > dico_resi:
                        dico["resistance"] = resistance

                    if dico["resistance"] == 4:
                        dico["resistance"] = "R"
                    elif dico["resistance"] == 2:
                        dico["resistance"] = "S"
                    elif dico["resistance"] == 3:
                        dico["resistance"] = "I"
                    elif dico["resistance"] == 1:
                        dico["resistance"] = "ND"

    with open(f"./data/genome/{genus[i]}.txt", 'a') as file:
        for x in liste_anti:
            file.write(f'{x["antibiotic"]}\t{x["resistance"]}\t{x["abricate"]}\n')


#genome 2
    x = methods.parse_argannot_results("./data/ARGANNOT.txt")
    z = methods.info_from_results(x)
    genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")
    genus, ids = methods.get_selected_genome("genomelist.txt")
    for i in range(len(ids)):
        liste_anti = []
        simple_list = []
        gid = ids[i]
        with open(f"./data/genome2/{genus[i]}.txt", 'w') as file:
            file.write(f"{genus[i]}\tPATRIC\tABRICATE\n")
        antibiotic_data = genome_objects[gid].get_antibiotics_data()
        print(antibiotic_data)
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
                print(dico)
                if dico["file"] == genus[i]:
                    list_tested = dico["product_resistance"]
                    if antibiotic in list_tested:
                        dico_anti["abricate"] = "R"
                    if antibiotic not in liste_antibio_2:
                        dico_anti["abricate"] = "ND"

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



        with open(f"./data/genome2/{genus[i]}.txt", 'a') as file:
            for x in liste_anti:
                file.write(f'{x["antibiotic"]}\t{x["resistance"]}\t{x["abricate"]}\n')

#genome 3
x = methods.parse_argannot_results("./data/ARGANNOT.txt")
    z = methods.info_from_results(x)
    genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")
    genus, ids = methods.get_selected_genome("genomelist.txt")
    for i in range(len(ids)):
        liste_anti = []
        simple_list = []
        gid = ids[i]
        with open(f"./data/genome3/{genus[i]}.txt", 'w') as file:
            file.write(f"{genus[i]}\tPATRIC\tABRICATE\n")
        antibiotic_data = genome_objects[gid].get_antibiotics_data()
        print(antibiotic_data)
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
                print(dico)
                if dico["file"] == genus[i]:
                    list_tested = dico["product_resistance"]
                    if antibiotic in list_tested:
                        dico_anti["abricate"] = "R"
                    if antibiotic not in liste_antibio_2:
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
                            dico["resistance"] = "ND"


        with open(f"./data/genome3/{genus[i]}.txt", 'a') as file:
            for x in liste_anti:
                file.write(f'{x["antibiotic"]}\t{x["resistance"]}\t{x["abricate"]}\n')