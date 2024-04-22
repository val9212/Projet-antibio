from requester import *
import re
import csv
import os
from Bio.Data import IUPACData
import pandas as pd

class methods():

    @classmethod
    def unique_genus(cls, genus_list):
        final_list = []
        for sublist in genus_list:
            for item in sublist:
                if item not in final_list:
                    final_list.append(item)
        return final_list

    @classmethod
    def devide_requester(cls, taxon_id, n_chunk = 10):
        i = 0
        genus_list = []
        n = int(len(taxon_id) / n_chunk)
        while i < len(taxon_id):
            w = request.requester(taxon_id[i:n + i])
            f = request.extract_genus(w)
            genus_list.append(f)
            i += n
        return genus_list

    @classmethod
    def devide_requesterV2(cls, taxon_id, n_chunk = 10):
        genus_dict = {}

        chunk_size = int(len(taxon_id) / n_chunk)
        for i in range(0, len(taxon_id), chunk_size):
            chunk_taxon_ids = taxon_id[i:i + chunk_size]  # Sélectionner le chunk actuel d'ID de taxon
            response = request.requester(chunk_taxon_ids)
            genus_chunk = request.extract_genus(response)
            for taxon, genus in zip(chunk_taxon_ids, genus_chunk):
                genus_dict[taxon] = genus

        return genus_dict

    @classmethod
    def get_header(cls, file):
        headers = []
        with open(file, 'r') as fasta_file:
            for line in fasta_file:
                if line.startswith('>'):
                    headers.append(line.strip())

        return headers

    @classmethod
    def replace_semicolon(cls, chaine):
        nchaine = chaine.replace('~~~', '~~~ ')
        return nchaine.replace(';', ' ')

    @classmethod
    def get_anti(cls, headers, param):
        antibiotic_classes = []
        if param == "argannot":
            antibiotic_class_regex = re.compile(r'>[^~]+~~~\(([^)]+)\)')
        if param == "card":
            antibiotic_class_regex = re.compile(r'~~~(\w+)(?:\s|$)')
        if param == "ecoh":
            antibiotic_class_regex = re.compile(r'~~~\S+~~~\S+\s(.+)')

        for header in headers:
            match = antibiotic_class_regex.search(header)
            if match:
                antibiotic_class = match.group(1)
                antibiotic_classes.append(antibiotic_class)
            else:
                headerspace = cls.replace_semicolon(header)
                anti = cls.trouver_mot_cle(headerspace)
                if anti is not None:
                    antibiotic_classes.append(anti)

        return antibiotic_classes

    @classmethod
    def normalize_antibiotic_names_s(cls, name):
        antibiotic_normalization_dict = {
            "amp": "ampicillin",
            "arma": "aminoglycoside",
            "arr": "rifamycin",
            "amox": "amoxicillin",
            "agly": "aminoglycoside",
            "agly_flqn": "aminoglycoside",
            "fcd": "fidaxomicin",
            "fcyn": "cephalosporin",
            "cep": "cephalexin",
            "flq": "fluoroquinole",
            "cfx": "cefuroxime",
            "gly": "glycopeptide",
            "mls": "macrolide",
            "phe": "phenylanin",
            "sul": "sulfamide",
            "col": "colistin",
            "tmt": "trimethoprim",
            "imi": "imipenem",
            "cfi": "imipenem",
            "cmr": "chloramphenicol",
            "pom": "beta-lactamas",
            "cph": "capreomycin",
            "fus": "fusidic_acid",
            "fusidic": "fusidic_acid",
            "pen": "penam",
            "TetracenomycinC": "tetracycline C",
            "amoxil": "amoxicillin",
            "penam": "penam",
            "antibacterial_free_fatty_acids": "antibacterial_free_fatty_acids",
            "16s": "16s",
            "23s": "23s",
            "beta-lactam": "beta-lactamas",
            "beta-lactamase": "beta-lactamas",
            "metallo-beta-lactamase": "beta-lactamas",
            "erm": "erythromycin",
            "fos": "fosfomycin",
            "fomb": "fosfomycin",
            "multidrug": "multiple antibiotics",
            "mupirocin-resistant": "mupirocin",
            "nim": "nitroimidazole",
            "rifamycin-inactivating": "rifamycin",
            "aac": "aminoglycoside",
            "aph": "aminoglycoside",
            "ceftaroline-resistant": "ceftaroline",
            "amg": "aminoglycosides",
            "amk": "amikacin",
            "dfr": "diaminopyrimidique",
            "amu": "aminocoumarin",
            "mcr": "colistin",
            "mdf": "tetracycline",
            "amx": "amoxicillin",
            "atm": "aztreonam",
            "avi": "avibactam",
            "azm": "azithromycin",
            "bdq": "bedaquiline",
            "bla": "beta-lactams",
            "cap": "capreomycin",
            "cef": "ceftazidime",
            "cza": "ceftazidime-Avibactam",
            "chl": "chloramphenicol",
            "cip": "ciprofloxacin",
            "cli": "clindamycin",
            "clr": "clarithromycin",
            "cst": "colistin",
            "dao": "dapsone",
            "dap": "daptomycin",
            "dcs": "d-cycloserine",
            "edn": "edeine",
            "elf": "elfamycin",
            "emb": "ethambutol",
            "emcm": "ethambutol/Capreomycin",
            "enc": "enacyloxin IIa",
            "enr": "enrofloxacin",
            "ery": "erythromycin",
            "eto": "ethionamide",
            "fa": "fusidic acid",
            "flo": "fluoroquinolones",
            "fof": "fosfomycin",
            "g418": "g418",
            "ge2a": "gE2270A",
            "gen": "gentamicin",
            "genc": "gentamicin C",
            "hgm": "hyrgomycin B",
            "inh": "isoniazid",
            "imp": "imipenem",
            "kan": "kanamycin",
            "kas": "kasugamicin",
            "kir": "kirromycin",
            "lys": "lysocin (E)",
            "lzd": "linezolid",
            "ntmdz": "nitroimidazole",
            "oxzln": "oxazolidinone",
            "mac": "macrolide",
            "mef": "macrolide",
            "mph": "macrolide",
            "msr": "macrolide",
            "mult": "multiple antibiotics",
            "mup": "mupirocin",
            "mtz": "metronidazole",
            "mxf": "moxifloxacin",
            "crp": "ciprofloxacine",
            "ros": "fluoroquinolones",
            "A": "SPIRAMYCINE/METRONIDAZOLE",
            "neo": "neomycin",
            "nit": "nitrofurantoin",
            "oxz": "oxazolidinone",
            "pac": "pactamycin",
            "par": "paromomycin",
            "pas": "para-aminosalicylic acid",
            "pcl": "perchlozone",
            "plm": "pleuromutilin",
            "plv": "pulvomycin",
            "pto": "prothionamide",
            "pza": "pyrazinamide",
            "rfb": "rifabutin",
            "qnr": "fluoroquinolone",
            "rmt": "aminoglycoside",
            "sal": "streptogramin",
            "rif": "rifampicin",
            "slf": "sulfonamides",
            "spt": "spectinomycin",
            "str": "streptomycin",
            "tmp": "trimethoprim",
            "tet": "tetracycline",
            "tob": "tobramycin",
            "trc": "triclosan",
            "tyl": "tylosin",
            "van": "vancomycin",

            "vio": "viomycin",
            "zol": "zoliflodacin",
        }
        name_lower = name.lower()
        if name_lower in antibiotic_normalization_dict:
            name_lower = antibiotic_normalization_dict[name_lower]

        return name_lower

    @classmethod
    def normalize_antibiotic_names(cls, antibiotic_list):
        antibiotic_normalization_dict = {
            "amp": "ampicillin",
            "arma": "aminoglycoside",
            "arr": "rifamycin",
            "amox": "amoxicillin",
            "agly": "aminoglycoside",
            "agly_flqn": "aminoglycoside",
            "bla": "beta-lactamas",
            "fcd": "fidaxomicin",
            "fcyn": "cephalosporin",
            "cef": "cefepime",
            "cep": "cephalexin",
            "flq": "fluoroquinole",
            "cfx": "cefuroxime",
            "gly": "glycopeptide",
            "mls": "macrolide",
            "phe": "phenylanin",
            "rif": "rifampicin",
            "sul": "sulfamide",
            "tet": "tetracycline",
            "col": "colistin",
            "tmt": "trimethoprim",
            "imi": "imipenem",
            "cfi": "imipenem",
            "cmr": "chloramphenicol",
            "pom": "beta-lactamas",
            "cph": "capreomycin",
            "fus": "fusidic_acid",
            "fusidic": "fusidic_acid",
            "pen": "penam",
            "TetracenomycinC": "tetracycline C",
            "amoxil": "amoxicillin",
            "penam": "penam",
            "antibacterial_free_fatty_acids": "antibacterial_free_fatty_acids",
            "16s": "16s",
            "23s": "23s",
            "beta-lactam": "beta-lactamas",
            "beta-lactamase": "beta-lactamas",
            "metallo-beta-lactamase": "beta-lactamas",
            "erm": "erythromycin",
            "fos": "fosfomycin",
            "fomb": "fosfomycin",
            "multidrug": "multiple antibiotics",
            "mupirocin-resistant": "mupirocin",
            "nim": "nitroimidazole",
            "rifamycin-inactivating": "rifamycin",
            "aac": "aminoglycoside",
            "aph": "aminoglycoside",
            "ceftaroline-resistant": "ceftaroline",
            "AMG": "aminoglycosides",
            "AMK": "amikacin",
            "dfr": "diaminopyrimidique",
            "AMU": "aminocoumarin",
            "mcr": "colistin",
            "mdf": "tetracycline",
            "AMX": "amoxicillin",
            "ATM": "aztreonam",
            "AVI": "avibactam",
            "AZM": "azithromycin",
            "BDQ": "bedaquiline",
            "BLA": "beta-lactams",
            "CAP": "capreomycin",
            "CEF": "ceftazidime",
            "CZA": "ceftazidime-Avibactam",
            "CHL": "chloramphenicol",
            "CIP": "ciprofloxacin",
            "CLI": "clindamycin",
            "CLR": "clarithromycin",
            "CST": "colistin",
            "DAO": "dapsone",
            "DAP": "daptomycin",
            "DCS": "d-cycloserine",
            "EDN": "edeine",
            "ELF": "elfamycin",
            "EMB": "ethambutol",
            "EMCM": "ethambutol & Capreomycin",
            "ENC": "enacyloxin IIa",
            "ENR": "enrofloxacin",
            "ERY": "erythromycin",
            "ETO": "ethionamide",
            "FA": "fusidic acid",
            "FLO": "fluoroquinolones",
            "FOF": "fosfomycin",
            "G418": "g418",
            "GE2A": "gE2270A",
            "GEN": "gentamicin",
            "GENC": "gentamicin C",
            "HGM": "hyrgomycin B",
            "INH": "isoniazid",
            "IPM": "imipenem",
            "KAN": "kanamycin",
            "KAS": "kasugamicin",
            "KIR": "kirromycin",
            "LYS": "lysocin (E)",
            "LZD": "linezolid",
            "ntmdz": "nitroimidazole",
            "oxzln": "oxazolidinone",
            "MAC": "macrolide",
            "mef": "macrolide",
            "mph": "macrolide",
            "msr": "macrolide",
            "MULT": "multiple antibiotics",
            "MUP": "mupirocin",
            "MTZ": "metronidazole",
            "MXF": "moxifloxacin",
            "NEO": "neomycin",
            "NIT": "nitrofurantoin",
            "OXZ": "oxazolidinone",
            "PAC": "pactamycin",
            "PAR": "paromomycin",
            "PAS": "para-aminosalicylic acid",
            "PCL": "perchlozone",
            "PLM": "pleuromutilin",
            "PLV": "pulvomycin",
            "PTO": "prothionamide",
            "PZA": "pyrazinamide",
            "RFB": "rifabutin",
            "ros": "fluoroquinolones",
            "crp": "ciprofloxacine",
            "qnr": "fluoroquinolone",
            "rmt": "aminoglycoside",
            "sal": "streptogramin",
            "RIF": "rifampicin",
            "SLF": "sulfonamides",
            "SPT": "spectinomycin",
            "STR": "streptomycin",
            "TMP": "trimethoprim",
            "TET": "tetracycline",
            "TOB": "tobramycin",
            "TRC": "triclosan",
            "TYL": "tylosin",
            "van": "vancomycin",
            "VAN": "vancomycin",
            "VIO": "viomycin",
            "ZOL": "zoliflodacin",
        }
        normalized_list = []
        for name in antibiotic_list:
            # Convertir en minuscules
            name_lower = name.lower()
            # Remplacer par le nom normalisé s'il existe dans le dictionnaire
            if name_lower in antibiotic_normalization_dict:
                name_lower = antibiotic_normalization_dict[name_lower]

            normalized_list.append(name_lower)
        return normalized_list

    @classmethod
    def normalize_antibiotic_names2(cls, name):
        normalized_name = IUPACData.ambiguous_dna_values.get(name, name)
        return normalized_name

    @classmethod
    def trouver_mot_cle(cls, chaine):
        mots_cles = ["amoxicillin", "Amikacin", "Clindamycin", "Cefepime", "Spiramycin", "Piperacillin+Tazobactam", "Ceftazidime", "Aztreonam", "Imipenem", "Tetracycline", "Gentamicin", "Dalfopristin", "doxycycline", "carbapenem", "Quinupristin", "Lincomycin", "rifamycin", "peptide", "cephalosporin", "D-lactate", "streptogramin", "mupirocin-resistant", "zorbamycin", "fomb", "fos", "phenicol", "polymyxin", "teicoplanin", "VanY-like", "tetracenomycin", "ciprofloxacin", "flavin", "streptomycin", "nitroimidazole", "ceftaroline-resistant", "aminocyclitol", "penicillinase", "fusidic", "peptidoglycan", "rifamycin-inactivating", "viomycin", "lincosamide", "beta-lactamase", "bleomycin", "kasugamycin", "aminoglycoside", "arma", "rifampin", "chloramphenicol", "glycopeptide", "trimethoprim", "erythromycin", "23s", "fosfomycin", "abc-f", "phosphoethanolamine--lipid", "beta-lactam", "macrolide", "oleandomycin", "multidrug", "quinolone", "serine", "d-ala-d-ala", "vancomycin", "alanine", "tetracycline", "sulfonamide", "streptothricin", "16s", "fluoroquinolone", "multidrug", "beta-lactam-resistant" , "methicillin", "macrolide", "metallo-beta-lactamase" ]
        for mot_cle in mots_cles:
            if mot_cle.lower() in chaine.lower():
                return mot_cle

        return cls.extraire_trois_caracteres_apres_tilde(chaine.lower())

    @classmethod
    def count_element(cls, liste):
        occurrences_dict = {}
        for element in liste:
            if element in occurrences_dict:
                occurrences_dict[element] += 1
            else:
                occurrences_dict[element] = 1
        return occurrences_dict

    @classmethod
    def present_elements(cls, liste):
        dictio = {}
        for element in liste:
            if element not in dictio:
                dictio[element] = 1

        return dictio

    @classmethod
    def extraire_trois_caracteres_apres_tilde(cls, chaine):
        # Trouver l'index du dernier "~~~"
        index_dernier_tilde = chaine.rfind("~~~")

        # Vérifier si le dernier "~~~" est présent dans la chaîne
        if index_dernier_tilde != -1:
            # Extraire les 3 caractères suivant l'index du dernier "~~~"
            trois_caracteres_apres_tilde = chaine[index_dernier_tilde + 5:index_dernier_tilde + 8]
            return trois_caracteres_apres_tilde
        else:
            # Si le dernier "~~~" n'est pas trouvé, retourner une chaîne vide
            return chaine

    @classmethod
    def get_anti_class(cls, normalized_list):
        class_list = []
        anti_class = {
                "16s": "ND",
                "23s": "ND",
                "abc-f": "ND",
                "acridine_dye": "ND",
                "alanine": "ND",
                "amikacin": "Aminoglycoside",
                "aminocoumarin": "ND",
                "aminocyclitol": "ND",
                "aminoglycoside": "Aminoglycoside",
                "amoxicillin": "Pénicilline",
                "ampicillin": "Pénicilline",
                "ant": "ND",
                "antibacterial_free_fatty_acids": "ND",
                "aztreonam": "aztreonam",
                "beta-lactamas": "ND",
                "bicyclomycin": "ND",
                "bleomycin": "ND",
                "capreomycin": "ND",
                "carbapenem": "Carbapénème",
                "carbomycin": "ND",
                "cefepime": "Céphalosporine",
                "ceftaroline": "ND",
                "ceftazidime": "Céphalosporine",
                "cefuroxime": "ND",
                "cephalexin": "Cephalosporine",
                "cephalosporin": "Céphalosporine",
                "cephamycin": "ND",
                "chloramphenicol": "ND",
                "ciprofloxacin": "Fluoroquinolone",
                "clindamycin": "Clindamycin",
                "colistin": "Polypeptides",
                "d-ala-d-ala": "ND",
                "d-lactate": "ND",
                "dalfopristin": "Quinuprsitin+Dalfospristin",
                "diaminopyrimidine": "ND",
                "diaminopyrimidique": "ND",
                "elfamycin": "ND",
                "erythromycin": "Macrolide",
                "fidaxomicin": "Macrolide",
                "flavin": "ND",
                "fluoroquinole": "Fluoroquinolone",
                "fluoroquinolone": "Fluoroquinolone",
                "fosfomycin": "Fosfomycin",
                "fusidic_acid": "ND",
                "gentamicin": "Aminoglycoside",
                "glycopeptide": "Glycopeptide",
                "hug": "ND",
                "imipenem": "Carbapénème",
                "kasugamycin": "ND",
                "lincomycin": "Lincosamide",
                "lincosamide": "Lincosamide",
                "macrolide": "Macrolide",
                "macrolides": "Macrolide",
                "methicillin": "Pénicilline",
                "multiple antibiotics": "ND",
                "mupirocin": "Mupirocin",
                "nitroimidazole": "ND",
                "nucleoside": "ND",
                "oleandomycin": "ND",
                "oxazolidinone": "Oxazolidinone",
                "penam": "Penicilline",
                "penicillinase": "Penicilline",
                "peptide": "ND",
                "peptidoglycan": "ND",
                "phenicol": "Phénicol",
                "phenylanin": "ND",
                "phosphoethanolamine--lipid": "ND",
                "piperacillin+tazobactam": "Penicilline",
                "pleuromutilin": "ND",
                "quinolone": "ND",
                "quinupristin": "Quinuprsitin+Dalfospristin",
                "rifampicin": "Rifamycine",
                "rifampin": "Rifamycin",
                "rifamycin": "Rifamycin",
                "serine": "ND",
                "spiramycin": "Macrolide",
                "streptogramin": "ND",
                "streptomycin": "Aminoglycoside",
                "streptothricin": "ND",
                "sulfamethoxazole": "Sulfamide",
                "sulfamide": "Sulfamide",
                "sulfonamide": "Sulfonamide",
                "teicoplanin": "Glycopeptide",
                "tetracenomycin": "ND",
                "tetracenomycinc": "ND",
                "tetracycline": "Tetracycline",
                "tiamulin": "ND",
                "tobramycin": "Aminoglycoside",
                "triclosan": "ND",
                "trimethoprim": "trimethoprim",
                "tylosin": "ND",
                "vancomycin": "Glycopeptide",
                "vany-like": "ND",
                "viomycin": "ND",
                "zorbamycin": "ND"
            }
        for name in normalized_list:
            if name in anti_class:
                name = anti_class[name]

                class_list.append(name)

        return class_list

    @classmethod
    def count_table(cls, fichier_csv):
        antibiotiques = {}

        with open(fichier_csv, 'r') as csv_file:
            lecteur_csv = csv.DictReader(csv_file)

            for ligne in lecteur_csv:
                for colonne, valeur in ligne.items():
                    if colonne != 'genome_id' and colonne != 'genome_name' and colonne != 'taxon_id':
                        if colonne in antibiotiques:
                            if valeur == 'Resistant' or valeur == "Nonsusceptible":
                                antibiotiques[colonne][0] += 1
                            elif valeur == 'Susceptible' or valeur == "IS" or valeur == "Susceptible-dose dependent":
                                antibiotiques[colonne][1] += 1
                            elif valeur == 'Intermediate' or valeur == "Reduced Susceptibility":
                                antibiotiques[colonne][2] += 1
                            elif valeur == 'DATA_ERROR':
                                antibiotiques[colonne][3] += 1

                        else:
                            if valeur == 'Resistant' or valeur == "Nonsusceptible":
                                antibiotiques[colonne] = [1, 0, 0, 0]
                            elif valeur == 'Susceptible' or valeur == "IS" or valeur == "Susceptible-dose dependent":
                                antibiotiques[colonne] = [0, 1, 0, 0]
                            elif valeur == 'Intermediate' or valeur == "Reduced Susceptibility":
                                antibiotiques[colonne] = [0, 0, 1, 0]
                            elif valeur == 'DATA_ERROR':
                                antibiotiques[colonne] = [0, 0, 0, 1]

        tableau_resultats = [["Antibiotique", "#R", "#S", "#I", "#DATA_ERROR"]]

        for antibiotique, compteurs in antibiotiques.items():
            tableau_resultats.append([antibiotique, str(compteurs[0]), str(compteurs[1]), str(compteurs[2]), str(compteurs[3])])

        with open("test.txt", 'w') as file:
            for data in tableau_resultats:
                file.write('\t'.join(data))
                file.write('\n')

    @classmethod
    def get_more_tested_genome(cls, genome_object, no_missing = False):
        genus_list = []
        genome_list = []
        nber_of_test_list = []
        for i in genome_object:
            genus = genome_object[i].get_genus()
            count = 0
            if genus not in genus_list and genus is not None:
                genus_list.append(genus)
                genome_list.append(genome_object[i].get_id())
                if no_missing:
                    antibiotics_data = genome_object[i].get_antibiotics_data()
                    for antibiotic in antibiotics_data:
                        if antibiotic["resistant_phenotype"] != '':
                            count += 1
                    nber_of_test_list.append(count)
                else:
                    nber_of_test_list.append(len(genome_object[i].get_antibiotics_data()))
            elif genus is not None:
                index = genus_list.index(genus)
                if no_missing:
                    antibiotics_data = genome_object[i].get_antibiotics_data()
                    for antibiotic in antibiotics_data:
                        if antibiotic["resistant_phenotype"] != '':
                            count += 1
                else:
                    count = len(genome_object[i].get_antibiotics_data())
                if nber_of_test_list[index] < count:
                    genome_list[index] = genome_object[i].get_id()
                    nber_of_test_list[index] = count


        with open("genomelist2.txt", 'w') as file:
            file.write(f"genus\tgenome_id\tnber_of_test\n")
            for i in range(len(genus_list)):
                file.write(f"{genus_list[i]}\t")
                file.write(f"{genome_list[i]}\t")
                file.write(f"{nber_of_test_list[i]}\t")
                file.write('\n')

    @classmethod
    def get_selected_genome(cls, file):
        liste_genus = []
        liste_id = []
        with open(file, 'r', newline='') as csvfile:
            # Créer un lecteur CSV
            lecteur_csv = csv.reader(csvfile, delimiter='\t')
            next(lecteur_csv)
            for ligne in lecteur_csv:
                # Extraire les données de la ligne
                genus = ligne[0]
                liste_genus.append(genus)
                genome_id = ligne[1]
                liste_id.append(genome_id)
        return liste_genus, liste_id
    @classmethod
    def remove_extension(cls,file_name):
        return os.path.splitext(file_name)[0]

    @classmethod
    def parse_argannot_results(cls, file_path):
        # Lecture du fichier CSV dans un DataFrame pandas
        resultats_df = pd.read_csv(file_path, sep='\t')

        # Initialisation de la liste de résultats
        results = []

        # Parcours des lignes du DataFrame
        for index, row in resultats_df.iterrows():
            if not row['#FILE'].startswith("#FILE") and not row['#FILE'].startswith("-"):
                current_entry = {}
                current_entry["file"] = cls.remove_extension(row['#FILE'].split("/")[2])
                current_entry["database"] = row['DATABASE']
                if row['DATABASE'] != "argannot":
                    if str(row['RESISTANCE']) != "nan":
                        current_entry["product_resistance"] = row['RESISTANCE'].lower().split(";")
                    else: current_entry["product_resistance"] = row['PRODUCT'].lower()
                else:
                    current_entry["product_resistance"] = row['PRODUCT'].lower()
                results.append(current_entry)

        return results

    @classmethod
    def info_from_results(cls, data):
        for dico in data:
            if dico["database"] == "argannot":
                antibio = re.search(r'\((.*?)\)', dico["product_resistance"]).group(1)
                antibio = cls.normalize_antibiotic_names_s(antibio)
                dico["product_resistance"] = antibio.split(";")
        return data