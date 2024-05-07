from requester import *
import re
import csv
import os
from Bio.Data import IUPACData
import pandas as pd

class methods():

    @classmethod
    def unique_genus(cls, genus_list):
        """
        Retourne une liste unique des noms de genre à partir d'une liste de listes de noms de genre.

        :param genus_list: Une liste de listes de noms de genre.
        :return: Une liste unique des noms de genre.
        """
        final_list = []
        for sublist in genus_list:
            for item in sublist:
                if item not in final_list:
                    final_list.append(item)
        return final_list

    @classmethod
    def devide_requester(cls, taxon_id, n_chunk = 10):
        """
        Divise les identifiants de taxon en chunks et récupère les noms de genre correspondants pour chaque chunk.

        :param taxon_id: Une liste d'identifiants de taxon.
        :param n_chunk: Le nombre de chunks à diviser les identifiants de taxon.
        :return: Une liste de listes de noms de genre correspondants à chaque chunk.
        """
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
        """
        Divise les identifiants de taxon en chunks et récupère les noms de genre correspondants pour chaque chunk.

        :param taxon_id: Une liste d'identifiants de taxon.
        :param n_chunk: Le nombre de chunks à diviser les identifiants de taxon.
        :return: Un dictionnaire où chaque identifiant de taxon est associé à son nom de genre correspondant.
        """
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
        """
        Extrait les en-têtes des séquences à partir d'un fichier FASTA.

        :param file: Le chemin vers le fichier FASTA.
        :return: Une liste des en-têtes des séquences.
        """
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
    def get_anti(cls, headers, param = "argannot"):
        """
        Extrait les noms des antibiotiques à partir des en-têtes des séquences.

        :param headers: Une liste des en-têtes des séquences.
        :param param: Le paramètre indiquant la manière d'extraire les noms des antibiotiques (par défaut: 'argannot').
        :return: Une liste des noms des antibiotiques extraits des en-têtes des séquences.
        """
        antibiotic_classes = []
        if param != "argannot":
            for header in headers:
                parties = header.rsplit("~~~", 3)
                derniers_trois_tildes = parties[3]
                nparties = derniers_trois_tildes.split(" ")
                final = nparties[0]
                liste_header = final.lower().split(";")
                for i in liste_header:
                    if i == '':
                        final = nparties[1]
                        tfinale = cls.normalize_antibiotic_names_s(final[:3])
                        antibiotic_classes.append(tfinale)
                    else: antibiotic_classes.append(i)
        else:
            for header in headers:
                parties = header.rsplit("~~~", 3)
                all = parties[1]
                antibio = re.search(r'\((.*?)\)', all).group(1)
                antibio = cls.normalize_antibiotic_names_s(antibio)
                antibiotic_classes.append(antibio)
        return antibiotic_classes

    @classmethod
    def normalize_antibiotic_names_s(cls, name):
        """
        Normalise un nom d'antibiotique en utilisant un dictionnaire de correspondance prédéfini.

        :param name: Le nom de l'antibiotique à normaliser.
        :return: Le nom normalisé de l'antibiotique.
        """
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
        """
        Normalise une liste de noms d'antibiotiques en utilisant un dictionnaire de correspondance prédéfini.

        :param antibiotic_list: La liste des noms d'antibiotiques à normaliser.
        :return: La liste des noms d'antibiotiques normalisés.
        """

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
        """
        Normalise le nom de l'antibiotique en utilisant la table de données IUPAC.

        :param name: Le nom de l'antibiotique à normaliser.
        :return: Le nom normalisé de l'antibiotique.
        """
        normalized_name = IUPACData.ambiguous_dna_values.get(name, name)
        return normalized_name

    @classmethod
    def trouver_mot_cle(cls, chaine):
        """
        Trouve un mot clé dans une chaîne donnée.

        :param chaine: La chaîne dans laquelle rechercher le mot clé.
        :return: Le premier mot clé trouvé dans la chaîne.
        """
        mots_cles = ["amoxicillin", "Amikacin", "Clindamycin", "Cefepime", "Spiramycin", "Piperacillin+Tazobactam", "Ceftazidime", "Aztreonam", "Imipenem", "Tetracycline", "Gentamicin", "Dalfopristin", "doxycycline", "carbapenem", "Quinupristin", "Lincomycin", "rifamycin", "peptide", "cephalosporin", "D-lactate", "streptogramin", "mupirocin-resistant", "zorbamycin", "fomb", "fos", "phenicol", "polymyxin", "teicoplanin", "VanY-like", "tetracenomycin", "ciprofloxacin", "flavin", "streptomycin", "nitroimidazole", "ceftaroline-resistant", "aminocyclitol", "penicillinase", "fusidic", "peptidoglycan", "rifamycin-inactivating", "viomycin", "lincosamide", "beta-lactamase", "bleomycin", "kasugamycin", "aminoglycoside", "arma", "rifampin", "chloramphenicol", "glycopeptide", "trimethoprim", "erythromycin", "23s", "fosfomycin", "abc-f", "phosphoethanolamine--lipid", "beta-lactam", "macrolide", "oleandomycin", "multidrug", "quinolone", "serine", "d-ala-d-ala", "vancomycin", "alanine", "tetracycline", "sulfonamide", "streptothricin", "16s", "fluoroquinolone", "multidrug", "beta-lactam-resistant" , "methicillin", "macrolide", "metallo-beta-lactamase" ]
        for mot_cle in mots_cles:
            if mot_cle.lower() in chaine.lower():
                return mot_cle

        return cls.extraire_trois_caracteres_apres_tilde(chaine.lower())

    @classmethod
    def count_element(cls, liste):
        """
        Compte le nombre d'occurrences de chaque élément dans une liste.

        :param liste: La liste contenant les éléments à compter.
        :return: Un dictionnaire contenant les éléments comme clés et leur nombre d'occurrences comme valeurs.
        """
        occurrences_dict = {}
        for element in liste:
            if element in occurrences_dict:
                occurrences_dict[element] += 1
            else:
                occurrences_dict[element] = 1
        return occurrences_dict

    @classmethod
    def present_elements(cls, liste):
        """
        Vérifie la présence de chaque élément dans une liste.

        :param liste: La liste contenant les éléments à vérifier.
        :return: Un dictionnaire indiquant la présence de chaque élément dans la liste.
        """
        dictio = {}
        for element in liste:
            if element not in dictio:
                dictio[element] = 1

        return dictio

    @classmethod
    def extraire_trois_caracteres_apres_tilde(cls, chaine):
        """
        Extrait les trois caractères suivant le dernier "~~~" dans une chaîne.

        :param chaine: La chaîne dans laquelle effectuer l'extraction.
        :return: Les trois caractères suivant le dernier "~~~" dans la chaîne.
        """
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
        """
        Associe chaque nom d'antibiotique à une classe d'antibiotique spécifique.

        :param normalized_list: Liste de noms normalisés d'antibiotiques.
        :return: Liste des classes d'antibiotiques associées aux noms d'antibiotiques.
        """
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
        """
        Compte le nombre d'occurrences de résistance/susceptibilité/intermédiaire/erreur de données
        pour chaque antibiotique dans un fichier CSV de résultats.

        :param fichier_csv: Le chemin vers le fichier CSV contenant les données de résistance.
        """
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
        """
        Obtient les informations des génomes avec le plus grand nombre de tests effectués.

        :param genome_object: Un dictionnaire contenant des objets de génome.
        :param no_missing: Un booléen indiquant s'il faut exclure les tests manquants.
        """
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
        """
        Obtient les informations des génomes sélectionnés à partir d'un fichier.

        :param file: Le chemin vers le fichier contenant les informations des génomes sélectionnés.
        :return: Une liste contenant les noms des genres et les identifiants des génomes sélectionnés.
        """
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
        """
        Supprime l'extension du nom de fichier.

        :param file_name: Le nom du fichier avec extension.
        :return: Le nom du fichier sans extension.
        """
        return os.path.splitext(file_name)[0]

    @classmethod
    def parse_argannot_results(cls, file_path):
        """
        Analyse les résultats provenant de la base de données Arg-Annot.

        :param file_path: Chemin vers le fichier CSV contenant les résultats de la recherche.
        :return: Une liste de dictionnaires contenant les informations analysées.
        """
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
        """
        Obtient des informations à partir des résultats analysés.

        :param data: Liste de dictionnaires contenant les résultats analysés.
        :return: Une liste de dictionnaires contenant les informations obtenues.
        """
        for dico in data:
            if dico["database"] == "argannot":
                antibio = re.search(r'\((.*?)\)', dico["product_resistance"]).group(1)
                antibio = cls.normalize_antibiotic_names_s(antibio)
                dico["product_resistance"] = antibio.split(";")
        return data