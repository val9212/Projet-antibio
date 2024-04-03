from requester import *
import re

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
    def normalize_antibiotic_names(cls, antibiotic_list):
        antibiotic_normalization_dict = {
            "amp": "ampicillin",
            "amox": "amoxicillin",
            "agly": "aminoglycoside",
            "agly_flqn": "aminoglycoside",
            "bla": "beta-lactamase",
            " bla": "beta-lactamase",
            "fcd": "fidaxomicin",
            "fcyn": "cephalosporin",
            " ceph": "cephalosporin",
            "cef": "cefepime",
            "flq": "fluoroquinole",
            "cfx": "cefuroxime",
            " cfx": "cefuroxime",
            "gly": "glycopeptide",
            "mls": "macrolides",
            "phe": "phenylanin",
            "rif": "rifampicin",
            "sul": "sulfamide",
            "tet": "tetracyclin",
            "col": "colistine",
            "mupirocin": "mupirocin",
            "ntmdz": "ntmdz",
            "oxzln": "oxzln",
            "tmt": "tmt",
            "TetracenomycinC": "tetracyclin",
            "amoxil": "amoxicillin",
            "amoxicillin": "amoxicillin",
            "aminoglycoside": "aminoglycoside",
            "cephalosporin": "cephalosporin",
            "penam": "penam",
            "fosfomycin": "fosfomycin",
            "fluoroquinole": "fluoroquinole",
            "phenicol": "phenicol",
            "peptide": "peptide",
            "rifamycin": "rifamycin",
            "acridine_dye": "acridine_dye",
            "antibacterial_free_fatty_acids": "antibacterial_free_fatty_acids",
            "bicyclomycin": "bicyclomycin",
            "carbapenem": "carbapenem",
            "cephamycin": "cephamycin",
            "diaminopyrimidine": "diaminopyrimidine",
            "fusidic_acid": 'fusidic_acid',
            "glycopeptide": "glycopeptide",
            "lincosamide": "lincosamide",
            "macrolide": "macrolide",
            "phenicole": "phenicole",
            "pleuromutilin": "pleuromutilin",
            "streptogramin": "streptogramin",
            "sulfonamide": "sulfonamide",
            "tetracycline": "tetracycline",
            "triclosan": "triclosan",
            "16s": "16s",
            "23s": "23s",
            "d-lactate": "d-lactate",
            "abc-f": "abc-f",
            "alanine": "alanine",
            "arma": "arma",
            "beta-lactam": "beta-lactamase",
            "beta-lactamase": "beta-lactamase",
            "metallo-beta-lactamase": "beta-lactamase",
            "bleomycin": "bleomycin",
            "d-ala-d-ala": "d-ala-d-ala",
            "erythromycin": "erythromycin",
            "erm": "erythromycin",
            " erm": "erythromycin",
            "flavin": "flavin",
            "fomb": "fomb",
            "fos": "fos",
            "fusidic": "fusidic",
            "methicillin": "methicillin",
            "multidrug": "multidrug",
            "mupirocin-resistant": "mupirocin",
            "nitroimidazole": "nitroimidazole",
            " nim": "nitoimidazole",
            "nim": "nitoimidazole",
            "peptidoglycan": "peptidoglycan",
            "phosphoethanolamine--lipid": "phosphoethanolamine",
            "quinolone": "quinolone",
            "rifampin": "rifampin",
            "rifamycin-inactivating": "rifamycin",
            "serine": "serine",
            "trimethoprim": "trimethoprim",
            "vancomycin": "vancomycin",
            "van": "vancomycin",
            "aac": "Aminoglycoside",
            "aph": "Aminoglycoside",
            " aph": "Aminoglycoside",
            "viomycin": "viomycin",
            "zorbamycin": "zorbamycin",
            "amikacin": "amikacin",
            "carbomycin": "Carbomycin",
            "ceftazidime": "Ceftazidime",
            "chloramphenicol": "chloramphenicol",
            "ciprofloxacin": "ciprofloxacin",
            "colistin": "colistin",
            "gentamicin": "gentamicin",
            "imipenem": "imipenem",
            "lincomycin": "lincomycin",
            "oleandomycin": "oleandomycin",
            "spiramycin": "spiramycin",
            "streptomycin": "streptomycin",
            "sulfamethoxazole": "sulfamethoxazole",
            "tylosin": "tylosin",
            "fluoroquinolone": "fluoroquinolone",
            "streptothricin": "streptothricin",
            "aminocoumarin": "aminocoumarin",
            "elfamycin": "elfamycin",
            "nucleoside": "nucleoside",
            "ceftaroline-resistant": "ceftaroline",
            "tetracenomycin": "tetracenomycin",
            "kasugamycin": "kasugamycin",
            "clindamycin": "clindamycin",





        }
        normalized_list = []
        for name in antibiotic_list:
            # Convertir en minuscules
            name_lower = name.lower()
            # Remplacer par le nom normalisé s'il existe dans le dictionnaire
            if name_lower in antibiotic_normalization_dict:
                name_lower = antibiotic_normalization_dict[name_lower]
            # Ajouter le nom normalisé à la liste
            normalized_list.append(name_lower)
        return normalized_list

    @classmethod
    def trouver_mot_cle(cls, chaine):
        mots_cles = ["metallo-beta-lactamase", "amoxicillin", "Amikacin", "Clindamycin", "Cefepime", "Spiramycin", "Piperacillin+Tazobactam", "Ceftazidime", "Aztreonam", "Imipenem", "Tetracycline", "Gentamicin", "Dalfopristin", "doxycycline", "carbapenem", "Quinupristin", "Lincomycin", "rifamycin", "peptide", "cephalosporin", "D-lactate", "streptogramin", "mupirocin-resistant", "zorbamycin", "fomb", "fos", "phenicol", "teicoplanin", "VanY-like", "tetracenomycin", "ciprofloxacin", "flavin", "streptomycin", "nitroimidazole", "ceftaroline-resistant", "aminocyclitol", "penicillinase", "fusidic", "peptidoglycan", "rifamycin-inactivating", "viomycin", "lincosamide", "beta-lactamase", "bleomycin", "kasugamycin", "aminoglycoside", "arma", "rifampin", "chloramphenicol", "glycopeptide", "trimethoprim", "erythromycin", "23s", "fosfomycin", "abc-f", "phosphoethanolamine--lipid", "beta-lactam", "macrolide", "oleandomycin", "multidrug", "quinolone", "serine", "d-ala-d-ala", "vancomycin", "alanine", "tetracycline", "sulfonamide", "streptothricin", "16s", "fluoroquinolone", "multidrug", "beta-lactam-resistant" , "methicillin", "macrolide" ]
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
            trois_caracteres_apres_tilde = chaine[index_dernier_tilde + 4:index_dernier_tilde + 8]
            return trois_caracteres_apres_tilde
        else:
            # Si le dernier "~~~" n'est pas trouvé, retourner une chaîne vide
            return chaine