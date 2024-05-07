import requests

from all import *
from requester import *
class GenomeData:
    def __init__(self, genome_id, genome_name, taxon_id):
        """
        Initialise un objet GenomeData avec un identifiant de génome, un nom de génome et un identifiant de taxon.

        :param genome_id: L'identifiant du génome.
        :param genome_name: Le nom du génome.
        :param taxon_id: L'identifiant du taxon.
        """
        self.__genome_id = genome_id
        self.__genome_name = genome_name
        self.__taxon_id = taxon_id
        self.__antibiotics_data = []
        self.__genus = None

    def add_antibiotic_data(self, antibiotic, measurement_value, measurement_unit, resistant_phenotype, laboratory_typing_method):
        """
        Ajoute des données sur les antibiotiques associées au génome.

        :param antibiotic: Le nom de l'antibiotique.
        :param measurement_value: La valeur de la mesure.
        :param measurement_unit: L'unité de mesure.
        :param resistant_phenotype: Le phénotype de résistance.
        :param laboratory_typing_method: La méthode de typage en laboratoire.
        """
        self.__antibiotics_data.append({
            'antibiotic': antibiotic,
            'measurement_value': measurement_value,
            'measurement_unit': measurement_unit,
            'resistant_phenotype': resistant_phenotype,
            'laboratory_typing_method': laboratory_typing_method
        })

    def get_antibiotics_data(self):
        """
        Renvoie les données sur les antibiotiques associées au génome.

        :return: Les données sur les antibiotiques.
        """
        return self.__antibiotics_data

    def get_id(self):
        """
        Renvoie l'identifiant du génome.

        :return: L'identifiant du génome.
        """
        return self.__genome_id

    def get_name(self):
        """
        Renvoie le nom du génome.

        :return: Le nom du génome.
        """
        return self.__genome_name

    def get_taxon_id(self):
        """
        Renvoie l'identifiant du taxon associé au génome.

        :return: L'identifiant du taxon.
        """
        if self.__taxon_id:
            return int(self.__taxon_id)
        else:
            return self.__genome_id

    def get_genus(self):
        """
        Renvoie le genre associé au génome.

        :return: Le genre du génome.
        """
        return self.__genus

    def add_genus(self, genus):
        """
        Ajoute le genre associé au génome.

        :param genus: Le genre du génome.
        """
        self.__genus = genus

    def __str__(self):
        """
        Renvoie une représentation sous forme de chaîne de l'identifiant du taxon.

        :return: L'identifiant du taxon sous forme de chaîne.
        """
        return (f"{self.__taxon_id}")

    def __repr__(self):
        """
        Renvoie une représentation sous forme de chaîne de l'objet GenomeData.

        :return: Une représentation de l'objet GenomeData.
        """
        return str(self)