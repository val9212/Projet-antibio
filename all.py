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

class db():
    @classmethod
    get




