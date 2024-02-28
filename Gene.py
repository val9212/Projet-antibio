class GenomeData:
    def __init__(self, genome_id, genome_name):
        self.__genome_id = genome_id
        self.__genome_name = genome_name
        self.__antibiotics_data = []

    def add_antibiotic_data(self, antibiotic, measurement_value, measurement_unit, resistant_phenotype, laboratory_typing_method):
        self.__antibiotics_data.append({
            'antibiotic': antibiotic,
            'measurement_value': measurement_value,
            'measurement_unit': measurement_unit,
            'resistant_phenotype': resistant_phenotype,
            'laboratory_typing_method': laboratory_typing_method
        })

    def get_antibiotics_data(self):
        return self.__antibiotics_data

    def get_id(self):
        return self.__genome_id

    def get_name(self):
        return self.__genome_name

    def __str__(self):
        return (f"{self.__genome_id}")

    def __repr__(self):
        return str(self)