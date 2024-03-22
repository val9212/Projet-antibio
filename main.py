from InOut import *

genome_objects, antibiotics = Inout.read_file("./data/genome_AMR.txt")

listes = []
for i in genome_objects:
    integer_part = i.split('.')[0]
    listes.append(int(integer_part))
y= set(listes)


def write_list_to_file(data_list, file_name='test.txt'):
    """
    Writes the elements of data_list to file_name, separated by commas.

    :param data_list: List of items to write to the file.
    :param file_name: The name of the file to create and write to.
    """
    try:
        # Open the file in write mode ('w') to create it or overwrite if it already exists
        with open(file_name, 'w') as file:
            # Join the list items into a single string separated by commas
            file.write(','.join(str(item) for item in data_list))
        print(f"Data successfully written to {file_name}")
    except IOError as e:
        # Handle the error and print an error message
        print(f"An error occurred: {e.strerror}")

write_list_to_file(y)
