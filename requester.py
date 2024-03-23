import requests
import xml.etree.ElementTree as ET

class request():

    @classmethod
    def requester(cls, taxon_id):
        # Define the base URL
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

        # Parameters
        db = 'taxonomy'
        id_list = taxon_id

        # EPost - Posts the UID list to the Entrez History server
        epost_url = f"{base_url}epost.fcgi?db={db}&id={id_list}"
        response = requests.post(epost_url)
        if not response.ok:
            raise Exception("Error posting epost request")

        # Parse WebEnv and QueryKey
        web_env = response.text.split('<WebEnv>')[1].split('</WebEnv>')[0]
        query_key = response.text.split('<QueryKey>')[1].split('</QueryKey>')[0]

        # ESummary - Retrieves summary information
        esummary_url = f"{base_url}esummary.fcgi?db={db}&query_key={query_key}&WebEnv={web_env}"
        summary_response = requests.get(esummary_url)
        if summary_response.ok:
            # Output the summaries
            return summary_response.text
        else:
            raise Exception("Error fetching esummary")

        # EFetch - Retrieves records in the requested format
        efetch_url = f"{base_url}efetch.fcgi?db={db}&query_key={query_key}&WebEnv={web_env}"
        efetch_url += "&rettype=fasta&retmode=text"
        fetch_response = requests.get(efetch_url)
        if fetch_response.ok:
            # Output the fetched data
            return fetch_response.text
        else:
            raise Exception("Error fetching efetch")

    @classmethod
    def extract_genus(cls, xml_document):
        # Parse the entire XML document
        root = ET.fromstring(xml_document)

        # Initialize an empty list to store all genus names
        genus_list = []

        # Iterate over all 'DocSum' elements in the document
        for doc_sum in root.findall('DocSum'):
            # Find the 'Item' element with the name 'Genus'
            genus_item = doc_sum.find(".//Item[@Name='Genus']")
            if genus_item is not None:
                # Append the genus name to the list
                genus_list.append(genus_item.text)

        return genus_list
