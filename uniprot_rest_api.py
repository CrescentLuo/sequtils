import requests, sys
import json
import xml.etree.ElementTree as ET

class UniProt:
    def __init__(self):
        self.requestURL = ''
        self.url_base = 'https://rest.uniprot.org'
        self.url_endpoint = 'uniprotkb'
        self.request = None

    def query(self, base='proteins',query='gene', query_value='AQP1',):
        url_base = self.url_base
        url_endpoint = self.url_endpoint
        url_query = 'search?&query={query}:{query_value}'.format(query=query, query_value=query_value)
        self.requestURL = '{}/{}/{}&size=1'.format(url_base, url_endpoint, url_query)

        self.request = requests.get(self.requestURL, headers={})

       
