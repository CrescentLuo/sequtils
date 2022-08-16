from ast import arguments
import requests, json
import re
from requests.adapters import HTTPAdapter, Retry

kegg_operation_set = {
    'entry',
    'pathway',
    'brite',
    'module',
    'network'}



'''
wget -c https://rest.kegg.jp/get/hsa:2538+mcf:102144494+ssc:100134959+ocu:100303766+rno:25634+mmu:14377
wget -c http://rest.kegg.jp/get/hsa:2538+mcf:102144494+ssc:100134959+ocu:100303766+rno:25634+mmu:14377/aaseq
wget -c http://rest.kegg.jp/get/hsa:2538+mcf:102144494+ssc:100134959+ocu:100303766+rno:25634+mmu:14377/ntseq
'''

def get_aa_multi_organism(kegg):
    start = False
    gene_entry = parse_entry(kegg.request)
    orth_id = gene_entry['ORTHOLOGY'][0]
    org_set = {
            'hsa':'Homo sapiens', # human
            'mcf':'Macaca fascicularis', # crab-eating macaque
            'mmu':'Mus musculus', # mouse
            'rno':'Rattus norvegicus', # rat
            'ocu':'Oryctolagus cuniculus', # rabbit
            'cfa':'Canis lupus familiaris', # dog
            'ssc':'Sus scrofa', # pig
            }
    orth_requst = kegg.operation_get(orth_id)
    orth_entry = parse_entry(orth_requst)
    seq_set = {}
    for org in org_set:
        query = org + ':' +orth_entry['GENES'][org][0]
        seq_id, seq = kegg.operation_get_aaseq(query)
        seq_set[seq_id] = seq

    return seq_set


def parse_entry(request):
    entry_dict = {}
    iterm = ''
    for line in request.split('\n'):
        sline = line.split(' ')
        if sline[0] == '':
            vline = [x for x in sline[1:] if x]
            if item == 'GENES':
                id, *name = vline[1].rstrip(')').split('(')
                if name:
                    name = name[0]
                else:
                    name = 'Null'
                entry_dict[item][vline[0].lower()[:-1]] = [id,name]
            else:
                entry_dict[item].append(vline)
        else:
            item = sline[0]
            vline = [x for x in sline[1:] if x]
            if item == 'GENES':
                id, *name= vline[1].rstrip(')').split('(')
                if name:
                    name = name[0]
                else:
                    name = 'Null'
                entry_dict[item] = {vline[0].lower()[:-1]:[id,name]}
            else:
                entry_dict[item]=vline
            
    return entry_dict

class KEGG:
    def __init__(self, query=None):
        self.POLLING_INTERVAL = 3
        self.API_URL = "https://rest.kegg.jp"
        self.retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=self.retries))
        print("session initiated")
        if query:
            print("fetching entries with id:{}".format(query))
            self.request = self.operation_get(query)
            self.orth_id = None
            self.entry = parse_entry(self.request)

    def operation_get(self, entrez_id):
        request = requests.post(
            f"{self.API_URL}/get/{entrez_id}",
        )
        while True:
            request = self.session.get(f"{self.API_URL}/get/{entrez_id}")
            request.raise_for_status()
            j = request.text
            if request.ok:
                return j
    def operation_get_aaseq(self, entrez_id):
        request = requests.post(
            f"{self.API_URL}/get/{entrez_id}/aaseq",
        )
        while True:
            request = self.session.get(f"{self.API_URL}/get/{entrez_id}/aaseq")
            request.raise_for_status()
            j = request.text
            sline = j.split('\n')
            tag = sline[0]
            seq = ''.join(sline[1:])
            if request.ok:
                return tag, seq



