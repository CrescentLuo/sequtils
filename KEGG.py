from ast import arguments
import requests, json

kegg_operation_set = {
    'entry',
    'pathway',
    'brite',
    'module',
    'network'}

org_set = {
    'hsa':'Homo sapiens', # human
    'mmu':'Mus musculus', # mouse
    'rno':'Rattus norvegicus', # rat
    'ocu':'Oryctolagus cuniculus' # rabbit
    }

'''
wget -c https://rest.kegg.jp/get/hsa:2538+mcf:102144494+ssc:100134959+ocu:100303766+rno:25634+mmu:14377
wget -c http://rest.kegg.jp/get/hsa:2538+mcf:102144494+ssc:100134959+ocu:100303766+rno:25634+mmu:14377/aaseq
wget -c http://rest.kegg.jp/get/hsa:2538+mcf:102144494+ssc:100134959+ocu:100303766+rno:25634+mmu:14377/ntseq
'''

#<org> = KEGG organism code


operation = ''
argument = ''
kegg_url = 'https://www.kegg.jp/{}/{}'.format(operation,argument)
user_passwd = ""
data = json.dumps('')

r = requests.post(kegg_url,)
