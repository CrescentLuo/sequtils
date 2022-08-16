import requests, sys
import json
import xml.etree.ElementTree as ET
import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry

def get_next_link(headers):
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

def get_batch(session, batch_response, file_format, compressed):
        batch_url = get_next_link(batch_response.headers)
        while batch_url:
            batch_response = session.get(batch_url)
            batch_response.raise_for_status()
            yield decode_results(batch_response, file_format, compressed)
            batch_url = get_next_link(batch_response.headers)

def combine_batches(all_results, batch_results, file_format):
        if file_format == "json":
            for key in ("results", "failedIds"):
                if key in batch_results and batch_results[key]:
                    all_results[key] += batch_results[key]
        elif file_format == "tsv":
            return all_results + batch_results[1:]
        else:
            return all_results + batch_results
        return all_results

def print_progress_batches(batch_index, size, total):
        n_fetched = min((batch_index + 1) * size, total)
        print(f"Fetched: {n_fetched} / {total}")

def decode_results(response, file_format, compressed):
        if compressed:
            decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
            if file_format == "json":
                j = json.loads(decompressed.decode("utf-8"))
                return j
            elif file_format == "tsv":
                return [line for line in decompressed.decode("utf-8").split("\n") if line]
            elif file_format == "xlsx":
                return [decompressed]
            elif file_format == "xml":
                return [decompressed.decode("utf-8")]
            else:
                return decompressed.decode("utf-8")
        elif file_format == "json":
            return response.json()
        elif file_format == "tsv":
            return [line for line in response.text.split("\n") if line]
        elif file_format == "xlsx":
            return [response.content]
        elif file_format == "xml":
            return [response.text]
        return response.text

def get_xml_namespace(element):
        m = re.match(r"\{(.*)\}", element.tag)
        return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)

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

class UniprotIDMapper:
    def __init__(self):
        self.POLLING_INTERVAL = 3
        self.API_URL = "https://rest.uniprot.org"
        self.retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=self.retries))
        print("session initiated")


    def submit_id_mapping(self, from_db, to_db, ids):
        request = requests.post(
            f"{self.API_URL}/idmapping/run",
            data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
        )
        request.raise_for_status()
        return request.json()["jobId"]

    def check_id_mapping_results_ready(self, job_id):
        while True:
            request = self.session.get(f"{self.API_URL}/idmapping/status/{job_id}")
            request.raise_for_status()
            j = request.json()
            if "jobStatus" in j:
                if j["jobStatus"] == "RUNNING":
                    print(f"Retrying in {self.POLLING_INTERVAL}s")
                    time.sleep(self.POLLING_INTERVAL)
                else:
                    raise Exception(request["jobStatus"])
            else:
                return bool(j["results"] or j["failedIds"])
    
    def get_id_mapping_results_link(self,job_id):
        url = f"{self.API_URL}/idmapping/details/{job_id}"
        request = self.session.get(url)
        request.raise_for_status()
        return request.json()["redirectURL"]
    
    def get_id_mapping_results_stream(self, url):
        if "/stream/" not in url:
            url = url.replace("/results/", "/stream/")
        request = self.session.get(url)
        request.raise_for_status()
        parsed = urlparse(url)
        query = parse_qs(parsed.query)
        file_format = query["format"][0] if "format" in query else "json"
        compressed = (
            query["compressed"][0].lower() == "true" if "compressed" in query else False
        )
        return decode_results(request, file_format, compressed)

    def get_id_mapping_results_search(self, url):
        parsed = urlparse(url)
        query = parse_qs(parsed.query)
        file_format = query["format"][0] if "format" in query else "json"
        if "size" in query:
            size = int(query["size"][0])
        else:
            size = 500
            query["size"] = size
        compressed = (
            query["compressed"][0].lower() == "true" if "compressed" in query else False
        )
        parsed = parsed._replace(query=urlencode(query, doseq=True))
        url = parsed.geturl()
        request = self.session.get(url)
        request.raise_for_status()
        results = decode_results(request, file_format, compressed)
        total = int(request.headers["x-total-results"])
        print_progress_batches(0, size, total)
        for i, batch in enumerate(get_batch(self.session, request, file_format, compressed), 1):
            results = combine_batches(results, batch, file_format)
            print_progress_batches(i, size, total)
        if file_format == "xml":
            return merge_xml_results(results)
        return results

    def map_id(self, from_db, to_db, ids):
        job_id = self.submit_id_mapping(
            from_db, to_db, ids)
        print(job_id)
        if self.check_id_mapping_results_ready(job_id):
            link = self.get_id_mapping_results_link(job_id)
            results = self.get_id_mapping_results_search(link)
            #results = self.get_id_mapping_results_stream(link)
            return results
        # Equivalently using the stream endpoint which is more demanding
        # on the API and so is less stable:
        # results = get_id_mapping_results_stream(link)

        
        # {'results': [{'from': 'P05067', 'to': 'CHEMBL2487'}], 'failedIds': ['P12345']}
