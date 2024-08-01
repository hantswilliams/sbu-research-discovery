from Bio import Entrez
import pandas as pd

# Provide your email address
Entrez.email = "hants.williams@stonybrook.edu"

### Search for a specific author
search_term = f"Hants Williams[Author]"
handle = Entrez.esearch(db="pubmed", term=search_term, retmax=100)
record = Entrez.read(handle)
handle.close()
pmids = record["IdList"]

### Fetch the records for the PMIDs
handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
records = Entrez.read(handle)
handle.close()

###
articles_keywords = {}
for record in records['PubmedArticle']:
    ## create a new dictionary with PMID as the key
    article = {}
    article['pmid'] = record['MedlineCitation']['PMID']
    article['data'] = record
    articles_keywords[article['pmid']] = article

len(articles_keywords)

articles_keywords.keys()
articles_keywords['38363587']['data']['MedlineCitation'].keys()

## other potential keys of interest:
### KeywordList
"""
'KeywordList': [ListElement([StringElement('PROM', attributes={'MajorTopicYN': 'N'}), StringElement('PSI', attributes={'MajorTopicYN': 'N'}), StringElement('Parsley Symptom Index', attributes={'MajorTopicYN': 'N'}), StringElement('chronic disease', attributes={'MajorTopicYN': 'N'}), StringElement('eHealth', attributes={'MajorTopicYN': 'N'}), StringElement('ePROM', attributes={'MajorTopicYN': 'N'}), StringElement('mHealth', attributes={'MajorTopicYN': 'N'}), StringElement('patient-reported outcome measure', attributes={'MajorTopicYN': 'N'}), StringElement('telehealth', attributes={'MajorTopicYN': 'N'}), StringElement('telemedicine', attributes={'MajorTopicYN': 'N'}), StringElement('validation', attributes={'MajorTopicYN': 'N'}), StringElement('web-based', attributes={'MajorTopicYN': 'N'})],
"""

### Retrieve KeywordList from each key in articles_keywords and print
article_keys_df = pd.DataFrame(columns=['pmid', 'keyword'])
for key in articles_keywords.keys():
    keys = articles_keywords[key]['data']['MedlineCitation']['KeywordList']
    for k in keys:
        ## find StringElement in k
        for i in k:
            print(f'{key}: {i}')
            article_keys_df = article_keys_df._append({'pmid': key, 'keyword': i}, ignore_index=True)



            