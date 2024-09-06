import os 
from dotenv import load_dotenv
import requests
import pandas as pd 

"""

API documentation: https://developer.clarivate.com/apis/wos-starter 
API Swagger: https://api.clarivate.com/swagger-ui/?url=https%3A%2F%2Fdeveloper.clarivate.com%2Fapis%2Fwos-starter%2Fswagger%3FforUser%3D73e459203fc3375aff081423b2349d8fdd4aa43d

"""

load_dotenv()

print(os.getenv('WEBSCIENCE_BASIC'))

endpoint = 'https://api.clarivate.com/apis/wos-starter/v1'

headers = {
    'X-ApiKey': os.getenv('WEBSCIENCE_BASIC')
}

def websci_byAuthor(author_lastname, author_firstname):
    print(f'Getting articles for {author_firstname} {author_lastname}')
    url = f'{endpoint}/documents?db=WOK&q=AU={author_lastname}, {author_firstname}&limit=50'
    response = requests.get(url, headers=headers)
    response_json = response.json()
    articles = response_json['hits']
    if len(articles) == 0:
        return None
    
    ## assume 50 records per page (total)
    totalPages = (response_json['metadata']['total'] // 50) + 1
    print(f'Total pages: {totalPages}')

    for page in range(2, totalPages+1):
        print(f'Getting page {page}')
        url = f'{endpoint}/documents?db=WOK&q=AU={author_lastname}, {author_firstname}&limit=50&page={page}'
        response = requests.get(url, headers=headers)
        response_json = response.json()
        articles += response_json['hits']

    return articles

## test with Lisa Muratori
articles = websci_byAuthor('Muratori', 'Lisa')
len(articles)





## for testing purposes, flatten and convert to pandas dataframe
df = pd.json_normalize(articles)
## add a column for author, make sure its the first column
df['sbu_shp_employee'] = 'Lisa Muratori'
cols = df.columns.tolist()
cols = cols[-1:] + cols[:-1]
df = df[cols]
df.to_csv('./webofscience/test_articles.csv', index=False)

