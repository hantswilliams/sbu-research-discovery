import os 
from dotenv import load_dotenv
import requests

load_dotenv()

print(os.getenv('WEBSCIENCE_BASIC'))

endpoint = 'https://api.clarivate.com/apis/wos-starter/v1'

headers = {
    'X-ApiKey': os.getenv('WEBSCIENCE_BASIC')
}

## documents endpoint; searching for AU = "Williams, Hants"
# and where db = 'WOK'
# and limit to 50 records

url = f'{endpoint}/documents?db=WOK&q=AU=Williams, Hants&limit=50'

response = requests.get(url, headers=headers)
response.text