import os
from dotenv import load_dotenv
import requests
import logging
import time

"""
API documentation: https://developer.clarivate.com/apis/wos-starter 
API Swagger: https://api.clarivate.com/swagger-ui/?url=https%3A%2F%2Fdeveloper.clarivate.com%2Fapis%2Fwos-starter%2Fswagger%3FforUser%3D73e459203fc3375aff081423b2349d8fdd4aa43d
"""

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Load environment variables
load_dotenv()

# Get the API key from environment variables
api_key = os.getenv('WEBSCIENCE_BASIC')
if not api_key:
    logging.error('API key not found. Please set the WEBSCIENCE_BASIC environment variable.')
    raise EnvironmentError('API key not found in environment variables')

# API endpoint
endpoint = 'https://api.clarivate.com/apis/wos-starter/v1'

# Set up headers with the API key
headers = {
    'X-ApiKey': api_key
}

def websci_byAuthor(author_lastname, author_firstname):
    # Validate input
    if not author_lastname or not author_firstname:
        logging.error('Both author_lastname and author_firstname must be provided')
        raise ValueError('Both author_lastname and author_firstname are required')

    logging.info(f'Getting articles for {author_firstname} {author_lastname}')
    
    try:
        # Initial request
        url = f'{endpoint}/documents?db=WOK&q=AU={author_lastname}, {author_firstname}&limit=50'
        response = requests.get(url, headers=headers)
        response.raise_for_status()  # Raise an HTTPError for bad responses
        response_json = response.json()

        # add in a 2 second delay between requests for rate limiting
        time.sleep(2)

        # Extract articles
        articles = response_json.get('hits', [])
        if not articles:
            logging.info(f'No articles found for {author_firstname} {author_lastname}')
            return None

        # Determine the total number of pages
        total_records = response_json.get('metadata', {}).get('total', 0)
        totalPages = (total_records // 50) + 1
        logging.info(f'Total pages: {totalPages}')

        # For testing purposes, keep only the first page
        if totalPages > 2:
            totalPages = 1

        # Loop through remaining pages if more exist
        for page in range(2, totalPages + 1):
            logging.info(f'Getting page {page}')
            url = f'{endpoint}/documents?db=WOK&q=AU={author_lastname}, {author_firstname}&limit=50&page={page}'
            response = requests.get(url, headers=headers)
            response.raise_for_status()  # Raise error if bad status
            response_json = response.json()
            articles += response_json.get('hits', [])

        # add in a key for the author name, the key is 'stonybrook_author' if len of articles > 0
        if len(articles) > 0:
            for article in articles:
                article['stonybrook_author'] = f'{author_firstname} {author_lastname}'

        return articles

    except requests.exceptions.HTTPError as http_err:
        logging.error(f'HTTP error occurred: {http_err}')
    except requests.exceptions.RequestException as req_err:
        logging.error(f'Request error occurred: {req_err}')
    except Exception as err:
        logging.error(f'An unexpected error occurred: {err}')

    return None
