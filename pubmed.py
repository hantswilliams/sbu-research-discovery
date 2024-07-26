from Bio import Entrez

def search_pubmed(author_name, email, retmax=100):
    # Set email for Entrez
    Entrez.email = email

    # Construct the search term
    search_term = f"{author_name}[Author]"

    # Perform the search
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()

    # Get the list of PubMed IDs (PMIDs)
    pmids = record["IdList"]

    return pmids

def fetch_article_details(pmids):
    # Fetch the details of each article
    handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    # Extract and print details
    articles = []
    for record in records['PubmedArticle']:
        article = {
            'Title': record['MedlineCitation']['Article']['ArticleTitle'],
            'Authors': [author['LastName'] + " " + author['ForeName'] for author in record['MedlineCitation']['Article']['AuthorList']],
            'Source': record['MedlineCitation']['Article']['Journal']['Title'],
            'PubDate': record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'],
            'Abstract': record['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [""])[0]
        }
        articles.append(article)
    return articles

def main():
    author_name = input("Enter the author's name: ")
    email = input("Enter your email address: ")

    pmids = search_pubmed(author_name, email)

    if not pmids:
        print(f"No articles found for author: {author_name}")
        return

    articles = fetch_article_details(pmids)

    for idx, article in enumerate(articles, start=1):
        print(f"Article {idx}:")
        print(f"Title: {article['Title']}")
        print(f"Authors: {', '.join(article['Authors'])}")
        print(f"Source: {article['Source']}")
        print(f"Publication Date: {article['PubDate']}")
        print(f"Abstract: {article['Abstract']}\n")

if __name__ == "__main__":
    main()
