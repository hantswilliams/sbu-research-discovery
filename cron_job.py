import os
import requests
from dotenv import load_dotenv
from Bio import Entrez
from app import create_app
from database import db, Article, FacultyMember, DeletedArticle

load_dotenv()

# Set the email for Entrez
Entrez.email = "hants.williams@stonybrook.edu"

def search_pubmed(author_name, retmax=200):
    search_term = f"{author_name}[Author]"
    try:
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Error searching PubMed: {e}")
        return []
    
    pmids = record.get("IdList", [])
    return pmids

def fetch_article_details(pmids):
    if not pmids:
        print("Received empty list of pmids")
        return []

    try:
        handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        articles = []
        for record in records['PubmedArticle']:
            try:
                pub_date_dict = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
                pub_date = f"{pub_date_dict.get('Year', '')} {pub_date_dict.get('Month', '')} {pub_date_dict.get('Day', '')}".strip()
            except Exception as e:
                print(f"Error fetching publication date: {e}")
                pub_date = ""

            try:
                abstract = record['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [""])[0]
                if isinstance(abstract, list):
                    abstract = " ".join(abstract)
                elif isinstance(abstract, dict):
                    abstract = abstract.get('text', '')
            except Exception as e:
                print(f"Error fetching abstract: {e}")
                abstract = ""

            authors = []
            for author in record['MedlineCitation']['Article']['AuthorList']:
                last_name = author.get('LastName', '')
                fore_name = author.get('ForeName', '')
                if last_name and fore_name:
                    authors.append(f"{last_name} {fore_name}")
                else:
                    authors.append(author.get('CollectiveName', ''))

            authors = ", ".join(authors)
            authors = authors[:450]  # Limit to first 450 characters
            pmid = record['MedlineCitation']['PMID']
            pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            article = {
                'Title': record['MedlineCitation']['Article']['ArticleTitle'],
                'Authors': authors,
                'Source': record['MedlineCitation']['Article']['Journal']['Title'],
                'PubDate': pub_date,
                'Abstract': str(abstract),
                'PMID': str(pmid),
                'PubMedLink': str(pubmed_link)
            }
            articles.append(article)
        return articles
    except Exception as e:
        print(f"Error fetching article details: {e}")
        return []

def send_email(new_articles, recipient_email):
    print("New articles found, sending email...")
    
    mailgun_api_key = os.getenv("MAILGUN_API_KEY")
    mailgun_domain = os.getenv("MAILGUN_DOMAIN")
    sender_email = os.getenv("MAILGUN_FROM_EMAIL")

    if not mailgun_api_key or not mailgun_domain or not sender_email:
        print("Mailgun configuration is missing. Please set MAILGUN_API_KEY, MAILGUN_DOMAIN, and MAILGUN_FROM_EMAIL environment variables.")
        return

    html_content = f"New articles found:<br><br>" + "<br><br>".join(
        f"<strong>{a['FacultyMember']}:</strong> <a href='{a['PubMedLink']}'>{a['Title']}</a>" for a in new_articles)

    response = requests.post(
        f"https://api.mailgun.net/v3/{mailgun_domain}/messages",
        auth=("api", mailgun_api_key),
        data={
            "from": sender_email,
            "to": recipient_email,
            "subject": "New PubMed Articles Found",
            "html": html_content
        }
    )

    if response.status_code == 200:
        print("Email sent successfully!")
    else:
        print(f"Failed to send email. Status code: {response.status_code}, response: {response.text}")

def update_articles():
    new_articles = []
    app = create_app()
    with app.app_context():
        faculty_members = FacultyMember.query.all()
        for faculty in faculty_members:
            print(f"Updating articles for {faculty.name}")
            pmids = search_pubmed(faculty.name, retmax=50)
            print(f"Found {len(pmids)} all time articles for {faculty.name}")
            articles = fetch_article_details(pmids)
            for article_data in articles:
                if not Article.query.filter_by(pmid=article_data['PMID']).first() and not DeletedArticle.query.filter_by(pmid=article_data['PMID']).first():
                    new_article = Article(
                        title=article_data['Title'],
                        authors=article_data['Authors'][:450],                         ## limit to first 450 characters
                        source=article_data['Source'],
                        pub_date=article_data['PubDate'],
                        abstract=article_data['Abstract'],
                        faculty_member=faculty.name,
                        pmid=article_data['PMID'],
                        pubmed_link=article_data['PubMedLink']
                    )
                    db.session.add(new_article)
                    article_data['FacultyMember'] = faculty.name  # Add faculty member name to article data
                    new_articles.append(article_data)
            db.session.commit()
            print(f"New articles found for {faculty.name} to be added: {len(new_articles)}")
    
    if new_articles:
        send_email(new_articles, "hants.williams@stonybrook.edu")

if __name__ == "__main__":
    update_articles()