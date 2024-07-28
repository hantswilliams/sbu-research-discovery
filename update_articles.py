import os
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from Bio import Entrez
import smtplib
from email.mime.text import MIMEText

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///articles.db'
db = SQLAlchemy(app)

Entrez.email = os.getenv("EMAIL")

class FacultyMember(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), nullable=False, unique=True)

class DeletedArticle(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    pmid = db.Column(db.String(50), nullable=False, unique=True)

class Article(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    title = db.Column(db.String(500), nullable=False)
    authors = db.Column(db.String(500), nullable=False)
    source = db.Column(db.String(500), nullable=False)
    pub_date = db.Column(db.String(50), nullable=False)
    abstract = db.Column(db.Text, nullable=False)
    faculty_member = db.Column(db.String(100), nullable=False)
    pmid = db.Column(db.String(50), nullable=False)
    pubmed_link = db.Column(db.String(500), nullable=False)

def search_pubmed(author_name, retmax=50):
    search_term = f"{author_name}[Author]"
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    pmids = record["IdList"]
    return pmids

def fetch_article_details(pmids):
    try:
        handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        articles = []
        for record in records['PubmedArticle']:
            pub_date_dict = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
            pub_date = f"{pub_date_dict.get('Year', '')} {pub_date_dict.get('Month', '')} {pub_date_dict.get('Day', '')}".strip()
            
            abstract = record['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [""])[0]
            if isinstance(abstract, list):
                abstract = " ".join(abstract)
            elif isinstance(abstract, dict):
                abstract = abstract.get('text', '')

            authors = []
            for author in record['MedlineCitation']['Article']['AuthorList']:
                last_name = author.get('LastName', '')
                fore_name = author.get('ForeName', '')
                if last_name and fore_name:
                    authors.append(f"{last_name} {fore_name}")
                else:
                    authors.append(author.get('CollectiveName', ''))

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
    sender_email = os.getenv("EMAIL")
    password = os.getenv("EMAIL_PASSWORD")

    message = MIMEText(f"New articles found:\n\n" + "\n\n".join(f"{a['Title']} - {a['PubMedLink']}" for a in new_articles))
    message['Subject'] = 'New PubMed Articles Found'
    message['From'] = sender_email
    message['To'] = recipient_email

    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
        server.login(sender_email, password)
        server.sendmail(sender_email, recipient_email, message.as_string())

def update_articles():
    with app.app_context():
        faculty_members = FacultyMember.query.all()
        for faculty in faculty_members:
            pmids = search_pubmed(faculty.name, retmax=50)
            articles = fetch_article_details(pmids)
            new_articles = []
            for article_data in articles:
                if not Article.query.filter_by(pmid=article_data['PMID']).first() and not DeletedArticle.query.filter_by(pmid=article_data['PMID']).first():
                    new_article = Article(
                        title=article_data['Title'],
                        authors=article_data['Authors'],
                        source=article_data['Source'],
                        pub_date=article_data['PubDate'],
                        abstract=article_data['Abstract'],
                        faculty_member=faculty.name,
                        pmid=article_data['PMID'],
                        pubmed_link=article_data['PubMedLink']
                    )
                    db.session.add(new_article)
                    new_articles.append(article_data)
            db.session.commit()
            if new_articles:
                send_email(new_articles, "recipient@example.com")  # Change to the recipient's email address

if __name__ == "__main__":
    update_articles()
