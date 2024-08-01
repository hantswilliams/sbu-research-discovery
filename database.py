from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class Article(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    title = db.Column(db.String(500), nullable=False)
    authors = db.Column(db.String(10000), nullable=False)
    source = db.Column(db.String(500), nullable=False)
    pub_date = db.Column(db.String(50), nullable=False)
    abstract = db.Column(db.Text, nullable=False)
    faculty_member = db.Column(db.String(100), nullable=False)
    pmid = db.Column(db.String(50), nullable=False, unique=True)  # Add unique constraint for PubMed ID
    pubmed_link = db.Column(db.String(500), nullable=False)
    keywords = db.relationship('ArticleKeyword', backref='article', cascade='all, delete-orphan', lazy=True)

class FacultyMember(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), nullable=False, unique=True)
class ArticleKeyword(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    pmid = db.Column(db.String(50), db.ForeignKey('article.pmid'), nullable=False)
    keyword = db.Column(db.String(500), nullable=False)

class DeletedArticle(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    pmid = db.Column(db.String(50), nullable=False, unique=True)
