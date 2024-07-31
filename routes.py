from flask import render_template, request, jsonify, Response, redirect, url_for, session
from functools import wraps
from Bio import Entrez
from database import db, Article, FacultyMember, DeletedArticle  # Importing the database and models

Entrez.email = "your.email@example.com"  # Replace with your email

def register_routes(app):

    app.secret_key = 'your_secret_key_here'  # Replace with your secret key

    def login_required(func):
        @wraps(func)
        def decorated_function(*args, **kwargs):
            if 'logged_in' not in session:
                return redirect(url_for('login'))
            return func(*args, **kwargs)
        return decorated_function

    def check_password(func):
        @wraps(func)
        def decorated_function(*args, **kwargs):
            auth = request.authorization
            if not auth or auth.password != 'SHP2024':
                return Response('Could not verify your access level for that URL.\nYou have to login with proper credentials', 401,
                                {'WWW-Authenticate': 'Basic realm="Login Required"'})
            return func(*args, **kwargs)
        return decorated_function

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

    @app.route('/')
    def view():
        articles = Article.query.all()
        ## check if user is logged in
        if 'logged_in' in session:
            return render_template('view.html', articles=articles, logged_in=True)
        else:
            return render_template('view.html', articles=articles, logged_in=False)
        
    @app.route('/login', methods=['GET', 'POST'])
    def login():
        if request.method == 'POST':
            username = request.form['username']
            password = request.form['password']
            if username == 'SHP' and password == 'SHP2024!!':  # Replace with your credentials
                session['logged_in'] = True
                return redirect(url_for('view'))
            else:
                return render_template('login.html', error='Invalid credentials')
        return render_template('login.html')

    @app.route('/logout')
    def logout():
        session.pop('logged_in', None)
        return redirect(url_for('view'))

    @app.route('/search', methods=['GET', 'POST'])
    @login_required
    def index():
        if request.method == 'POST':
            author_name = request.form['author_name']
            article_limit = int(request.form['article_limit'])  # Get the article limit from the form
            pmids = search_pubmed(author_name, retmax=article_limit)
            print('PMIDs:', pmids)  # Debug print

            articles = fetch_article_details(pmids)

            # Fetch already saved articles
            saved_articles = Article.query.filter_by(faculty_member=author_name).all()
            saved_articles_set = set((a.title, a.authors, a.pub_date) for a in saved_articles)
            
            rendered_html = render_template('select.html', articles=articles, author_name=author_name, saved_articles_set=saved_articles_set)
            
            return rendered_html
        return render_template('index.html')

    @app.route('/get_article/<int:article_id>', methods=['GET'])
    def get_article(article_id):
        article = Article.query.get(article_id)
        if article:
            return jsonify({
                'id': article.id,
                'title': article.title,
                'authors': article.authors,
                'source': article.source,
                'pub_date': article.pub_date,
                'abstract': article.abstract,
                'faculty_member': article.faculty_member
            })
        return jsonify(error="Article not found"), 404

    @app.route('/edit_article', methods=['POST'])
    def edit_article():
        data = request.form
        article_id = data.get('id')
        article = Article.query.get(article_id)
        if article:
            article.title = data.get('title')
            article.authors = data.get('authors')
            article.source = data.get('source')
            article.pub_date = data.get('pub_date')
            article.abstract = data.get('abstract')
            article.faculty_member = data.get('faculty_member')
            db.session.commit()
            return jsonify(success=True)
        return jsonify(success=False, error="Article not found")


    @app.route('/add_article', methods=['POST'])
    @check_password
    def add_article():
        data = request.form
        try:
            new_article = Article(
                title=data['title'],
                authors=data['authors'],
                source=data['source'],
                pub_date=data['pub_date'],
                abstract=data['abstract'],
                faculty_member=data['faculty_member'],
                pmid=data['pmid'],
                pubmed_link=data['pubmed_link']
            )
            db.session.add(new_article)
            db.session.commit()
            return jsonify(success=True)
        except Exception as e:
            print(f"Error adding article: {e}")
            return jsonify(success=False, error=str(e))


    @app.route('/save', methods=['POST'])
    @login_required
    def save():
        print('Saving articles')  # Debug print
        data = request.get_json()
        faculty_member = data.get('faculty_member')
        selected_articles = data.get('selected_articles', [])
        print('faculty_member:', faculty_member)  # Debug print
        print('selected_articles:', selected_articles)  # Debug print

        # Ensure faculty member is saved
        faculty = FacultyMember.query.filter_by(name=faculty_member).first()
        if not faculty:
            faculty = FacultyMember(name=faculty_member)
            db.session.add(faculty)
            db.session.commit()

        for article_data in selected_articles:
            # Skip articles that are marked as deleted
            if DeletedArticle.query.filter_by(pmid=article_data['pmid']).first():
                continue

            try:
                new_article = Article(
                    title=article_data['title'],
                    authors=article_data['authors'][:450], ## limit authors to first 450 characters
                    source=article_data['source'],
                    pub_date=article_data['pub_date'],
                    abstract=article_data['abstract'],
                    faculty_member=faculty_member,
                    pmid=article_data['pmid'],
                    pubmed_link=article_data['pubmed_link']
                )
                db.session.add(new_article)
            except Exception as e:
                print(f"Error converting data: {e}")
                return jsonify(success=False)

        db.session.commit()
        return jsonify(success=True)

    @app.route('/saved')
    @login_required
    def saved():
        articles = Article.query.all()
        return render_template('saved.html', articles=articles)

    @app.route('/manage')
    @login_required
    def manage():
        articles = Article.query.all()
        return render_template('manage.html', articles=articles)

    @app.route('/delete_multiple', methods=['POST'])
    @login_required
    def delete_multiple():
        try:
            data = request.get_json()
            article_ids = data.get('article_ids', [])
            for article_id in article_ids:
                article = Article.query.get(article_id)
                if article:
                    # Add the PMIDs of deleted articles to the DeletedArticle table
                    deleted_article = DeletedArticle(pmid=article.pmid)
                    db.session.add(deleted_article)
                    db.session.delete(article)
            db.session.commit()
            return jsonify(success=True)
        except Exception as e:
            print(f"Error deleting articles: {e}")
            return jsonify(success=False, error=str(e))
