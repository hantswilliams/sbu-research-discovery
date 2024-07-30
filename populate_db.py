import pandas as pd
from app import create_app
from database import db, FacultyMember

# Read the CSV file
df = pd.read_csv('faculty_list.csv')

# Keep only the Name_Abbreviated column
df = df[['Name_Abbreviated']]

# Rename the column to name
df = df.rename(columns={'Name_Abbreviated': 'name'})

# Create an app instance and initialize the database
app = create_app()

with app.app_context():
    # Drop and recreate the FacultyMember table
    db.drop_all()  # This will drop all tables, including the FacultyMember table
    db.create_all()  # Recreate all tables
    
    # Populate the FacultyMember table
    db.session.query(FacultyMember).delete()  # Clear the existing data
    for index, row in df.iterrows():
        faculty_member = FacultyMember(name=row['name'])
        db.session.add(faculty_member)
    db.session.commit()

print("Faculty members have been successfully populated.")
