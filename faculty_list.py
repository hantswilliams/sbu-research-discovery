import requests
from bs4 import BeautifulSoup
import pandas as pd

# Fetch the webpage content
url = "https://healthprofessions.stonybrookmedicine.edu/faculty/alpha"
response = requests.get(url)
response.raise_for_status()  # Check if the request was successful

# Parse the HTML content
soup = BeautifulSoup(response.content, 'html.parser')

# Find the div container with class "region region-content"
content_div = soup.find('div', class_='region region-content')

# Find all list items <li> within this div
list_items = content_div.find_all('li')


# Extract data
data = []
for item in list_items:
    link = item.find('a')
    if link:
        name = link.text.strip()
        href = link['href']
        data.append({'Name': name, 'Link': href})

# Convert to DataFrame
df = pd.DataFrame(data)

## print each Name
for name in df['Name']:
    print(name)

## create another name_abbreviated column that only keeps what comes before the first comma
df['Name_Abbreviated'] = df['Name'].str.split(',').str[0]

## create another column just for the last name, and another for the first name based on the new Name_Abbreviated column where we take the first word as the first name and the last word as the last name
df['Last_Name'] = df['Name_Abbreviated'].str.split().str[-1]
df['First_Name'] = df['Name_Abbreviated'].str.split().str[0]

## create another column called terminal degree, where if the name contains "PhD", "Ph.D", "PHD", "DPT", "MD", "DDS", "DVM", "DO", "EdD", "DrPH" we will put YES, otherwise NO
df['Terminal_Degree'] = df['Name'].str.contains('PhD|Ph.D|PHD|DPT|MD|DDS|DVM|DO|EdD|DrPH', case=False, regex=True).map({True: 'YES', False: 'NO'})

## save the data to a CSV file
df.to_csv('faculty_list.csv', index=False)
