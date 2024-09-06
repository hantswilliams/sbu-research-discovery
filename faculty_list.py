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
df['SHP_Website_Terminal_Degree'] = df['Name'].str.contains('PhD|Ph.D|PHD|DPT|MD|DDS|DVM|DO|EdD|DrPH', case=False, regex=True).map({True: 'YES', False: 'NO'})

## create a new column that just keeps the faculty first and last name, but in the format "Last, First"
df['Last_First'] = df['Last_Name'] + ', ' + df['First_Name']

## save the data to a CSV file
df.to_csv('faculty_list.csv', index=False)
df.columns

## Can then bring in the Employees_SR_HR_Employees_Roster.csv 
# which is a tableau generate file that contains all employees at Stony Brook University for SHP

tableau_df = pd.read_csv('Employees__SR_HR_Employees_Roster.csv')

## create a merged table that brings together the employee roster on the left and the faculty list on the right, merging 
## on the left by 'Employee' and on the right by 'Last_First'
merged_df = pd.merge(tableau_df, df, left_on='Employee', right_on='Last_First', how='left')
merged_df.rename(columns={'Name': 'SHP_Website_Name', 'Link': 'SHP_Website_Link'}, inplace=True)
merged_df.drop(columns=['Name_Abbreviated'], inplace=True)

## add any rows from the faculty list that did not match with the employee roster
## these are the faculty that are not in the employee roster
missing_faculty = df[~df['Last_First'].isin(merged_df['Last_First'])]
missing_faculty.rename(columns={'Name': 'SHP_Website_Name', 'Link': 'SHP_Website_Link'}, inplace=True)
missing_faculty = missing_faculty[['SHP_Website_Name', 'SHP_Website_Link', 'Last_First', 'SHP_Website_Terminal_Degree']]
missing_faculty = missing_faculty.to_dict(orient='records')

## loop through the missing faculty and add them to the merged_df
for faculty in missing_faculty:
    merged_df = merged_df._append(faculty, ignore_index=True)

## loop through each row, and if Last_First is null, then set the value to what Employee is
for index, row in merged_df.iterrows():
    if pd.isnull(row['Last_First']):
        merged_df.at[index, 'Last_First'] = row['Employee']

## make the Last_First column the first column
cols = merged_df.columns.tolist()
cols = cols[-1:] + cols[:-1]
merged_df = merged_df[cols]

## then sort the dataframe by Last_First
merged_df.sort_values(by='Last_First', inplace=True)

## change Last_Name and First_Name to SHP_Website_Last_Name and SHP_Website_First_Name
merged_df.rename(columns={'Last_Name': 'SHP_Website_Last_Name', 'First_Name': 'SHP_Website_First_Name'}, inplace=True)

## if column name does not start with SHP_, or is Last_First, add in Tableau_ to the beginning of the column name
for col in merged_df.columns:
    if not col.startswith('SHP_') and col != 'Last_First':
        merged_df.rename(columns={col: 'Tableau_' + col}, inplace=True)

## save the merged data to a CSV file
merged_df.to_csv('merged_faculty_list.csv', index=False)