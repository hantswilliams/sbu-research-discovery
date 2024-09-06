import pandas as pd 
from webofscience.request import websci_byAuthor

faculty = pd.read_csv('merged_faculty_list.csv')

articles_list = []
notfound_list = []

for index, row in faculty.iterrows():
    if row['SHP_Website_Terminal_Degree'] == 'YES':
        facultyName = row['Last_First']
        last_name = row['Last_First'].split()[0]
        last_name = last_name.replace(',', '')
        first_name = row['Last_First'].split()[1]
        print(f'working on {last_name}, {first_name} ...')
        try:
            articles = websci_byAuthor(last_name, first_name)
            print(f'Found {len(articles)} articles for {first_name} {last_name}')
            articles_list.append(articles)
        except Exception as e:
            print(f'Error or not articles for {first_name} {last_name}: {e}')
            notfound_list.append(facultyName)
            continue
    else:
        print(f'{row["Last_First"]} does not have a terminal degree')

len(articles_list)

## loop through the articles_list and flatten it
## then convert to a pandas dataframe

df_list = []

for articles in articles_list:
    df = pd.json_normalize(articles)
    df_list.append(df)

## flatten the list of dataframes
df = pd.concat(df_list, ignore_index=True)
        
## save the data to a CSV file
df.to_csv('./webofscience/webofscience_articles_loop.csv', index=False)

