#!/usr/bin/env python
# coding: utf-8

# ## Module 1
# 
# Pubmed crawler module tht can collect paper title, author list, publication time, and abstract from PUBMED for a given keyword (i.e., SARS-CoV-2) within a pre-specified time window (that is, 01/01/2020 - 12/31/2020). Save the output to csv for further analysis later.

# In[230]:


import pandas as pd
import datetime
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sqlite3
import numpy as np


# In[127]:


#run once after connecting to install and import tools
#!pip install biopython


# In[240]:


#write Pubmed search query
from Bio import Entrez

def search(query):
    Entrez.email = 'bdbacik@gmail.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='1000',
                            retmode='xml',
                            datetype='pdat',
                            mindate='2020/03/01',
                            maxdate='2020/12/31',
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'isela.delacerda@gmail.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results


if __name__ == '__main__':
    results = search('SARS-CoV-2')
    id_list = results['IdList']
    papers = fetch_details(id_list)
#print(papers)


# In[241]:


#pull relevant fields from query into pandas dataframe

searchoutput = {"Title":[], "DateCompleted":[], "PublicationDate": [], "Authors": [], "Abstract": []}
for i, paper in enumerate(papers['PubmedArticle']): 
    try:
        Title = paper['MedlineCitation']['Article']['ArticleTitle']
        DateCompleted = paper['MedlineCitation']['Article']['ArticleDate'] #['Journal']['JournalIssue']['PubDate']['Month']
        PublicationDate = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month']
        Authors = paper['MedlineCitation']['Article']['AuthorList']
        Abstract = paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
    except KeyError as e:
        continue
    searchoutput["Title"].append(Title)
    searchoutput['DateCompleted'].append(DateCompleted)
    searchoutput["PublicationDate"].append(PublicationDate)
    searchoutput["Authors"].append(Authors)
    searchoutput["Abstract"].append(str(Abstract))

df = pd.DataFrame(searchoutput)

df.head()


# In[242]:


# Clean Date Completed and Publication Date Columns
df["PublicationDate"] = df["PublicationDate"].replace(["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],["01","02","03","04","05","06","07","08","09","10","11","12"])
df['DateCompleted'] = df.explode('DateCompleted')['DateCompleted'].apply(pd.Series).loc[:,['Month']]
print(df["PublicationDate"].value_counts())
print(df['DateCompleted'].value_counts())


# In[243]:


#create new dataframe with columns for first and last name of each author of each paper
df2 = df.explode('Authors')['Authors'].apply(pd.Series).loc[:,['LastName','ForeName','AffiliationInfo']]
df2["Author"] = df2["ForeName"].str.cat(df2["LastName"], sep=" ")
#merge original df and df2 with author names
df3 = pd.merge(df,df2,left_index=True,right_index=True)
df3 = df3[["Title",'DateCompleted',"PublicationDate", "Abstract", "Author", 'AffiliationInfo']]
df3.head()


# In[244]:


#save new dataframe with author names to csv
df3.to_csv(r'trialrun2.csv', index=False)


# # Module 2
# Database module that imports the CSV file and build SQLite database. Then query SQL table to find authors with the most publications with keyword 'SARS-CoV-2' in 2020 and pull all publications for a given author.

# In[245]:


conn = sqlite3.connect('pubmed_crawl.db')  # create a new database or connect to database if already exists
c = conn.cursor() # create connection object


# In[246]:


# Create table if it doesn't yet exist
try:
  c.execute('''CREATE TABLE pubmed_table2
              ([Title] text, [DateCompleted] text, [PublicationDate] text, [Abstract] text, [Author] text, [AffiliationInfo] text)''')
  conn.commit()
  print('SQL table created successfully!')

except:
  print('SQL table already exists!')


# In[214]:


df = pd.read_csv('trialrun2.csv') #read csv file created in Task 1 (could also access pandas df directly)
df.to_sql('pubmed_table2', conn, if_exists='replace', index=False) #add data to sql table


# In[248]:


#SQL query to find top 10 authors by number of 'SARS-CoV-2' publications in 2020
c.execute('SELECT DISTINCT Author, count() OVER(PARTITION BY Author) as Count FROM pubmed_table2 WHERE Author IS NOT NULL ORDER BY Count DESC LIMIT 10').fetchall()


# In[266]:


#query sql table to return all publications for a given author using input from previous step
def count_titles():
    #prompt user to input author name they wish to search for
    full_name = input('Enter author first name last name(separated by space): ')
    count = c.execute('SELECT count(*) FROM pubmed_table2 WHERE Author = ?', (full_name,)).fetchone()
    titles = c.execute('SELECT Title FROM pubmed_table2 WHERE Author = ?', (full_name,)).fetchall()
    affiliation = c.execute('SELECT AffiliationInfo FROM pubmed_table2 WHERE Author = ?', (full_name,)).fetchone()
    print(full_name, "has" , count , "articles." ,"\n")
    print(full_name, 'is affiliated with: ', affiliation, '\n')
    for row in titles:
        print(row)
    
count_titles()


# # Module 3 
# Read CSV file and visualize trend in number of publications by month

# In[263]:


#Take a subset of the titles only once
publications = pd.read_csv('trialrun2.csv')
publications = publications.drop_duplicates(subset=['Title'])
pubs_by_month = publications["DateCompleted"].value_counts()
#get summary statistics for publications by month
pubs_by_month_summary = pubs_by_month.describe()


# In[270]:


#create dashboard to visualize publications by month, trend over time, and summary statistics

#define subplot figure contents
fig = make_subplots(
    rows=3, cols=1,
    subplot_titles=("Bar Chart of # Publications by Month", "Boxplot of # Publications by Month", "Summary Statistics of # Publications by Month"),
    shared_xaxes=False,
    vertical_spacing=0.1,
    specs=[[{"type": "bar"}],
           [{"type": "box"}],
           [{"type": "table"}]])

#create histogram in first row
fig.add_trace(go.Bar(x=pubs_by_month.index, y=pubs_by_month.values,text=pubs_by_month.values,
                     textposition='auto',name='Bar Chart of # Publications by Month'), row=1, col=1 )
fig.update_layout(xaxis = dict(tickmode = 'linear',tick0 = 1,dtick = 1 ))

#create boxplot in second row
fig.add_trace(go.Box(x=pubs_by_month), row=2,col=1 )

#create summary statistics table in thirs row
fig.add_trace(go.Table(header=dict(values=['Stat', 'Value']), 
                       cells=dict(values=[pubs_by_month_summary.index, pubs_by_month_summary.values])), row=3, col=1)

#define layout
fig.update_layout(
    height=800,
    showlegend=False,
    title_text="Dashboard of Publications by Month")

fig.show()


# In[ ]:




