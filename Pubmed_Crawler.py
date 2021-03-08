import pandas as pd
import datetime
from datetime import date
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sqlite3
import numpy as np
import boto3

# # Module 1 Pubmed Crawler ####################################################################################
# Pubmed crawler module tht can collect paper title, author list, publication time, and abstract from PUBMED for a given keyword (i.e., COVID-19 Vaccine) within a specified time window.

#get date ranges to use in pubmed query
today = date.today()
cur_date = today.strftime("%Y/%m/%d")

last_week = (today-datetime.timedelta(days=7)).strftime("%Y/%m/%d")
print(cur_date, last_week)


#write Pubmed search query
from Bio import Entrez

#returns list of article ids that contain keyword 
def search(query):
    Entrez.email = 'bdbacik@gmail.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='100000',
                            retmode='xml',
                            datetype='pdat',
                            mindate=last_week,
                            maxdate=cur_date,
                            term=query)
    search_results = Entrez.read(handle)
    return search_results

#users article ids as input and returns dictionary with article details
def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'bdbacik@gmail.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    fetch_results = Entrez.read(handle)
    return fetch_results


if __name__ == '__main__':
    search_results = search('COVID-19 vaccine')
    id_list = search_results['IdList']
    fetch_results = fetch_details(id_list)


print('length of esearch output: ', len([item for item in search_results['IdList']]))
print('length of efetch output: ', len([key for key in fetch_results['PubmedArticle']]))


#pull relevant fields from query into pandas dataframe
searchoutput = {"Title":[], "Keywords":[], "PublicationDate": [], "Authors": [], 
                "Abstract": [], 'Country':[]}
for i, paper in enumerate(fetch_results['PubmedArticle']): 
    try:
        Title = paper['MedlineCitation']['Article']['ArticleTitle']
        Keywords = paper['MedlineCitation']['KeywordList']
        PublicationDate = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']#['Month']
        Authors = paper['MedlineCitation']['Article']['AuthorList']
        Abstract = paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
        Country = paper['MedlineCitation']['MedlineJournalInfo']['Country']
    except KeyError as e:
        continue
    searchoutput["Title"].append(Title)
    searchoutput["Keywords"].append(Keywords)
    searchoutput["PublicationDate"].append(PublicationDate)
    searchoutput["Authors"].append(Authors)
    searchoutput["Abstract"].append(str(Abstract))
    searchoutput["Country"].append(Country)

df = pd.DataFrame(searchoutput)

#get publication date in YYYYMM 
df['Pub_Year'] = [str(d.get('Year')) for d in df['PublicationDate']]
df['Pub_Month'] = [str(d.get('Month')) for d in df['PublicationDate']]
df["Pub_Month"] = df["Pub_Month"].replace(["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],["01","02","03","04","05","06","07","08","09","10","11","12"])
df['Pub_Date'] = pd.to_datetime(df.Pub_Year+'-'+df.Pub_Month, errors='ignore')
df = df.drop(columns=['PublicationDate','Pub_Year','Pub_Month'])

#make separate df with multi value attribute columns 
keywords = df.explode('Keywords').explode('Keywords').iloc[:,:2]
keywords['Keywords'] = keywords['Keywords'].str.lower()


#create new dataframe with columns for first and last name of each author of each paper
authors = df.explode('Authors')['Authors'].apply(pd.Series).loc[:,['LastName','ForeName','AffiliationInfo']]
authors["Author"] = authors["ForeName"].str.cat(authors["LastName"], sep=" ")
#merge original df and df2 with author names
df3 = pd.merge(df,authors,left_index=True,right_index=True)
df3 = df3[["Title", "Author", 'AffiliationInfo']]
df3['AffiliationInfo'] = df3['AffiliationInfo'].astype(str)
df3['AffiliationInfo'] = df3.apply(lambda x: x['AffiliationInfo'].split('Affiliation')[-1], axis=1).str[3:-2]

#remove duplicate columns from main df
df = df.drop(columns=['Keywords','Authors'])


# # Module 2 - Database Development ####################################################################
# 1. Builds SQLite database to store our pubmed crawler output. 
# 2. Create SQL tables to represent data from module 1.
# 3. Query tables to get key statistics and setup for visualization.

# ## Create database
conn = sqlite3.connect('pubmed.db')  # create a new database or connect to database if already exists
c = conn.cursor() # create connection object


# ## Create tables and add data to them
# We have three tables in the database:
# 1. pub_info - this contains the main information about the publication, including author, abstract, country published in, and publication date
# 2. authors - this table contains author name and affiliaition, with title as foreign key
# 3. keywords - this table contains keywords for each article, with title as foreign key

# Create pub_info table
try:
  c.execute('''CREATE TABLE pub_info
              ([Title] text, [Abstract] text, [Country] text, [Pub_Date] text)''')
  conn.commit()
  print('SQL table created successfully!')

except:
  print('SQL table already exists!')

#add data to sql table from our pandas dataframe
df.to_sql('pub_info', conn, if_exists='append', index=False)

#query database to get the number of publications by country
pubs_by_country = c.execute('''SELECT DISTINCT Country, count() OVER(PARTITION BY Country) as Num_Publications
                            FROM pub_info 
                            ORDER BY Num_Publications DESC 
                            ''').fetchall()

#query database to get the number of publications by month
pubs_by_month = c.execute('''SELECT Pub_date, count(Pub_Date) as Pubs_by_month
                            FROM pub_info 
                            GROUP BY Pub_date
                            ORDER BY Pubs_by_month DESC 
                            ''').fetchall()

# Create authors table
try:
  c.execute('''CREATE TABLE authors
              ([Title] text, [Author] text, [AffiliationInfo] text)''')
  conn.commit()
  print('SQL table created successfully!')

except:
  print('SQL table already exists!')


#add data to sql table from our pandas dataframe
df3.to_sql('authors', conn, if_exists='append', index=False)


#query database to find top 10 authors by number of publications
top_authors = c.execute('''SELECT DISTINCT Author, count(Author) as Count, AffiliationInfo
                        FROM authors 
                        WHERE Author IS NOT NULL
                        GROUP BY Author
                        ORDER BY Count DESC 
                        LIMIT 10''').fetchall()

# Create keywords table
try:
  c.execute('''CREATE TABLE keywords
              ([Title] text, [Keywords] text)''')
  conn.commit()
  print('SQL table created successfully!')

except:
  print('SQL table already exists!')


#add data to sql table from our pandas dataframe
keywords.to_sql('keywords', conn, if_exists='append', index=False)


#query database to find top 10 keywords
top_keywords = c.execute('''SELECT Keywords, count(Keywords) as Count
                        FROM keywords
                        GROUP BY Keywords
                        ORDER BY Count DESC 
                        LIMIT 20''').fetchall()
top_keywords[0:20]


# # Module 3 - Visualization ################################################################################
# Use SQL query output to produce a dashboard summarizing key statistics related to publications

# ## 1. Publications by month
#convert SQL query output to pandas df for visualization
df_pubs_by_month = pd.DataFrame(pubs_by_month, columns=['Pub_Date', 'Count'])
df_pubs_by_month['Pub_Date'] = df_pubs_by_month['Pub_Date'].astype(str)
df_pubs_by_month.head()


# ## 2. Publications by country
#convert SQL query output to pandas df for visualization
df_pubs_by_country = pd.DataFrame(pubs_by_country, columns=['Country', 'Num_Publications'])
total_pubs = df_pubs_by_country['Num_Publications'].sum()
print(total_pubs)
df_pubs_by_country['Percent_Pubs'] = df_pubs_by_country['Num_Publications']/total_pubs
df_pubs_by_country.head()


# ## 3. Publications by author
df_top_authors = pd.DataFrame(top_authors,columns=['Author', 'Count', 'Affiliation'])
df_top_authors.head()


# ## Create Dashboard
# using plotly

#create dashboard to visualize publications by month, trend over time, and summary statistics

#define subplot figure contents
fig = make_subplots(
    rows=3, cols=1,
    subplot_titles=("Publications by Month", "Publications by Country", "Top Authors by Number of Publications"),
    shared_xaxes=False,
    vertical_spacing=0.1,
    specs=[[{"type": "bar"}],
           [{"type": "bar"}],
           [{"type": "table"}]])

#create histogram in first row
fig.add_trace(go.Bar(x=df_pubs_by_month.loc[:12,'Pub_Date'], y=df_pubs_by_month.loc[:12,'Count'],text=df_pubs_by_month.loc[:12,'Count'],
                     textposition='auto',name='Publications by Month'), row=1, col=1 )
fig.update_xaxes(
    dtick=1,
    tick0=1,
    tickformat="%m\n%Y", type='category')
#fig.update_layout(xaxis = dict(tickformat="%b\n%Y" ))

#create boxplot in second row
fig.add_trace(go.Bar(x=df_pubs_by_country.loc[:10,'Country'], y=df_pubs_by_country.loc[:10,'Num_Publications'],text=df_pubs_by_country.loc[:10,'Num_Publications'],
                     textposition='auto',name='Publications by Country'), row=2, col=1 )
fig.update_layout(xaxis = dict(tickmode = 'linear',tick0 = 1,dtick = 1 ))

#create summary statistics table in thirs row
fig.add_trace(go.Table(header=dict(values=['Author', 'Publications', 'Affiliations']), 
                       cells=dict(values=[df_top_authors.loc[:2,'Author'], df_top_authors.loc[:2,'Count'], 
                                          df_top_authors.loc[:2,'Affiliation']])), row=3, col=1)

#define layout
fig.update_layout(
    height=1000,
    showlegend=False,
    title_text="PubMed Publications for Research Related to Covid-19 Vaccine")

fig.show()

#save image to local machine
fig.write_image("pubmed_dashboard2.png", scale=1, width=1200, height=1200)

#upload image to AWS
s3 = boto3.resource('s3')
s3.meta.client.upload_file('pubmed_dashboard2.png', 'pubmedcrawler', 'pubmed_dashboard2.png', 
                           ExtraArgs={'ACL':'public-read'})




