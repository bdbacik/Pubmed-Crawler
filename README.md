## Pubmed-Crawler for Covid-19 Research Publications
* Pubmed crawler module that collects paper title, author list, publication time, and abstract from PubMed for a given keyword (i.e., COVID-19 Vaccine) within a specified time window.
* Database module that creates a SQLite database, populates tables, and queries tables for data needed for visualization.
* Visualization module that produces a dashboard summarizing key statistics related to COVID-19 vaccine publications.
* Automation module to automatically run script and update visualization weekly using AWS and cron.


### Code and Resources Used
* Python Version: 3.7
* Packages: biopython, pandas, numpy, sqlite3, plotly, orca, boto3, datetime
* To install requirements: 'pip install -r requirements.txt'
* Automation with AWS and cron: https://www.freecodecamp.org/news/how-to-create-auto-updating-data-visualizations-in-python-with-matplotlib-and-aws/
 

### Module 1: Pubmed Crawler
* Query pubmed for all publications for given search term and time period.  
  * search term used: 'COVID-19 Vaccine'
  * original search time period: 03/01/2020 - 03/01/2021
* Parse results and add data to pandas dataframes
* Note: Now that database is populated with results from 3/1/20-3/1/21, date range has been updated to pull most recent 7 days data.  Program is setup to run every sunday morning and add new results to the database.


### Module 2: Database Module
* Build SQLite database to store pubmed crawler output.
* Create SQL tables to represent data from module 1.
* Query tables to get key statistics and setup for visualization.

### Module 3: Visualization
* Produce a dashboard summarizing key statistics related to Covid-19 vaccine publications
* Save dashboard as .png to local machine
 
![pubmed dashboard](https://github.com/bdbacik/Pubmed-Crawler/blob/main/images/pubmed_dashboard2.png)

### Module 4: Automation
* Create AWS S3 bucket to store dashboard image that we are updating weekly
* Save dashboard .png file from local machine to AWS S3 bucket
* Embed dashboard image on Github pages site using URL of the image from AWS S3 bucket
* Create AWS EC2 instance to schedule our Python script to run automatically every week
* Schedule Python script to run every Sunday morning at 7am using cron: 
  * Create .cron file: 'vim pubmed_crawler.cron'
  * Schedule to run every sunday at 7am: '00 7 * * 7 python3 Pubmed_Crawler.py'
