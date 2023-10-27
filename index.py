#####Import libraries
import gc
import base64
import re
from time import sleep
import uuid
from array import array

import numpy as np
import pandas as pd
import streamlit as st
from scipy.stats.distributions import chi2
from stqdm import stqdm
#Done importing relevant libraries


##### Make function for save download of files without reloading (a streamlit bug)
def download_button(object_to_download, download_filename, button_text):
    
    # This function has been created as an alternative to streamlit's st.download_button
    # because of bug number 3832 

    object_to_download = object_to_download.to_csv(index=None)
    
    try:
        b64 = base64.b64encode(object_to_download.encode()).decode()
    
    except AttributeError as e:
        b64 = base64.b64encode(object_to_download).decode()

    button_uuid = str(uuid.uuid4()).replace('-', '')
    button_id = re.sub('\d+', '', button_uuid)

    custom_css = f""" 
        <style>
            #{button_id} {{
                background-color: rgb(19,23, 32);
                color: inherit;
                padding: 0.5rem 0.75rem;
                position: relative;
                text-decoration: none;
                border-radius: 0.25rem;
                border-width: 1px;
                border-style: solid;
                border-color: rgb(250, 250, 250,0.2);
                border-image: initial;
            }} 
            #{button_id}:hover {{
                border-color: rgb(246, 51, 102);
                color: rgb(246, 51, 102);
            }}
            #{button_id}:active {{
                box-shadow: none;
                background-color: rgb(246, 51, 102);
                color: white;
                }}
        </style> """

    dl_link = custom_css + f'<a download="{download_filename}" id="{button_id}" href="data:file/txt;base64,{b64}">{button_text}</a><br></br>'

    return dl_link
#Done with download button function


# Set name for Webpage*******************************
st.set_page_config(
    page_title="Methylation Source File Generator"
)
#Done with setting Name for Webpage


##### Modify config.toml
# Done

##### Load template files
sample_tem_df=pd.read_csv("SampleDataTemplate.csv")
avg_tem_df=pd.read_csv("AverageDataTemplate.csv")


##### Set Title************************************************
st.title("Welcome!")
""

##### Make upload button for Methylation data File********************
st.subheader("Methylation Data File")

colA, colB =st.columns(2)
with colA:
    st.write("You may make use of the following template.")
with colB:
    sample_template_btn = download_button(sample_tem_df,'SamplesDataTemplate.csv',"Samples Data Template")
    st.markdown(sample_template_btn, unsafe_allow_html=True)

st.caption("The first column should contain the Target ID (cgIDs, etc) followed by the data for each sample.")

data_files= st.file_uploader("Upload your File here !", type=["csv","excel","xlsx"],key=0)
""
# Done with taking in data file


##### Make Template for Average Beta Values
st.subheader("Average Data File")
# Make Download Button for Avg Beta Values Template ********************
colA, colB =st.columns(2)
with colA:
    st.write("You may make use of the following template.")
with colB:
    sample_template_btn = download_button(avg_tem_df,'AverageDataTemplate.csv',"Average Data Template")
    st.markdown(sample_template_btn, unsafe_allow_html=True)

st.caption("The first column should contain the Target ID (cgIDs, etc) followed by the Average for each ID.")

# Make Upload Button For Avg Values that follow template *****************
avg_files= st.file_uploader("Upload your File here !", type=["csv","excel","xlsx"], key=1)
""
#Done with taking in average file


##### Make sure 450K array Manifest file is included in folder
# Done

##### Submit button *******************
submit =st.button("Start Analysis!")
#Done with submit button
""
loader=st.empty()


#/******************* First Part Done ***********************/
