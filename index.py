#####Import libraries
import base64
import re
from time import sleep
import uuid

import numpy as np
import pandas as pd
import streamlit as st
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


import io
import zipfile

##### Make function for save download of zip files without reloading (a streamlit bug)
def download_button_zip(data_frame, download_filename, button_text):
    # Create a BytesIO buffer to hold the zipped data
    zip_buffer = io.BytesIO()
    
    # Write the data frame to a CSV file within the zip buffer
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED, False) as zipf:
        # You can change the file name within the zip here if needed
        zipf.writestr(download_filename, data_frame.to_csv(index=False))
    
    # Seek to the beginning of the zip buffer
    zip_buffer.seek(0)
    
    # Encode the zipped data to base64
    zip_data = base64.b64encode(zip_buffer.read()).decode()
    
    # Create the download link
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

    dl_link = custom_css+f'''
        <a download="{download_filename}.zip" id="{button_id}" href="data:application/zip;base64,{zip_data}">{button_text}</a><br></br>
    '''
    
    return dl_link
#Done with download button for zip file


# Set name for Webpage*******************************
st.set_page_config(
    page_title="Methylation Source File Generator"
)
#Done with setting Name for Webpage


##### Modify config.toml
# Done

##### Load template files
sample_tem_df=pd.read_csv("DataFileTemplate.csv")


##### Set Title************************************************
st.title("Welcome!")
""

##### Make upload button for Methylation data File********************
st.subheader("Methylation Data File")

colA, colB =st.columns(2)
with colA:
    st.write("You may make use of the following template.")
with colB:
    sample_template_btn = download_button(sample_tem_df,'DataFileTemplate.csv',"Data File Template")
    st.markdown(sample_template_btn, unsafe_allow_html=True)

st.caption("The first column should contain the Target ID (cgIDs, etc) followed by their corresponding averages and then the methylation data for each sample.")
#st.caption("Your file should follow the format as provided in the template, which is: ")
#st.caption("Column A: Target ID")
#st.caption("Column B: Average")
#st.caption("Column C, D, E, etc: Sample Methylation Values")

data_files= st.file_uploader("Upload your File here !", type=["csv","excel","xlsx"],key=0)
""
# Done with taking in data file



##### Make sure 450K array Manifest file is included in folder
# Done

##### Submit button *******************
submit =st.button("Start Analysis!")
#Done with submit button
""
loader=st.empty()
status=st.empty()

#/******************* First Part Done ***********************/


# If Submit Button clicked
if submit and data_files is not None:

    status.empty()
    data_df=pd.DataFrame()
    mani_df=pd.DataFrame()



    # Run Loader for activity
    with st.spinner('Reading Methylation Sample Data File'):

        # Load Sample File based on format
        if data_files.type=='text/csv':
            data_df=pd.read_csv(data_files,low_memory=False)
        else:
            data_df=pd.read_excel(data_files)
        
        data_df=data_df.fillna("")

    # Display status: Sample Data File Fetched
    status.info("Methylation Sample Data File Read Successfully")
    
    

    # Run Loader for Manifest File Fetching
    with st.spinner('Fetching Manifest File'):

        # Get Manifest file
        mani_df=pd.read_csv("450K Methylation Manifest file.csv",low_memory=False)
        mani_df=mani_df.fillna("")
        mani_df['UCSC_REFGENE_NAME'] = mani_df['UCSC_REFGENE_NAME'].astype("string")

    # Display status: Manifest Data File Fetched
    status.info("Manifest File Fetched Successfully")

    

    # Run loader for Mapping: Avg-> Manifest 
    with st.spinner('Mapping Data File to Manifest File'):
        # Map Averages to Manifest File
        mani_data_df=pd.merge(mani_df,data_df,how='inner',on = 'TargetID')
        #Sort
        mani_data_df=mani_data_df.sort_values(by=['CHR','MAPINFO'])
        # Download for testing (optional)
        mani_data_btn = download_button_zip(mani_data_df,'ManifestWithInput.csv',"Download Combined Data File (ZIP)")
        st.markdown(mani_data_btn, unsafe_allow_html=True)
    

        
    # Display status: Base File created
    status.success("Files Merged Successfully")


if submit and (data_files is None):
    status.error("Either one or both of Data File and Average File is not uploaded.")


#/*********************Second Part Done**********************************************/