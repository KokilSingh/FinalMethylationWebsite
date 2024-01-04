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


# Custom Sort for CHR
def chr_sort_key(value):
    try:
        return (float(value), '')
    except ValueError:
        return (float('inf'), str(value))


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
st.subheader("Methylation Data File ")

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
        
        # Warning label for Missing Avg Values
        #st.write(f"{data_df.loc[:,data_df.columns[1]].isna() }")
        if data_df.loc[:,data_df.columns[1]].isna().values.sum()>0 :
            st.warning("Alert:- Missing Average Values Found !")

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
        data_df['TargetID'] = data_df['TargetID'].astype("string")
        for i in range(1,data_df.shape[1]):
            data_df[data_df.columns[i]] = data_df[data_df.columns[i]].apply(pd.to_numeric, downcast='float', errors='coerce')
        mani_data_df=pd.merge(mani_df,data_df,how='inner',on = 'TargetID')

        # Check if certain cgIDs are not recognized by the Manifest File
        if len(data_df['TargetID']) != len(mani_data_df['TargetID']) :
            st.warning('Alert:- TargetIDs in Data File not found in Manifest File !')
            
            ### Thinks that we need to work on : Handling !series_matrix_end

        #Sort and ntype assignment of CHR and MAPINFO
        mani_data_df=mani_data_df.sort_values(by=['CHR', 'MAPINFO'], key=lambda x: x.map(chr_sort_key))
        
        # Download for testing (optional)
        mani_data_btn = download_button_zip(mani_data_df,'ManifestWithInput.csv',"Download Combined Data File (ZIP)")
        st.markdown(mani_data_btn, unsafe_allow_html=True)

        
    

        
    # Display status: Base File created
    status.info("Files Merged Successfully")

    #/*********************Second Part Done**********************************************/

    #######The third part includes Step 2: Delta Beta Generation


    # Initializing global variables 
    sample_list=None        # Stores the list of samples
    avg_idx=1               # Stores the loci of the Average column



    # Run loader for fetching Average and Sample columns
    with st.spinner('Identifying Columns...'):
        #Get average_column index position
        avg_idx=mani_data_df.columns.get_loc('STRAND') + 1 #We find strands loci because if there is any diff in the spelling of Average it will cause an issue
        #Get Sample Names
        sample_list=mani_data_df.columns[avg_idx+1:]
    
    #Display status: Done Identifying which columns are what
    status.info('Done Identifying Columns')


    
    #### Calculating Delta Beta Values

    # Run loader for creating and calculating Delta Beta Values
    with st.spinner('Calculating Delta Beta Values'):

        for col_name in sample_list:
            
            #Get index number of sample column
            index_no = mani_data_df.columns.get_loc(col_name)
            #Insert Delta Beta Column Right Next to Sample 
            mani_data_df.insert(index_no+1,'Delta_Beta_'+col_name,'')
            #Calculate Delta Beta
            mani_data_df.loc[:,'Delta_Beta_'+col_name] =mani_data_df.loc[:,col_name]-mani_data_df.loc[:,'Average']

        # Download button for delta_beta file
        mani_data_df=mani_data_df.fillna("")
        delta_beta_btn = download_button_zip(mani_data_df,'DeltaBeta.csv',"Download Input Data With Delta Beta File (ZIP)")
        st.markdown(delta_beta_btn, unsafe_allow_html=True)
    
    status.success('Delta Beta Values Generated')

    #/*********************Third Part Done**********************************************/
    # Commit to update Streamlit App
    


if submit and (data_files is None):
    status.error("Either one or both of Data File and Average File is not uploaded.")
