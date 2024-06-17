#####Import libraries
import base64
import re
from time import sleep
import uuid

import numpy as np
import pandas as pd
import streamlit as st
from scipy import stats
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



##### Adding Funtions pertaining to Fourth part


# Identifies hypomethylated consecutive ID blocks
def identify_consecutive_blocks_hypo(col_delta_beta, df, hypo_parameter):
    
    # hypoparameter shd always be <0 (is typically -0.1)
    # col_delta_beta holds name of delta beta column
    # df - The dataframe being used
    
    # Turn Delta Beta to np array
    dif=np.array(df[col_delta_beta].iloc[:])
    # Create filter to identify hypomethylated IDs
    in1=np.where(dif<=hypo_parameter)
    dif[in1[0]]=1

    bf=[]
    i=-1
    
    maxv=np.where(dif==1)[0]
    # Check if there are any occurrences of 1
    if len(maxv) > 0:
        # Retrieve the last index where dif is equal to 1
        maxv = maxv[-1]
    else:
        # Handle the case where there are no occurrences of 1
        maxv = -1  # or any other appropriate value

    while i<maxv:
        bn=[]
        i=i+1
        if dif[i]==1 and i<=maxv:
            bn.append(i)
            chrp=df['CHR'].iloc[i]
            i=i+1
            while i<maxv and dif[i]==1 and chrp==df['CHR'].iloc[i]:
                bn.append(i)
                i=i+1
            
        if i==maxv:
            bn.append(maxv)
        
        if len(bn)>=3:
            bf.append(bn)
    
    return bf


# Identifies hypermethylated consecutive ID blocks
def identify_consecutive_blocks_hyper(col_delta_beta, df, hyper_parameter):
    
    # hyperparameter shd always be >0 (is typically 0.1)
    # col_delta_beta holds name of delta beta column
    # df - The dataframe being used
    
    # Turn Delta Beta to np array
    dif=np.array(df[col_delta_beta].iloc[:])
    
    # Create filter to identify hypomethylated IDs
    in1=np.where(dif>=hyper_parameter)
    dif[in1[0]]=1

    bf=[]
    i=-1
    maxv=np.where(dif==1)[0]
    # Check if there are any occurrences of 1
    if len(maxv) > 0:
        # Retrieve the last index where dif is equal to 1
        maxv = maxv[-1]
    else:
        # Handle the case where there are no occurrences of 1
        maxv = -1  # or any other appropriate value

    while i<maxv:
        bn=[]
        i=i+1
        if dif[i]==1 and i<=maxv:
            bn.append(i)
            chrp=df['CHR'].iloc[i]
            i=i+1
            while i<maxv and dif[i]==1 and chrp==df['CHR'].iloc[i]:
                bn.append(i)
                i=i+1
            
        if i==maxv:
            bn.append(maxv)
        
        if len(bn)>=3:
            bf.append(bn)
    
    return bf


# Function to calculate means of averages, sample values, their diff, p-value and organise them for final df
def organise_blocks(bf,df,avg_idx,sample_col_name):

    fd=[]

    for t in bf:
        f=[]
        
        # Add Sample
        cna =sample_col_name.split('.')
        f.append(cna[0])

        #First blocks chromosome
        chv_first=df['CHR'].iloc[t[0]]
        #Last blocks chromosome
        chv_last =df['CHR'].iloc[t[-1]]

        if chv_first==chv_last:
            cvff=str(chv_first)+':'+str(int(df['MAPINFO'].iloc[t[0]]))+'-'+str(int(df['MAPINFO'].iloc[t[-1]]))
        
        f.append(cvff)
        f.append(df['MAPINFO'].iloc[t[-1]]-df['MAPINFO'].iloc[t[0]])
        f.append(len(t))
        f.append((df['MAPINFO'].iloc[t[-1]]-df['MAPINFO'].iloc[t[0]])/len(t))

        for j in range(3,avg_idx):
            f.append(df[df.columns[j]].iloc[t[0]])
        
        # Mean Controls_average
        f.append(np.mean(df[df.columns[avg_idx]].iloc[t]))
        # Mean Sample Beta
        f.append(np.mean(df[sample_col_name].iloc[t]))
        # Methylation Diff
        f.append(np.mean(df[sample_col_name].iloc[t])-np.mean(df[df.columns[avg_idx]].iloc[t]))
        t_value,p_value=stats.ttest_rel(df[df.columns[avg_idx]].iloc[t],df[sample_col_name].iloc[t])
        # P-value
        f.append(p_value)
        fd.append(f)
    return fd




##### Adding Funtions pertaining to Fourth part


# Identifies hypomethylated consecutive ID blocks
def identify_consecutive_blocks_hypo(col_delta_beta, df, hypo_parameter):
    
    # hypoparameter shd always be <0 (is typically -0.1)
    # col_delta_beta holds name of delta beta column
    # df - The dataframe being used
    
    # Turn Delta Beta to np array
    dif=np.array(df[col_delta_beta].iloc[:])
    # Create filter to identify hypomethylated IDs
    in1=np.where(dif<=hypo_parameter)
    dif[in1[0]]=1

    bf=[]
    i=-1
    
    maxv=np.where(dif==1)[0]
    # Check if there are any occurrences of 1
    if len(maxv) > 0:
        # Retrieve the last index where dif is equal to 1
        maxv = maxv[-1]
    else:
        # Handle the case where there are no occurrences of 1
        maxv = -1  # or any other appropriate value

    while i<maxv:
        bn=[]
        i=i+1
        if dif[i]==1 and i<=maxv:
            bn.append(i)
            chrp=df['CHR'].iloc[i]
            i=i+1
            while i<maxv and dif[i]==1 and chrp==df['CHR'].iloc[i]:
                bn.append(i)
                i=i+1
            
        if i==maxv:
            bn.append(maxv)
        
        if len(bn)>=3:
            bf.append(bn)
    
    return bf


# Identifies hypermethylated consecutive ID blocks
def identify_consecutive_blocks_hyper(col_delta_beta, df, hyper_parameter):
    
    # hyperparameter shd always be >0 (is typically 0.1)
    # col_delta_beta holds name of delta beta column
    # df - The dataframe being used
    
    # Turn Delta Beta to np array
    dif=np.array(df[col_delta_beta].iloc[:])
    
    # Create filter to identify hypomethylated IDs
    in1=np.where(dif>=hyper_parameter)
    dif[in1[0]]=1

    bf=[]
    i=-1
    maxv=np.where(dif==1)[0]
    # Check if there are any occurrences of 1
    if len(maxv) > 0:
        # Retrieve the last index where dif is equal to 1
        maxv = maxv[-1]
    else:
        # Handle the case where there are no occurrences of 1
        maxv = -1  # or any other appropriate value

    while i<maxv:
        bn=[]
        i=i+1
        if dif[i]==1 and i<=maxv:
            bn.append(i)
            chrp=df['CHR'].iloc[i]
            i=i+1
            while i<maxv and dif[i]==1 and chrp==df['CHR'].iloc[i]:
                bn.append(i)
                i=i+1
            
        if i==maxv:
            bn.append(maxv)
        
        if len(bn)>=3:
            bf.append(bn)
    
    return bf




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
st.title("Welcome EPIC ARRAY V2!")
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


st.subheader("Parameters List ")
##### Get hypo and hyper methylation values**************************
colC,colD =st.columns(2)
with colC:
    st.write("Degree of Hyper/Hypomethylation selection")
    st.caption("If you enter 0.1, then only those instances having delta_beta values >=0.1 or <=-0.1, will be checked for consecutive hypermethylation.")
with colD:
    hyper_param=st.number_input("Label",min_value=0.0,max_value=1.0,value=1e-1,format="%.3f",key='hyper')

hypo_param = -hyper_param

##### Get Maximum Ratio value *****************************************
""
agree_row = st.checkbox('**Do you want to filter the DMRs selected by setting a maximum ratio value?**')

colE,colF =st.columns(2)
with colE:
    st.write("Maximum Ratio value ")
    st.caption("Ratio value attributes to the maximum value (Length of DMR/No. of Probes) that should be used to filter DMRs")
with colF:
    r_max =st.number_input("Label2",min_value=0,value=182,key='r_param',disabled=not agree_row)



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
    with st.spinner('Fetching EPIC V2 Manifest File'):

        # Get Manifest file
        mani_df=pd.read_csv("EPIC v2 Manifest File.csv",low_memory=False)
        mani_df=mani_df.fillna("")
        mani_df['UCSC_RefGene_Name'] = mani_df['UCSC_RefGene_Name'].astype("string")
        mani_df['UCSC_RefGene_Accession'] = mani_df['UCSC_RefGene_Accession'].astype("string")
        mani_df['UCSC_RefGene_Group'] = mani_df['UCSC_RefGene_Group'].astype("string")


    # Display status: Manifest Data File Fetched
    status.info("Manifest File Fetched Successfully")

    

    # Run loader for Mapping: Avg-> Manifest 
    with st.spinner('Mapping Data File to Manifest File'):

        # Map Averages to Manifest File
        #data_df['TargetID'] = data_df['TargetID'].astype("string")
        
        for i in range(1,data_df.shape[1]):
            data_df[data_df.columns[i]] = data_df[data_df.columns[i]].apply(pd.to_numeric, downcast='float', errors='coerce')
        
        
        #mani_data_df=pd.merge(mani_df,data_df,how='inner',on = 'TargetID')
        #Alternate strategy
        #print(mani_df['TargetID'].dtype)
        #print(data_df['TargetID'].dtype)
        mani_df['TargetID'] = mani_df['TargetID'].astype(str)
        data_df['TargetID'] = data_df['TargetID'].astype(str) 
        mani_df['TargetID'] = mani_df['TargetID'].str.strip()
        data_df['TargetID'] = data_df['TargetID'].str.strip()

        #common_ids = set(mani_df['TargetID']).intersection(set(data_df['TargetID']))
        #print(f"Number of common TargetIDs: {len(common_ids)}")
        
        # Print some common TargetIDs
        #print(list(common_ids)[:10])


        mani_df.set_index('TargetID', inplace=True)
        data_df.set_index('TargetID', inplace=True)

        # Perform the merge
        mani_data_df = mani_df.join(data_df, how='inner')

        # Reset the index if needed
        mani_data_df.reset_index(inplace=True)
        data_df.reset_index(inplace=True)

        # Check if certain cgIDs are not recognized by the Manifest File
        if len(data_df['TargetID']) != len(mani_data_df['TargetID']) :
            diff=len(data_df['TargetID'])-len(mani_data_df['TargetID'])
            st.warning(f'Alert: {diff} TargetIDs in Data File not found in Manifest File !')
            
            ### Thinks that we need to work on : Handling !series_matrix_end

        #Sort and ntype assignment of CHR and MAPINFO
        #mani_data_df=mani_data_df.sort_values(by=['CHR', 'MAPINFO'], key=lambda x: x.map(chr_sort_key))
        
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
        avg_idx=mani_data_df.columns.get_loc('Strand_CO') + 1 #We find strands loci because if there is any diff in the spelling of Average it will cause an issue
        #Get Sample Names
        sample_list=mani_data_df.columns[avg_idx+1:]
    
    #Display status: Done Identifying which columns are what
    status.info('Done Identifying Columns')


    
    #### Calculating Delta Beta Values*********************************************************
    
    # Run loader for creating and calculating Delta Beta Values
    with st.spinner('Calculating Delta Beta Values'):

        for i in range(17,mani_data_df.shape[1]):
            mani_data_df[mani_data_df.columns[i]] = mani_data_df[mani_data_df.columns[i]].apply(pd.to_numeric, downcast='float', errors='coerce')

        for col_name in sample_list:
            
            #Get index number of sample column
            index_no = mani_data_df.columns.get_loc(col_name)
            #Insert Delta Beta Column Right Next to Sample 
            mani_data_df.insert(index_no+1,'Delta_Beta_'+col_name,'')
            #Calculate Delta Beta
            mani_data_df.loc[:,'Delta_Beta_'+col_name] =mani_data_df.loc[:,col_name]-mani_data_df.loc[:,mani_data_df.columns[avg_idx]]

        # Download button for delta_beta file
        mani_data_df=mani_data_df.fillna("")
        mani_data_df=mani_data_df.sort_values(by=['CHR','MAPINFO'], key=lambda x: x.map(chr_sort_key))
        delta_beta_btn = download_button_zip(mani_data_df,'DeltaBeta.csv',"Download Input Data With Delta Beta File (ZIP)")
        st.markdown(delta_beta_btn, unsafe_allow_html=True)
    
    status.info('Delta Beta Values Generated')



    #/*********************Third Part Done**********************************************/
    
    # The Fourth part deals with some hard core calculations

    # First get hypo and hyperparameter values through form
    # ------> It is already stored in their respective holders

    # Importing stats is done
    # Add functions identify_consecutive_blocks_hypo, identify_consecutive_blocks_hyper and organise_blocks
    # Done (Check top)

    # Run loader to prepare file for next step
    with st.spinner("Setting Things Ready Before Next Step"):
        
        #Initializing Dataframe to store new results
        col=['Sample','Coordinate','Len','Number of Probes','Ratio',mani_data_df.columns[3],mani_data_df.columns[4],mani_data_df.columns[5],mani_data_df.columns[6],mani_data_df.columns[7],mani_data_df.columns[8],mani_data_df.columns[9],mani_data_df.columns[10],mani_data_df.columns[11],mani_data_df.columns[12],mani_data_df.columns[13],mani_data_df.columns[14],mani_data_df.columns[15],mani_data_df.columns[16],'Mean Controls_Average','Mean Sample Beta Value','Methylation Difference','P-value']
        Final_df=pd.DataFrame(columns=col)

        # To handle NaN values
        for i in range(17,mani_data_df.shape[1]):
            mani_data_df[mani_data_df.columns[i]] = mani_data_df[mani_data_df.columns[i]].apply(pd.to_numeric, downcast='float', errors='coerce')

    # Get Progress bar ready
    #pbar = stqdm("Identifying DMRs",total=len(sample_list))

    for col_name in stqdm(sample_list, desc="Identifying DMRs"):
        #print(col_name+"\n")

        #First Store Hypomethylated blocks identified 
        hypo_blocks =identify_consecutive_blocks_hypo('Delta_Beta_'+col_name,mani_data_df,hypo_param)
        #if len(hypo_blocks)>0:
        #    print(hypo_blocks)

        #Second Organise the hypoblocks
        first_block=organise_blocks(hypo_blocks,mani_data_df,17,col_name)

        #Third Store Hypermethylated blocks identified
        hyper_blocks=identify_consecutive_blocks_hyper('Delta_Beta_'+col_name,mani_data_df,hyper_param)
        #if len(hyper_blocks)>0:
        #    print(hypo_blocks)
        
        # Fourth Organise the Hyperblocks
        second_block=organise_blocks(hyper_blocks,mani_data_df,17,col_name)

        col_consecutive = first_block + second_block
        col_df =pd.DataFrame(col_consecutive,columns=col)
        Final_df =pd.concat([Final_df, col_df], axis=0)

        sleep(0.01)

    

    #*************************Fourth Part*******************************8

    
    if agree_row:
        #Select only those DMR with Ratio<= Cutoff
        with st.spinner('Filtering DMRS'):
            Final_df =Final_df[Final_df['Ratio']<=r_max]

    
    with st.spinner("Getting Your File Ready"):
        # Download button for dmr file
        Final_df['Coordinate'] = 'chr'+Final_df['Coordinate']
        Final_df=Final_df.fillna("")
        dmr_btn = download_button_zip(Final_df,'DMR.csv',"Download DMRs File (ZIP)")
        st.markdown(dmr_btn, unsafe_allow_html=True)

    status.success('All DMRs Identified')

    


if submit and (data_files is None):
    status.error("Either one or both of Data File and Average File is not uploaded.")
