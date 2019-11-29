import pandas, os
from collections import Counter

# Download urls for MarRef and MarDB metadata (.tsv)
urls = {'MarRef': 'https://s1.sfb.uit.no/public/mar/MarRef/Metadatabase/Current.tsv', 
        'MarDB': 'https://s1.sfb.uit.no/public/mar/MarDB/Metadatabase/Current.tsv'}

# Comprehensive mapping of metadata column headers to inputdata (pandas) coulmn headers to avoid invalid couln names
strheadermap = {
        'mmp_ID': 'mmp_ID',
        'analysis_project_type': 'analysis_project_type',
        'full_scientific_name': 'full_scientific_name',
        'rrnas': 'rrnas'
        }

floatheadermap = {
        'depth': 'depth',
        'env_salinity': 'env_salinity',
        'env_temp': 'env_temp',
        'num_replicons': 'num_replicons',
        'genes': 'genes',
        'cds': 'cds',
        'pseudo_genes': 'pseudo_genes',
        #'rrnas': 'rrnas', (This is str before reformatting below.)
        'tRNAs': 'trnas',
        'frameshifted_genes': 'frameshifted_genes',
        'optimal_temperature': 'optimal_temperature',
        'plasmids': 'plasmids',
        'Assembly_length': 'genome_length',
        'GC': 'gc_content',
        'sequencing_depth': 'sequencing_depth',
        'contigs': 'contigs',
        'Completeness': 'estimated_completeness',
        'Contamination': 'estimated_contamination',
        'Strain_heterogeneity': 'strain_heterogeneity',
        'QS': 'qs'

        }

complete_df = pandas.DataFrame()

# Iterate over DB's and urls to s1.sfb.uit.no
for db, url in urls.items():
    # Set tmp df, to store curent data iteration in
    tmp = pandas.DataFrame()
    # Download and parse mmp metadta into pandas df
    metadata = pandas.read_csv(url, sep='\t')
    # Create tmpdf columns with str, ensuring correct headernames
    for new, old in strheadermap.items():
        tmp[new] = metadata[old].astype(str)

    # Create tmpdf using headermap, ensuring correct headernames 
    for new, old in floatheadermap.items():
        tmp[new] = pandas.to_numeric(metadata[old], errors='coerce')
        
    # Make DB and db columns used for Tooltip (CamelCase) and url (lowercase) respectively
    tmp['DB'] = str(db)
    tmp['db'] = tmp['DB'].str.lower()

    # Fix rRNA triplet bars (Metadata has one rrnas column, visualization template should have 3, rRNA5S, rRNA16S and rRNA23S)
    rrnas = tmp['rrnas'].str.split(pat=',', n=-1,  expand=True)
    tmp['rRNA5S'] = pandas.to_numeric(rrnas[0], errors='coerce')
    tmp['rRNA16S'] = pandas.to_numeric(rrnas[1], errors='coerce')
    tmp['rRNA23S'] = pandas.to_numeric(rrnas[2], errors='coerce')

    # Concatonate this df with that df
    complete_df = pandas.concat([complete_df, tmp], ignore_index=True)

# Calculate Draft quality based on data from dataframe
quality = []
for index, row in complete_df.iterrows():
    if row['DB'] == 'MarRef':
        quality.append('Finished')
    elif row['Completeness'] == 'missing' or row['Contamination'] == 'missing':
        quality.append('NA')
    elif float(row['Completeness']) > 90 and float(row['Contamination']) < 5:
        quality.append('High Quality Draft')
    elif float(row['Completeness']) >= 50 and float(row['Contamination']) < 10:
        quality.append('Medium Quality Draft')
    elif float(row['Completeness']) < 50 and float(row['Contamination']) < 10:
        quality.append('Low Quality Draft')
    elif float(row['Contamination']) > 10:
        quality.append('Very Low Quality Draft')
    else:
        quality.append('NA')
complete_df['quality'] = quality

# Set indivival colors for discrete databases
colorarray = []
for DB, apt in zip(complete_df['DB'], complete_df['analysis_project_type']):
        if DB == "MarRef" and apt == "missing":
                colorarray.append("grey")
        elif DB == "MarDB" and apt == "missing":
                colorarray.append("grey")
        elif DB == "MarRef":
                colorarray.append("blue")
        elif DB == "MarDB":
                colorarray.append("green")
        else:
                colorarray.append("grey")
complete_df['color'] = colorarray

# Tune alphas for analysis_project_type
alpha = []
for analysis_type in complete_df['analysis_project_type']:
        if analysis_type == "Whole genome sequencing (WGS)":
                alpha.append(1)
        elif analysis_type == "Metagenome assembled genome (MAG)":
                alpha.append(0.6)
        elif analysis_type == "Single amplified genome (SAG)":
                alpha.append(0.4)
        else:
                alpha.append(1)
complete_df['alpha'] = alpha

# Make labels for legend
label= []
for DB, apt in zip(complete_df['DB'], complete_df['analysis_project_type']):
        if DB == "MarRef" and apt == "missing":
                label.append("MarRef (Unknown type)")
        elif DB == "MarDB" and apt == "missing":
                label.append("MarDB (Unknown type)")
        elif DB == "MarRef" and apt == "Whole genome sequencing (WGS)":
                label.append("MarRef (Whole genome sequencing (WGS))")
        elif DB == "MarDB" and apt == "Whole genome sequencing (WGS)":
                label.append("MarDB (Whole genome sequencing (WGS))")
        elif DB == "MarRef" and apt == "Metagenome assembled genome (MAG)":
                label.append("MarRef (Metagenome assembled genome (MAG))")
        elif DB == "MarDB" and apt == "Metagenome assembled genome (MAG)":
                label.append("MarDB (Metagenome assembled genome (MAG))")
        elif DB == "MarRef" and apt == "Single amplified genome (SAG)":
                label.append("MarRef (Single amplified genome (SAG))")
        elif DB == "MarDB" and apt == "Single amplified genome (SAG)":
                label.append("MarDB (Single amplified genome (SAG))")
        else:
                label.append("Unknown")
complete_df['label'] = label

complete_df.to_csv(os.path.join(os.getcwd(), "mmp_interactive/inputdata/input.tsv"), sep="\t")

