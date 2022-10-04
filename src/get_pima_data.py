import pandas as pd
import numpy as np
#
# 	THESE FUNCTIONS WERE INITIALLY DEFINED HERE
#	https://gist.githubusercontent.com/richardtjornhammar/1b9f5742391b1bcf30f4821a00f30b6a/raw/9542e711ef98fb53da042bd90f4968f1fe7aa056/pima.py
#
from bs4 import BeautifulSoup

def prune_whitespaces ( s ):
    return ( ' '.join([s_ for s_ in s.split(' ') if len(s_)>0]) )

def get_sample_data ( content , subject_id = [ 'Pima Indian (' , ')' ] , label = 'description' ,
                             references = ['sample-ref','platform-ref'],
                             structure_on = {'Subject:':'Subject' , 'Array type:':'Array' , 'Keywords =':'Types'} ) :
    #
    samples = [c.get('ref') for c in content.find_all(references[0])]
    platforms = [c.get('ref') for c in content.find_all(references[1])]
    sample_info = [ str(c).split('>')[1].split('<')[0].split('\n') for c in content.find_all(label) if 'Keywords' in str(c) ]
    sample_info = [ [ s for s in sample if len(s.replace(' ',''))>0 ] for sample in sample_info if 'Subject' in ''.join(sample) ]
    structured_samples = []
    #
    def prune_whitespaces ( s ):
        return ( ' '.join([s_ for s_ in s.split(' ') if len(s_)>0]) )
    #
    common_types = None
    for sample , id, platform in zip ( sample_info,samples,platforms ) :
        sample_ledger = { v:[] for v in structure_on.values() }
        kv = structure_on.keys()
        for s in sample :
            for k in kv :
                if k in s :
                    sample_ledger[ structure_on[k] ].append( prune_whitespaces(s.split(k)[-1].split('.')[0]) )
        if 'list' in str( type(subject_id)) :
            sample_ledger['Name'] = sample_ledger['Subject'][0].split(subject_id[0])[-1].split(subject_id[-1])[0]
            [ sample_ledger['Types'] .append( prune_whitespaces(t) ) for t in sample_ledger['Subject'][0].split(subject_id[0])[0].split(' ') ]
        sample_ledger['Platform'] = platform
        structured_samples.append ( sample_ledger )
        if common_types is None :
            common_types = set(sample_ledger['Types'])
        else :
            common_types = common_types & set(sample_ledger['Types'])
    sample_dictionary = { id:s for s,id in zip(structured_samples,samples) }
    for sample in sample_dictionary.items() :
        T = [ t for t in sample[1]['Types'] if not t in common_types ]
        sample_dictionary[ sample[0] ] ['Types'] = '-'.join(T)
        for t in range(len(T)):
            sample_dictionary[ sample[0] ] ['Type'+str(t)] = T[t]
    df_ = pd.DataFrame(sample_dictionary)
    df_ .loc[ 'Array' ] = [ ''.join(v) for v in df_.loc['Array'].values ]
    #
    return ( df_ , common_types )


def get_transcripts ( df , work_dir , data_file_suffix='-tbl-1.txt' , platform_label='Platform' ,skip=[] ) :
    suf_ = data_file_suffix
    pos_ = 9
    extra_sample_information = {}
    #  present (P), absent (A), marginal (M), or no call (NC)
    #  WE ALSO RETURN THE AMOUNT OF EACH TO THE LEDGER
    common_index = None
    all_data = []
    for c in df.columns :
        analytes = pd.read_csv(work_dir+c+suf_,header=None,index_col=0,sep='\t')
        which    = df.loc[platform_label,c]
        if which in set(skip):
            continue
        lookup   = pd.read_csv( work_dir+ which  +suf_,index_col=0,sep='\t').iloc[:,pos_]
        rd       = { i:prune_whitespaces(str(v)) for i,v in zip(lookup.index,lookup.values ) }
        extra_sample_information[c] = { 'Marginal': np.sum( analytes.iloc[:,1]=='M') ,
         'Present': np.sum( analytes.iloc[:,1]=='P') ,
         'Absent': np.sum( analytes.iloc[:,1]=='A') ,
         'NoCall': np.sum( analytes.iloc[:,1]=='NC') }
        adf = analytes .rename( index=rd )
        bUse = [ not 'nan' in str(v).lower() for v in adf.index.values]
        adf = adf.iloc[bUse,[0]].rename(columns={1:c})
        all_data.append( adf )
        if common_index is None :
            common_index = set ( adf.index.values )
        else :
            common_index = common_index & set(adf.index.values)
    return ( pd.concat(all_data,1) , pd.DataFrame(extra_sample_information) )


if __name__ == '__main__' :
    #
    # WE DOWNLOADED : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2508
    # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2508/miniml/GSE2508_family.xml.tgz
    import os
    os.system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2508/miniml/GSE2508_family.xml.tgz')
    os.system('gunzip GSE2508_family.xml.tgz')
    os.system('tar xvf GSE2508_family.xml.tar')
    #
    work_dir  = './'
    filename  = 'GSE2508_family.xml'
    full_path = work_dir + filename
    with open ( full_path , 'r' ) as input :
        content = BeautifulSoup ( input , features="html.parser")
    #
    # AND NOW WE PRUNE THE DATA TO GET A JOURNAL AND ANALYTE FRAME
    df_ , common_set = get_sample_data( content )
    print ( common_set )
    #
    ccat = {'Array':'C(Array)','Types':'C(Types)','Type0':'C(Type0)','Type1':'C(Type1)','Platform':'C(Platform)'}
    journal_df  = df_ .loc[ ccat.keys() ].rename(index=ccat)
    print ( journal_df )
    #
    skip_ = ['GPL91','GPL92','GPL93','GPL94','GPL95']
    analyte_df,einf_df = get_transcripts( df_ , work_dir , skip=skip_ )
    journal_df = pd.concat( [ journal_df.loc[:,einf_df.columns],einf_df] ,0 )
    formula = 'f~'+'+'.join(journal_df.index.values)
    #
    print ( formula )
    #
    from impetuous.quantification import multifactor_evaluation, multivariate_aligned_pca
    analyte_df.to_csv( '../data/pima_analytes_df.tsv' , sep='\t' )
    journal_df.to_csv( '../data/pima_journal_df.tsv'  , sep='\t' )
    multifactor_results = multifactor_evaluation ( analyte_df , journal_df , formula )

    journal_df.loc['Sample ID'] = journal_df.columns.values

    results = multivariate_aligned_pca ( analyte_df , journal_df , sample_label='Sample ID' ,
					 align_to = 'C(Types)' , add_labels = ['C(Platform)','C(Array)'] ,
					 n_components=15 )
    check_level = 0.95
    sig_analytes  = results[0].iloc[ [ v>check_level for v in results[0].loc[:,'Corr,r'].values ] , : ]
    sig_analytes2 = multifactor_results.iloc[[ v >check_level for v in multifactor_results.loc[:,'obese,r'].values ],:]
    print ( 'OBESITY TRANSCRIPTS:', set( sig_analytes.index.values ) & set( sig_analytes2.index.values ) )
    #
    if True :
        os.system("rm *.txt")
        os.system("rm *.tar")
        os.system("rm *.zip")
        os.system("rm *.xml")
