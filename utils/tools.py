import numpy as np 
import pandas as pd
import sys
import gzip
import urllib.request
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from matplotlib import pyplot as plt
import io
import os


###########################################################################################
###########################################################################################

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

###########################################################################################
###########################################################################################

def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)

###########################################################################################
###########################################################################################

def download_file(url):
    
    """
    download .gz file and decompress locally                    
    -> function used to acces the feature profile from the main 'Train_Data' file
    
    """
    out_file = 'ile'

       # Download archive
    try:
          # Read the file inside the .gz archive located at url
        with urllib.request.urlopen(url) as response:
            with gzip.GzipFile(fileobj=response) as uncompressed:
                file_content = uncompressed.read()

          # write to file in binary mode 'wb'
        with open(out_file, 'wb') as f:
            f.write(file_content)
            return 0

    except Exception as e:
        print(e)
        return 1
    
###########################################################################################
###########################################################################################

def load_bed (file):

        """
        load a .bed file in the right format

        """
        header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd',
              'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
        bed_file = pd.read_csv(file, sep = '\t', comment = 't', header = None)
        bed_file.columns = header[:len(bed_file.columns)]
        return bed_file 
    
###########################################################################################
###########################################################################################

def load_narrowPeak (file):
    
        """
        load a .narrowPeak (BED4+6) file in the right format

        """
        header = ['chrom', 'peakStart', 'peakEnd', 'name', 'peakScore',
                         'strand', 'PeakSignalValue', 'p-value', 'qvalue', 'peak',
                         'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
                         'blockSizes', 'blockStarts']
        narrow_file = pd.read_csv(file, sep = '\t', comment = 't', header = None)
        narrow_file.columns = header[:len(narrow_file.columns)]
        return narrow_file 
    
    
###########################################################################################
###########################################################################################

def load_out_bed (file, index = False):
    
        """
        load a .bed.out (output of DeepSEA) file in the right format
        
        """
        out_file = pd.read_csv(file , index_col=index)
        out_file.rename(columns = {'Unnamed: 0' : "seq-count", 'chr' : 'chrom', 'start' : 'chromStart', 
                                   'end' : 'chromEnd'}, inplace = True)
        out_file.drop(columns = ['seq-count'], inplace = True)
        return out_file 
    
###########################################################################################
###########################################################################################

def isinPeak (seq_start_idx, seq_end_idx, peak_start_idx, peak_end_idx) -> bool:

        """
        Checks if more than half of a smaller DNA sequence (typically 200-bp) is 
        contained in a larger chromatin peak region of the chromosome 
        
        -> Returns True if this is the case, False otherwise 

        """ 
        
        if seq_start_idx > seq_end_idx or peak_start_idx > peak_end_idx:
            raise ValueError('start index should be smaller than end index')
            
        seq_start_idx += 400
        seq_end_idx -= 400    
        
        if seq_end_idx <= peak_end_idx and seq_start_idx >= peak_start_idx:
            return True 
        if seq_start_idx < peak_start_idx:
            if seq_end_idx < peak_start_idx:
                return False
            elif seq_end_idx >= peak_start_idx:
                if seq_end_idx - peak_start_idx >= np.rint(0.5*(seq_end_idx - seq_start_idx)):
                    return True 
                else:
                    return False
        if seq_end_idx > peak_end_idx:
            if seq_start_idx > peak_end_idx:
                return False 
            elif seq_start_idx <= peak_end_idx:
                if peak_end_idx - seq_start_idx >= np.rint(0.5*(seq_end_idx - seq_start_idx)):
                    return True
                else:
                    return False
                
###########################################################################################
###########################################################################################
                
def to_DeepSEAbed (input, file_name):
    
    """
    generate an input .bed file accepted by DeepSEA 
    
    input is a dataframe containing all of the sequences of 1000-bp
    
    """
    #df = input[input[input.columns[pd.Series(input.columns).str.startswith('chrom')]]==1]
    
    #out = input.drop(columns = ['TFbind'], inplace = True)

    input.to_csv(file_name, sep = '\t', index=False, header=False)
    
    
###########################################################################################
###########################################################################################
    
    
def generate_chro_profile_from_file(input, chr):

    """
    generates the full chromatin profile results  for a given pair of chromosome/feature
    
    requires full path to a .narrowPeak file where the profile is stored
    
    """
    profile = load_narrowPeak(input)
    profile = profile.loc[profile['chrom'] == chr]  # extract information for chromosome 8
    profile = profile.sort_values(by = 'peakStart')
    profile['peakLength'] = profile.loc[:, 'peakEnd'] - profile.loc[:,'peakStart']
    profile = profile.drop_duplicates(['peakStart']) # remove duplicates 
    
    return profile

###########################################################################################
###########################################################################################

def cut_profile(input, frac):
    
    """
    returns a percentage of the profile in terms of number of called peaks 
    
    """
    cutoff = np.int(frac*input.shape[0])
    return input[:cutoff]

###########################################################################################
###########################################################################################

def generate_TF_pos(profile_start, range, chr, nb_samples):
                    #, full = False ):

    """
    generates a set of random 200-bp sequences within a specified range of the chromosome profile 
    
    different ranges are possible 
    
        -> full : the range of generated sequences spans the entire chromatin profile (all peaks are considered)
        -> if Not full : a specified range  of base pairs starting from start of profile
    
    """
    #profile_start = profile.iloc[0]['peakStart']
    #profile_end = profile.iloc[profile.shape[0]-1]['peakEnd']
    #window_size = 200   # we focus on 200-bp sequences 
    #if full == False:
    
    rg = profile_start + range  # range of analysis from start position    
    #elif full == True:
        #rg = profile_end
                
    start_id  = np.random.randint(profile_start, profile_start + range, size=nb_samples)
    end_id = np.add(start_id, 1000)
    name = np.full(np.size(start_id), chr)
    TF_pos = pd.DataFrame({'chrom':name, 'chromStart':start_id, 'chromEnd':end_id} )
    
    return TF_pos

###########################################################################################
###########################################################################################

def generate_mini_test_set(reads, profile):

    reads['TFbind'] = np.full(reads.shape[0], 0)

    for j in profile.index:
    
        peak_start = profile.at[j, 'peakStart']
        peak_end = profile.at[j, 'peakEnd']
        
        for i in reads.index:
            seq_start = reads.at[i, 'chromStart']
            seq_end = reads.at[i, 'chromEnd']

            if isinPeak (seq_start, seq_end, peak_start, peak_end): #and profile.iloc[j]['peakScore']: # be more stringent on peakscore to reduce background ?
                reads.at[i, 'TFbind'] = 1
    return reads
                

###########################################################################################
###########################################################################################


def extract_features(all=False, dnase=False, tf=False, histone_mark=False):
    features = load_out_bed('..\\tests\\outputs\\test.bed.out').drop(columns = ['chrom', 'chromStart', 'chromEnd'])
    feature_names = list(features.columns)

    if all:
        return feature_names
    elif tf:
        return feature_names[125:815]  # 690 TF features
    elif dnase:
        return feature_names[1:125]   # 125 dnase features ( we don't keep the first because no data for it it ENCODE file)
    elif histone_mark:
        return feature_names[815:919]  # 104 histone features 

###########################################################################################
###########################################################################################    
    
def extract_ref(extracted_features, feature_index, dnase=False, tf=False, hk=False):
    
    df = pd.read_csv('..\\data\\pos.bed', sep = '\t')
    
    propos_dnase = df[0:124]
    propos_tf = df[124:814]
    propos_hk = df[814:918]

    idx = extracted_features.index(feature_index)   
    my_feature = extracted_features[idx].replace('|', '_')
    
    if dnase:
        prop = propos_dnase.iloc[idx]['Positive Proportion']
        return [idx, my_feature, prop]
    elif tf:
        prop = propos_tf.iloc[idx]['Positive Proportion']
        return [idx, my_feature, prop]
    elif hk:
        prop = propos_hk.iloc[idx]['Positive Proportion']
        return [idx, my_feature, prop]

     
###########################################################################################
########################################################################################### 

def predict_from_DeepSEA(feature_refs, test_input, test_output, thresh, color):
    
    target = test_input['TFbind'].to_numpy()
    ones = np.count_nonzero(target == 1)
    yhat   = np.empty(len(target))

    # exceptions 
    if 'SK-N-SH_RA' in feature_refs[1]:
        yhat = test_output[rreplace(feature_refs[1], '_', '|', 2)].to_numpy()
        
    elif 'BE2_C' in feature_refs[1]:
        yhat = test_output[rreplace(feature_refs[1], '_', '|', 2)].to_numpy()
        
    elif 'HSMM_emb' in feature_refs[1]:
        yhat = test_output[rreplace(feature_refs[1], '_', '|', 2)].to_numpy()
    
    elif 'Adult_CD4_Th0' in feature_refs[1]:
        yhat = test_output[rreplace(feature_refs[1], '_', '|', 2)].to_numpy()
        
    elif 'CD34+_Mobilized' in feature_refs[1]:
        yhat = test_output[rreplace(feature_refs[1], '_', '|', 2)].to_numpy()
        
    elif 'Monocytes-CD14+_RO01746' in feature_refs[1]:
        yhat = test_output[rreplace(feature_refs[1], '_', '|', 2)].to_numpy()
        
    elif 'Monocytes-CD14+_RO01746' in feature_refs[1]:
        yhat = test_output[rreplace(feature_refs[1], '_', '|', 2)].to_numpy()
    
    else:
        yhat = test_output[feature_refs[1].replace('_', '|', 2)].to_numpy()
    
    c_train = feature_refs[2]
    
    yhat_norm = (1/(1+np.exp(-(np.log(yhat/(1-yhat))+np.log(5/100/(1-5/100))-np.log(c_train/(1-c_train ))))))
    
    print('feature {:d}\ntag : {}\nchromosome : {} [{:d} - {:d}] \n positives = {:d}'
          
                  .format(feature_refs[0]+1, feature_refs[1], test_input.loc[0]['chrom'], test_input['chromStart'].min(), 
                          test_input['chromEnd'].max(), ones))
    print(' \n')
    
    # calculate roc curves
    fpr, tpr, thresholds = roc_curve(target, yhat_norm)

    score = 0

    if ones >= thresh:
    
        # plot the roc curve for the model
        plt.plot(fpr, tpr, color=color, linewidth=0.8)
        score = roc_auc_score(target, yhat_norm) 

    return target, yhat, yhat_norm, ones, score 
    

###########################################################################################
###########################################################################################     
    
def unwrap_profiles(all=False, dnase=False, tf=False, hm=False):

    """
    unwrap the zipped url links from .xls files
    -> all : unwrap all feature links
    -> or unwrap specific feature type if specified 
    
    """
    df = pd.read_excel('..\\data\\TrainData.xlsx')
    
    if all:
        return df
    elif dnase:
        df = df[0:124]
    elif tf:
        df = df[124:814]
    elif hm:
        df = df[814:918]
        
    return df

###########################################################################################
###########################################################################################

def download_file(url, id, dnase=False, tf=False, hk=False):
    
    out_path = None
    extension = '.bed'
    if dnase:
        out_path = '..\data\Dnase_chromatin_profiles\\'
    
    elif tf:
        out_path = '..\data\TF_chromatin_profiles\\'
        
    elif hk:
        out_path = '..\HK_chromatin_profiles\\'
        
    out_file = out_path + id + extension
    
       # Download archive
    try:
          # Read the file inside the .gz archive located at url
        with urllib.request.urlopen(url) as response:
            with gzip.GzipFile(fileobj=response) as uncompressed:
                file_content = uncompressed.read()

          # write to file in binary mode 'wb'
        with open(out_file, 'wb') as f:
            f.write(file_content)
            return 0

    except Exception as e:
        print(e)
        return 1
    
###########################################################################################
###########################################################################################  

def generate_input(start_idx, range_, chromosome, n_test, file_name):

    input_seqs = generate_TF_pos(start_idx, range_, chromosome, n_test)
    to_DeepSEAbed(input_seqs, '..\\data\\BED_input\\' + file_name )


def boxplot(data, title='TF', marker='r.'):

    """ 
    plots a single boxplot for a single model type
    
    INPUT : roc auc scores 
    """


    boxdict = dict(linestyle='-', linewidth=2, color='black')
    whiskerdict = dict(linestyle='-', linewidth=2, color='black')
    mediandict = dict(linestyle='--', linewidth=1.5, color='red')

    fig1, ax1 = plt.subplots(1,1,figsize=(10,7))

    ax1.set_title(title)

    bplot = ax1.boxplot(data, patch_artist=False, showfliers=True, showcaps=False, boxprops=boxdict, whiskerprops=whiskerdict, medianprops=mediandict)
    plt.xticks([1], [''])

    y = data
    x = np.random.normal(1, 0.04, size=len(y))
    ax1.plot(x, y, marker, alpha=0.7)

    ax1.yaxis.grid(True)
    ax1.set_xlabel('')
    ax1.set_ylabel('AUC score (%)')