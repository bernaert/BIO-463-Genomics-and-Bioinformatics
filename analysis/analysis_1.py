import numpy as np 
import pandas as pd
import sys
import gzip
import urllib.request
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from matplotlib import pyplot
sys.path.append('..')
from utils import tools as tl

##################################
###### Transciption Factors ######
##################################

h = ['chrom', 'chromStart', 'chromEnd']
TF_seqs_ = pd.read_csv('..\\data\\BED_test_data\\chrom9_5000.bed', sep = '\t', comment = 't', header = None)
TF_seqs_.columns = h[:len(TF_seqs_.columns)]
TF_seqs_sorted = TF_seqs_.sort_values(by=['chromStart'], ascending=True).reset_index()
start_idx = TF_seqs_sorted['chromStart'].min()
end_idx = TF_seqs_sorted['chromEnd'].max()
TFs_output = tl.load_out_bed('..\\data\\BED_test_outputs\\chrom9_5000.bed.out')
all_tf_features = tl.extract_features(tf=True )

pyplot.plot([0,1], [0,1], linestyle='--',  color='k', linewidth=0.5) # plot the 'random model' lines

all_scores = []

for i in range (len(all_tf_features)):

        tf_ref = tl.extract_ref(all_tf_features, all_tf_features[i], tf = True)
        id = tf_ref[1]
        profile = tl.generate_chro_profile_from_file('..\\data\\TF_chromatin_profiles\\' + id + '.bed', 'chr9')
        filter_profile = profile[(profile['peakStart'] >= start_idx) & (profile['peakEnd'] <= end_idx) ]
        test_input = tl.generate_mini_test_set(TF_seqs_sorted, filter_profile)
        target, yhat, yhat_norm, ones, score = tl.predict_from_DeepSEA(tf_ref, test_input, TFs_output, 20, color='red')

        all_scores.append(score)


all_scores = list(filter(lambda x: x!= 0, all_scores))
# axis labels
pyplot.xlabel('False Positive Rate')
pyplot.ylabel('True Positive Rate')
pyplot.title('')
# pyplot.legend()
# show the plot
pyplot.show()

tl.boxplot(all_scores, title='', marker='r.')


#############################################
####### DNase hypersensitive sites ##########
#############################################

h = ['chrom', 'chromStart', 'chromEnd']
DNASE_seqs_ = pd.read_csv('..\\data\\BED_test_data\\chrom9_5000.bed', sep = '\t', comment = 't', header = None)
DNASE_seqs_.columns = h[:len(DNASE_seqs_.columns)]
DNASE_seqs_sorted = DNASE_seqs_.sort_values(by=['chromStart'], ascending=True).reset_index()
start_idx = DNASE_seqs_sorted['chromStart'].min()
end_idx = DNASE_seqs_sorted['chromEnd'].max()
output = tl.load_out_bed('..\\data\\BED_test_outputs\\chrom9_5000.bed.out')
all_dnase_features = tl.extract_features(dnase=True )

pyplot.plot([0,1], [0,1], linestyle='--',  color='k', linewidth=0.5) # plot the 'random model' lines

all_scores_2 = []


for i in range (len(all_dnase_features)):

        dnase_ref = tl.extract_ref(all_dnase_features, all_dnase_features[i], dnase = True)
        id = dnase_ref[1]
        profile = tl.generate_chro_profile_from_file('..\\data\\DNASE_chromatin_profiles\\' + id + '.bed', 'chr9')
        filter_profile = profile[(profile['peakStart'] >= start_idx) & (profile['peakEnd'] <= end_idx) ]
        test_input = tl.generate_mini_test_set(DNASE_seqs_sorted, filter_profile)
        
        target, yhat, yhat_norm, ones_2, score = tl.predict_from_DeepSEA(dnase_ref, test_input, output, 20, color='green')


        all_scores_2.append(score)


all_scores_2 = list(filter(lambda x: x!= 0, all_scores_2))
# axis labels
pyplot.xlabel('False Positive Rate')
pyplot.ylabel('True Positive Rate')
pyplot.title(' ')
# pyplot.legend()
# show the plot
pyplot.show()

tl.boxplot(all_scores_2, title='  ', marker='g.')


# axis labels
pyplot.xlabel('False Positive Rate')
pyplot.ylabel('True Positive Rate')
pyplot.title('  ')
# pyplot.legend()
# show the plot
pyplot.show()



##############################
####### Histone Marks ########
##############################

h = ['chrom', 'chromStart', 'chromEnd']
HM_seqs_ = pd.read_csv('..\\data\\BED_test_data\\chrom9_5000.bed', sep = '\t', comment = 't', header = None)
HM_seqs_.columns = h[:len(HM_seqs_.columns)]
HM_seqs_sorted = HM_seqs_.sort_values(by=['chromStart'], ascending=True).reset_index()
start_idx = HM_seqs_sorted['chromStart'].min()
end_idx = HM_seqs_sorted['chromEnd'].max()
output = tl.load_out_bed('..\\data\\BED_test_outputs\\chrom9_5000.bed.out')
all_hm_features = tl.extract_features(histone_mark=True )

pyplot.plot([0,1], [0,1], linestyle='--',  color='k', linewidth=0.5) # plot the 'random model' lines

all_scores_3 = []

for i in range (len(all_hm_features)):

        hm_ref = tl.extract_ref(all_hm_features, all_hm_features[i], hk=True)
        id = hm_ref[1]
        profile = tl.generate_chro_profile_from_file('..\\data\\HM_chromatin_profiles\\' + id + '.bed', 'chr9')
        filter_profile = profile[(profile['peakStart'] >= start_idx) & (profile['peakEnd'] <= end_idx) ]
        test_input = tl.generate_mini_test_set(HM_seqs_sorted, filter_profile)
        
        target, yhat, yhat_norm, ones_3, score = tl.predict_from_DeepSEA(hm_ref, test_input, output, 20, color='blue')

        all_scores_3.append(score)


all_scores_3 = list(filter(lambda x: x!= 0, all_scores_3))
# axis labels
pyplot.xlabel('False Positive Rate')
pyplot.ylabel('True Positive Rate')
pyplot.title(' ')
# pyplot.legend()
# show the plot
pyplot.show()

tl.boxplot(all_scores_3, title=' ', marker='b.')

