from __future__ import division
import numpy as np
import pandas as pd
from collections import OrderedDict
from scipy import sparse
import timeit
from mNMF01_Tools import msconvert_reader
from mNMF02_Insilicodigestion_Yeast import insilico_detect, param_autoread
from mNMF03_RefmatConstruction import refMS1_construction
from mNMF04_RefNoiseConstruction import noise_subspace
from mNMF05_ExpmatConstruction import expmat_construction
from mNMF06_RTPrediction import rt_prediction
from mNMF07_WeightEstimation import weight_estimation
from mNMF08_HConstruction import h_construction
from mNMF09_NMFrun import kl_nmf
from mNMF10_NMFReport import nmf_identification, nmf_savemat

# define global param
defParams_ = True
if defParams_:
    # mz and time parameters
    mm, tt = 100, 200 # default mz bin = 1/100 per m/z and tt bin = 1/200 per minute
    # smoothing parameters
    window_slideminute = 12/60 # 0.2 minute
    window = int(window_slideminute*tt) 
    shift = int(window/2) # default shift = 50% of a window = 0.1 minute            
    # amino acid dictionary
    aa_monomass = {
        'A': 71.03711,'C': 103.00919,'D': 115.02694,
        'E': 129.04259,'F': 147.06841,'G': 57.02146, 
        'H': 137.05891,'I': 113.08406,'K': 128.09496,
        'L': 113.08406,'M': 131.04049,'N': 114.04293,
        'O': 255.15829,'P': 97.05276,'Q': 128.05858,
        'R': 156.10111,'S': 87.03203,'T': 101.04768,
        'U': 150.95364,'V': 99.06841,'W': 186.07931,
        'Y': 163.06333,   
        }   
    aa_composition = {
        'A':   ({'H': 5, 'C': 3, 'O': 1, 'N': 1}),
        'C':   ({'H': 5, 'C': 3, 'O': 1, 'N': 1, 'S': 1}),
        'D':   ({'H': 5, 'C': 4, 'O': 3, 'N': 1}),
        'E':   ({'H': 7, 'C': 5, 'O': 3, 'N': 1}),
        'F':   ({'H': 9, 'C': 9, 'O': 1, 'N': 1}),
        'G':   ({'H': 3, 'C': 2, 'O': 1, 'N': 1}),
        'H':   ({'H': 7, 'C': 6, 'O': 1, 'N': 3}),
        'I':   ({'H': 11, 'C': 6, 'O': 1, 'N': 1}),
        'K':   ({'H': 12, 'C': 6, 'O': 1, 'N': 2}),
        'L':   ({'H': 11, 'C': 6, 'O': 1, 'N': 1}),
        'M':   ({'H': 9, 'C': 5, 'O': 1, 'N': 1, 'S': 1}),
        'N':   ({'H': 6, 'C': 4, 'O': 2, 'N': 2}),
        'O':   ({'H': 19, 'C': 12, 'O': 2, 'N': 3}),
        'P':   ({'H': 7, 'C': 5, 'O': 1, 'N': 1}),
        'Q':   ({'H': 8, 'C': 5, 'O': 2, 'N': 2}),
        'R':   ({'H': 12, 'C': 6, 'O': 1, 'N': 4}),
        'S':   ({'H': 5, 'C': 3, 'O': 2, 'N': 1}),
        'T':   ({'H': 7, 'C': 4, 'O': 2, 'N': 1}),
        'U':   ({'H': 5, 'C': 3, 'O': 1, 'N': 1, 'Se' : 1}),
        'V':   ({'H': 9, 'C': 5, 'O': 1, 'N': 1}),
        'W':   ({'H': 10, 'C': 11, 'O': 1, 'N': 2}),
        'Y':   ({'H': 9, 'C': 9, 'O': 2, 'N': 1}),
        'H-':  ({'H': 1}),
        '-OH': ({'O': 1, 'H': 1}),
        }   
    # modification dictionary
    aa_fixmodmass = {'C':160.03065} # Cys_CAM
    aa_fixmodcompos = {'C':({'H': 7, 'C': 5, 'O': 2, 'N': 2, 'S': 1})}
    aa_varmodmass = OrderedDict([('Oxidation (M)'  ,(+15.99491,'all',{'O': +1},'M',1))]) 
    if len(aa_fixmodmass) > 0:
        for f_k, f_v in aa_fixmodmass.items():
            aa_monomass[f_k] = f_v
        for f_k, f_v in aa_fixmodcompos.items():
            aa_composition[f_k] = f_v
    # other dictionary
    other_mass = {'Proton':1.00728, 'Oxygen':15.99491, 'H20':18.01056,'Phospho':79.9663}
    el_abundance = {
        'H':   OrderedDict([('H1', 0.9999),('H2', 0.0001)]),
        'C':   OrderedDict([('C12', 0.9889), ('C13', 0.0111)]),
        'O':   OrderedDict([('O16', 0.9976), ('O17', 0.0004), ('O18' , 0.0020)]),
        'N':   OrderedDict([('N14', 0.9964), ('N15', 0.0036)]),
        'P':   OrderedDict([('P31', 1)]),
        'S':   OrderedDict([('S32', 0.9500), ('S33', 0.0076), ('S34', 0.0422), ('S35', 0), ('S36', 0.0002)]),
        }
    dictionary_list = [aa_monomass, aa_composition, aa_fixmodmass, aa_varmodmass, \
                        other_mass, el_abundance]
    eps = 1e-15

# define data
defData_ = True
if defData_:
    msconvert_file = 'Example_yeast500ng_rep1.mzML' #txt or pkl
    maxquant_file = 'Example_yeast500ng_MaxQuantnotMBR_toNMFdictionary.xlsx'
    insilico_file = maxquant_file[:-5]+'_insilico.xlsx'
    output_file = 'Output_' + msconvert_file[:-4] #filename

### 001 msconvert read tool
print(' report time: ', timeit.default_timer())
print('1 msconvert processing')
if msconvert_file[-3:] == 'pkl':
    Vdata_file = msconvert_file
else:
    Vdata_file = msconvert_reader(msconvert_file)
print(' report time: ', timeit.default_timer())
#### 2 in silico digestion
print('2 in silico digestion')
# def in silico parameters
misclvge_from, misclvge_to = 0,2
len_from, len_to = 7,50
charge_from, charge_to = 2,4
insilicoparam_list = [misclvge_from, misclvge_to, len_from, len_to, charge_from, charge_to]
# def params as default
WithIsotopicPeak_ = True
if WithIsotopicPeak_:
    isopeak_number = 6
else:
    isopeak_number = 1
# run
try:
    insilico_df = pd.read_excel(insilico_file)
except:
    insilico_df = insilico_detect(maxquant_file, dictionary_list, insilicoparam_list, isopeak_number)
    insilico_df.to_excel(insilico_file)
mz_header, iso_header, charge_list, globalparam_list = param_autoread(insilico_df, isopeak_number, mm, tt, window, shift)

# print(mz_header,charge_list,iso_header)
print(' report time: ', timeit.default_timer())
### 3 W_mat construction
print('3 W construction')
# run
W_mat, finalMS1_df, prot_peptcount, pept_ioncount \
= refMS1_construction(insilico_df, mz_header, iso_header, charge_list, isopeak_number, globalparam_list, eps)
del(insilico_df)
print(' report W: ', W_mat.shape)
print(' report time: ', timeit.default_timer())
# finalMS1_df.to_excel(maxquant_file[:-5]+'_finalinsilico_df.xlsx')

### 4 noise subspace
print('4 noise subspacing')
# def
WithNoiseSubspace_ = True
if WithNoiseSubspace_:
    noise_extractfrom, noise_extractto = 100, 110
    noise_number = 2
    W_mat, noise_mean \
    = noise_subspace(Vdata_file, W_mat, noise_number, globalparam_list, charge_list, noise_extractfrom, noise_extractto)
print(noise_mean)
sparse.save_npz(output_file + '_Wmat_init.npz', W_mat)
print(' report time: ', timeit.default_timer())

### 5 V_mat Construction
print('5 V Construction')
V_mat, globalparam_list, Vmat_mean\
= expmat_construction(Vdata_file, globalparam_list, charge_list)
sparse.save_npz(output_file + '_Vmat_init.npz', V_mat)
print(' report time: ', timeit.default_timer())

### 6 RT prediction
print('6 RT prediction')
# def
WithRTPrediction_ = True
if WithRTPrediction_:
    gaussian_peakwidth = 3 # minutes
else:
    gaussian_peakwidth = 0
RT_equation = None

if WithRTPrediction_:
    initRT_tuple, initRT_width\
    = rt_prediction(maxquant_file, finalMS1_df, RT_equation)
else:
    initRT_tuple, initRT_width = None, None  
print(' report time: ', timeit.default_timer())

### 7 Penalty and weight estimation from sparse measure
print('7 Weight estimation')
# WithL1_ = True
# if WithL1_:
#     weight_est = weight_estimation(V_mat, prot_peptcount, globalparam_list)
# else:
#     weight_est = 0
weight_est = 1
print(' report weight:', weight_est)
print(' report time: ', timeit.default_timer())

### 8 H_mat Construction  
print('8 H Construction')
# def
WithCosineScreen_ = False
if WithCosineScreen_:
    cosine_cutoff = 0.5
else:
    cosine_cutoff = 0
# run
H_mat, initRT, cos_mat = h_construction(WithRTPrediction_, WithCosineScreen_, initRT_tuple, initRT_width, \
                            V_mat, W_mat, noise_number, globalparam_list, isopeak_number, gaussian_peakwidth, cosine_cutoff, output_file)
sparse.save_npz(output_file + '_Hmat_init.npz', H_mat)
print(' report time: ', timeit.default_timer())

### 9 NMF 
print('9 NMF run and 10 report')
# def
WithStandardize_ = True
max_iter = 50 # number of max iteration
peak_report = 1 # number of best peaks for report
kl_nmf(V_mat, W_mat, H_mat, cos_mat, prot_peptcount,\
    weight_est, max_iter, eps, noise_number, noise_mean, Vmat_mean, globalparam_list, isopeak_number,\
    initRT, output_file, finalMS1_df, peak_report, WithStandardize_)
print(' report time: ', timeit.default_timer())


