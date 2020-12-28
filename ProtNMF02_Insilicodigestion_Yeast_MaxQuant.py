from __future__ import division
import numpy as np
import pandas as pd
from collections import Counter
import math
import itertools
import gc
from mNMF01_Tools import mz_scaling, time_scaling

### MAIN ###
def insilico_detect(mascot_file, dictionary_list, param_list, iso_maxnumber):
	gc.disable()
	print(' from:', str(mascot_file))
	maxquant_df = pd.read_excel(mascot_file)
	aa_monomass, aa_composition, aa_fixmodmass, aa_varmodmass, other_mass, el_abundance = dictionary_list
	misclvge_from, misclvge_to, len_from, len_to, charge_from, charge_to = param_list
	prot_list, seq_list, modseq_list, mod_list, \
	charge_list, mass_list, mz_list, \
	rt_list, iso_list \
	= [], [], [], [], [], [], [], [], []
	
	maxquant_df = maxquant_df.sort_values(by='Score', ascending=False)
	maxquant_df = maxquant_df.drop_duplicates(subset=['Leading proteins','Modified sequence','Charge'], keep='first')
	maxquantuniquekey_dict = {}
	maxquantuniquepept_dict = {}

	for index, row in maxquant_df.iterrows():			
		leadproteins = row['Leading proteins']
		leadproteins = leadproteins.split(';')
		sequence = row['Sequence']
		modseq = str(row['Modified sequence'])
		mod = str(row['Modifications'])

		charge = row['Charge']
		mass = row['Mass']
		mz = row['m/z']
		rt = row['Retention time']

		comp_ans = comp_calc(sequence, aa_composition)
		ox_number = row['Oxidation (M)']
		comp_ans['O'] += ox_number
		iso_ans = iso_calc(comp_ans, el_abundance, iso_maxnumber)

		for protein in leadproteins:
			# protein must be swissprot
			if protein[:3] == 'sp|':
				key = (protein, modseq, charge)
				if key not in maxquantuniquekey_dict:					
					prot_list.append(protein)
					seq_list.append(sequence)
					modseq_list.append(modseq)
					mod_list.append(mod)
					mass_list.append(mass)
					mz_list.append(mz)
					charge_list.append(charge)
					# comp_list.append(comp_ans)
					rt_list.append(rt)
					iso_list.append(iso_ans)
					maxquantuniquekey_dict[key] = 1

	M_df = pd.DataFrame({'prot': prot_list,
			'pept': seq_list,
			'modpept': modseq_list,
			'mod': mod_list,
			'[M]': mass_list,
			'mz' : mz_list,
			'charge': charge_list,
			'rtpeak': rt_list			
			})

	charge_from = np.amin(charge_list)
	charge_to = np.amax(charge_list)
	print('re charges:', charge_from, charge_to)
	M_df = M_df[['prot','pept','modpept','mod','[M]', 'mz', 'charge', 'rtpeak']]
	# M_df = mz_calc(M_df, aa_monomass, other_mass, ms1param_list, iso_maxnumber, charge_from, charge_to)
	# M_df = iso_calc(M_df, el_abundance, iso_maxnumber)
	rt_df = pd.pivot_table(M_df, values=['rtpeak'], index=['prot','pept','modpept','mod','[M]'],
                    columns=['charge'])
	rt_df.columns = rt_df.columns.droplevel()
	rt_df = rt_df.loc[:, np.arange(charge_from, charge_to+1)].mean(axis=1)
	rt_df = rt_df.reset_index().rename(columns={0:'rtpeak'})
	
	M_df = pd.pivot_table(M_df, values=['mz'], index=['prot','pept','modpept','mod','[M]'],
                    columns=['charge'])
	M_df.columns = M_df.columns.droplevel()
	M_df = M_df.reset_index().rename(columns={x:'Mcharge'+str(x) for x in np.arange(charge_from, charge_to+1)})
	M_df = pd.concat([M_df, rt_df.loc[:,'rtpeak']], axis=1)
	print(M_df.columns)
	iso_header = ['isoab' + str(x) for x in np.arange(iso_maxnumber)]
	print(iso_header)
	iso_df = pd.DataFrame(iso_list, columns = iso_header)
	M_df = pd.concat([M_df, iso_df], axis=1)
	print(M_df.columns)

	return M_df
 
### FUNCTIONS ###
def mass_calc(sequence, aa_monomass, other_mass):
	return sum(aa_monomass[aa] for aa in sequence) + other_mass['H20']

# def mz_calc(M_df, aa_monomass, other_mass, globalparam_list, iso_maxnumber, charge_from, charge_to):
# 	mrange_min = globalparam_list[1][globalparam_list[0].index('mrange_min')]
# 	mrange_max = globalparam_list[1][globalparam_list[0].index('mrange_max')]

# 	mrange_maxbound = int((mrange_max - (iso_maxnumber/charge_from))*charge_to)
	# for charge in range(charge_from, charge_to+1):
		# head_mz = 'Mcharge'+str(charge)
		# calc_mz = (M_df[['[M]']] + (charge*other_mass["Proton"]))/charge
		# M_df[head_mz] = calc_mz
		# # remove every charge of M_df [each row, column:head_mz] = np.NaN
		# M_df.loc[(M_df[head_mz] <= mrange_min, head_mz)] = np.NaN
		# M_df.loc[M_df[head_mz] >= mrange_max-(iso_maxnumber/charge), head_mz] = np.NaN

	# head_mz = 'Mcharge'+str(charge)
	# calc_mz = (M_df[['[M]']] + (charge*other_mass["Proton"]))/charge
	# M_df[head_mz] = calc_mz
	# # remove every charge of M_df [each row, column:head_mz] = np.NaN
	# M_df.loc[(M_df[head_mz] <= mrange_min, head_mz)] = np.NaN
	# M_df.loc[M_df[head_mz] >= mrange_max-(iso_maxnumber/charge), head_mz] = np.NaN
		
	# M_df = M_df.groupby('prot').apply(pd.DataFrame.sort_values, '[M]').reset_index(drop=True)
	# charge = charge_to - charge_from
	# # clean blank mz
	# mz_header = []
	# for charge in range(charge_from, charge_to):
	# 	mz_header.append('Mcharge'+str(charge))
	# M_df = M_df.dropna(subset = [mz_header], how='all')
	# M_df = M_df.reset_index(drop=True)
	# return M_df

def comp_calc(sequence, aa_composition):
	# apply each row function and set up a new column
	comp_dict = Counter()
	for aa in sequence:
		comp_dict.update(aa_composition[aa]) 
	return comp_dict

def iso_calc(comp_each, el_abundance, iso_maxnumber):
	def dist_coef(el, sub_num):
		if comp_each[el] >= sub_num: #combinatorics of way to sub of each el
			coef = math.factorial(comp_each[el])/(math.factorial(comp_each[el]-sub_num)*math.factorial(sub_num))
			return coef			
		else:
			return 0
	def dist_calc(possible_el, space_tosub, iso_of_el_tosub):
		abundance = 0
		for el in possible_el: #combinatorics * new iso abundance / core iso abundance
			abundance += dist_coef(el,space_tosub)*(list(el_abundance[el].values())[iso_of_el_tosub]**space_tosub/list(el_abundance[el].values())[0]**space_tosub)
		return abundance
	def dist_ofisopeak(space_tosub):
		#define el and space able to sub
		el_space = [('S',4),('S',2),('S',1),('O',2),('O',1),('H',1),('C',1),('N',1),('x',0)]
		#combinatorics of el upto number of space_tosub. 0 is necessary
		combi_all = [p for p in itertools.combinations_with_replacement(el_space, space_tosub)]
		combi_dist_calc = 0
		for combi_each in combi_all:
			if np.sum([x[1] for x in combi_each]) == space_tosub: #iso of el to sub must == space_tosub
				combi_each = [i for i in combi_each if i[1] > 0] #delete (x,0)
				combi_each = Counter(combi_each) #find how many item per el
				component_dist_calc = 1
				for k, v in combi_each.items():
					component_dist_calc *= dist_calc([k[0]],v,k[1]) #multiply each component
				combi_dist_calc += component_dist_calc
		return combi_dist_calc

	# add water
	iso_eachlist = []
	comp_each['H'] += 2
	comp_each['O'] += 1
	# calc core as product(abundance power number in comp) of each comp
	core_abundance = 1
	for comp_k, comp_v in comp_each.items():
		core_abundance *= list(el_abundance[comp_k].values())[0] ** comp_v	
	# locate core
	# M_df.loc[idx,'isoab0'] = core_abundance
	iso_eachlist.append(core_abundance)
	# calc and locate iso
	for isopeak in np.arange(iso_maxnumber-1):
		isopeak += 1 #start from 1, not calc core_abundance
		iso_eachlist.append(dist_ofisopeak(isopeak) * core_abundance)

	# iso_header = [head for head in M_df.columns if 'isoab' in head]
	# M_df['iso_sum'] = M_df[iso_header].sum(axis=1)
	# M_df = M_df.drop('comp', axis=1).reset_index(drop=True)
	return iso_eachlist

def param_autoread(M_df, iso_maxnumber, mm, tt, window, shift):
	mrange = M_df.loc[:, [x for x in M_df.columns if x[:7] == 'Mcharge']].values
	mrange = mrange[np.logical_not(np.isnan(mrange))]
	mrange_min = int(mrange.min())
	mrange_max = int(mrange.max() + 5)
	
	gradient = M_df.loc[:, 'rtpeak']
	gradient_starttime = int(np.amin(gradient))
	gradient_endtime = int(np.amax(gradient) + 2) 
	mz_range, gradient_time = mrange_max - mrange_min, gradient_endtime - gradient_starttime
	MZ_SCALE = mz_scaling(mz_range, mm)
	TIME_SCALE = time_scaling(gradient_time, tt)
	globalparam_list = [['mm', 'mrange_min', 'mrange_max', 
	                    'mz_range', 'MZ_SCALE',
	                    'tt', 'gradient_starttime', 'gradient_endtime',
	                    'gradient_time', 'TIME_SCALE', 
	                    'window', 'shift'],
	                 [mm, mrange_min, mrange_max, 
	                    mz_range, MZ_SCALE,
	                    tt, gradient_starttime, gradient_endtime,
	                    gradient_time, TIME_SCALE, 
	                    window, shift]] 

	crange = [int(x[7:]) for x in M_df.columns if x[:7] == 'Mcharge']
	charge_from = np.amin(crange)
	charge_to = np.amax(crange)
	charge_to += 1
	mz_header, charge_list = [], []
	for charge in np.arange(charge_from, charge_to):
		mz_header.append('Mcharge'+str(charge))
		charge_list.append(charge)
	iso_header = ['isoab' + str(x) for x in np.arange(iso_maxnumber)]
	print('param_autoread:', globalparam_list)
	return mz_header, iso_header, charge_list, globalparam_list


