from __future__ import division
import numpy as np
import pandas as pd
from collections import Counter
import math
import itertools
# import timeit
import gc

### MAIN ###
def insilico_detect(mascot_file, dictionary_list, param_list, ms1param_list, iso_maxnumber):
	gc.disable()
	print ' from:', str(mascot_file)
	mascot_df = pd.read_excel(mascot_file)
	aa_monomass, aa_composition, aa_fixmodmass, aa_varmodmass, other_mass, el_abundance = dictionary_list
	misclvge_from, misclvge_to, len_from, len_to, charge_from, charge_to = param_list
	# mm, mrange_min, mrange_max, mz_range, MZ_SCALE, tt, gradient_starttime, gradient_endtime, \
	#	 gradient_time, TIME_SCALE, window, shift = ms1param_list[1]
	mm = ms1param_list[1][ms1param_list[0].index('mm')]
	mass_decimal = int(np.log10(mm))
	mod_mass = aa_varmodmass.values()[0][0]
	prot_list, seq_list, mod_list, modseq_list, mass_list, rt_list, comp_list = [], [], [], [], [], [], []
	mascotunique_dict = {}
	for index, row in mascot_df.iterrows():
		protein = row['AccNum'] 
		sequence = row['Seq']
		comp_ans = comp_calc(sequence, aa_composition)
		mod = str(row['Mod'])
		try:
			rt = row['PepRtimePeak correct']
		except:
			rt = row['PepRtimePeak']
		key = (protein, sequence, mod)
		if key not in mascotunique_dict:
			if 'Oxidation (M)' in str(mod):
				mod_mass = 15.99491
				ox_position = str(row['ModDetail']).split(',')
				ox_position = [x for x in ox_position if '@M:' in x]
				ox_position = [int(x.split('@M:')[1]) for x in ox_position]
				ox_number = len(ox_position)
				comp_ans['O'] += ox_number
				start_section, sequence_section_keep = 0, ''			
				for ox_position_each in ox_position:
					ox_position_each -= 1
					# end_section = start_section+ox_position_each
					sequence_section = sequence_section_keep + sequence[start_section:ox_position_each] + 'x'
					sequence_section_keep = sequence_section
					start_section = ox_position_each
				sequence_section_keep = sequence_section_keep + sequence[ox_position[-1]-1:]
				modseq_list.append(sequence_section_keep)				
				mass_list.append(mass_calc(sequence, aa_monomass, other_mass) + (mod_mass*ox_number))
				mod_list.append(mod[:-4]+'@M:'+str(ox_position)[1:-1])
			else:
				mod_list.append('')
				modseq_list.append(sequence)
				mass_list.append(mass_calc(sequence, aa_monomass, other_mass))
			prot_list.append(protein)
			seq_list.append(sequence)
			rt_list.append(rt)
			comp_list.append(comp_ans)
			mascotunique_dict[key] = 1

	M_df = pd.DataFrame({'prot': prot_list,
			'pept': seq_list,
			'modpept': modseq_list,
			'mod': mod_list,
			'[M]': mass_list,
			'comp': comp_list,
			'rtpeak':rt_list
			})
	# print(mascot_df.shape, M_df.shape)
	M_df = M_df[['prot','pept','modpept','mod','[M]','comp','rtpeak']]
	M_df = mz_calc(M_df, aa_monomass, other_mass, ms1param_list, iso_maxnumber, charge_from, charge_to)
	M_df = iso_calc(M_df, el_abundance, iso_maxnumber)
	mz_header, iso_header, charge_list = label_df(iso_maxnumber, charge_from, charge_to)
	return M_df, mz_header, iso_header, charge_list
 
### FUNCTIONS ###
def mass_calc(sequence, aa_monomass, other_mass):
	return sum(aa_monomass[aa] for aa in sequence) + other_mass['H20']

def mz_calc(M_df, aa_monomass, other_mass, globalparam_list, iso_maxnumber, charge_from, charge_to):
	mrange_min = globalparam_list[1][globalparam_list[0].index('mrange_min')]
	mrange_max = globalparam_list[1][globalparam_list[0].index('mrange_max')]

	mrange_maxbound = int((mrange_max - (iso_maxnumber/charge_from))*charge_to)
	for charge in range(charge_from, charge_to+1):
		head_mz = 'Mcharge'+str(charge)
		calc_mz = (M_df[['[M]']] + (charge*other_mass["Proton"]))/charge
		M_df[head_mz] = calc_mz
		# remove every charge of M_df [each row, column:head_mz] = np.NaN
		M_df.loc[(M_df[head_mz] <= mrange_min, head_mz)] = np.NaN
		M_df.loc[M_df[head_mz] >= mrange_max-(iso_maxnumber/charge), head_mz] = np.NaN
		
	M_df = M_df.groupby('prot').apply(pd.DataFrame.sort_values, '[M]').reset_index(drop=True)
	charge = charge_to - charge_from
	# clean blank mz
	mz_header = []
	for charge in range(charge_from, charge_to):
		mz_header.append('Mcharge'+str(charge))
	M_df = M_df.dropna(subset = [mz_header], how='all')
	M_df = M_df.reset_index(drop=True)
	return M_df

def comp_calc(sequence, aa_composition):
	# apply each row function and set up a new column
	comp_dict = Counter()
	for aa in sequence:
		comp_dict.update(aa_composition[aa]) 
	return comp_dict

def iso_calc(M_df, el_abundance, iso_maxnumber):
	def dist_coef(el, sub_num):
		if comp_each[el] >= sub_num: #combinatorics of way to sub of each el
			coef = math.factorial(comp_each[el])/(math.factorial(comp_each[el]-sub_num)*math.factorial(sub_num))
			return coef			
		else:
			return 0
	def dist_calc(possible_el, space_tosub, iso_of_el_tosub):
		abundance = 0
		for el in possible_el: #combinatorics * new iso abundance / core iso abundance
			abundance += dist_coef(el,space_tosub)*(el_abundance[el].values()[iso_of_el_tosub]**space_tosub/el_abundance[el].values()[0]**space_tosub)
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
				for k, v in combi_each.iteritems():
					component_dist_calc *= dist_calc([k[0]],v,k[1]) #multiply each component
				combi_dist_calc += component_dist_calc
		return combi_dist_calc
	
	for idx, comp_each in M_df.loc[:,['comp']].itertuples():
		# add water
		comp_each['H'] += 2
		comp_each['O'] += 1
		# calc core as product(abundance power number in comp) of each comp
		core_abundance = 1
		for comp_k, comp_v in comp_each.iteritems():
			core_abundance *= el_abundance[comp_k].values()[0] ** comp_v	
		# locate core
		M_df.loc[idx,'isoab0'] = core_abundance
		# calc and locate iso
		for isopeak in range(iso_maxnumber-1):
			isopeak += 1 #start from 1, not calc core_abundance
			M_df.loc[idx,'isoab'+str(isopeak)] = dist_ofisopeak(isopeak) * core_abundance

	iso_header = [head for head in M_df.columns if 'isoab' in head]
	M_df['iso_sum'] = M_df[iso_header].sum(axis=1)
	M_df = M_df.drop('comp', axis=1).reset_index(drop=True)
	return M_df

def label_df(iso_maxnumber, charge_from, charge_to):
	charge_to += 1
	mz_header, charge_list = [], []
	for charge in range(charge_from, charge_to):
		mz_header.append('Mcharge'+str(charge))
		charge_list.append(charge)
	iso_header = ['isoab' + str(x) for x in range(iso_maxnumber)]
	return mz_header, iso_header, charge_list


