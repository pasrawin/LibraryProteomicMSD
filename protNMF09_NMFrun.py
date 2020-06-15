from __future__ import division
import numpy as np
import pandas as pd
from scipy import sparse
from scipy import special
import timeit

### FUNCTIONS ###
from mNMF05_ExpmatConstruction import smoothingtime_mat
from mNMF10_NMFReport import nmf_identification, nmf_savemat

def tol_test(it, convergence, tracking_it_all):
	tol_value = 1e-4		
	if it > 8: # error tolerance (run at least 4 + 1 it)
		tol = (tracking_it_all[-2] - tracking_it_all[-1])/tracking_it_all[-2]
		if tol <= tol_value:
			convergence = True #go to last it
			print ' report convergence at it: ', it, tol, tracking_it_all
		else:
			print ' report at it: ', it, tol, tracking_it_all
	return convergence

def kl_calc(Vmat, Wmat, Hmat, Vmat_mean, peptcount_eachsqrt, eachprot_sum_it, weight, eps):
	if weight != 0:
		penaltyH_it = np.sum(np.multiply(peptcount_eachsqrt,np.log(np.asarray(eachprot_sum_it) + eps))) #log L1 group
		penaltyH_it = weight*penaltyH_it
		# print ' chk main penaltyH_it:', penaltyH_it
	else:
		penaltyH_it = 0

	V_approx = Wmat.dot(Hmat) + eps
	cost_it = special.kl_div(Vmat.todense(),V_approx)
	cost_it = np.sum(cost_it) + penaltyH_it
	del(V_approx)	
	return cost_it

def kl_nmf(Vmat, Wmat, Hmat, cos_mat, prot_peptcount,\
        weight, max_iter, eps, noise_number, noise_mean, Vmat_mean, globalparam_list, iso_maxnumber,\
        initRT, outfile, finalMS1_df, peak_report, WithStandardize_):
	print ' process NMF'
	for param in globalparam_list:
		mm, mrange_min, mrange_max, mz_range, MZ_SCALE, \
		tt, gradient_starttime, gradient_endtime, gradient_time, TIME_SCALE, \
		window, shift =  globalparam_list[1]

	print ' shapes:', Vmat.shape, Wmat.shape, Hmat.shape
	peptcount_each = prot_peptcount.values
	# ioncount_each = pept_ioncount.values
	# print ioncount_each
	if WithStandardize_:
		peptcount_eachsqrt = peptcount_each**0.5 # for inside standardize
		# ioncount_eachsqrt = ioncount_each**0.5 
	else:
		peptcount_eachsqrt = [1 for x in peptcount_each]
	peptcount_cumsum = np.insert(np.cumsum(peptcount_each), 0, 0)
	# ioncount_cumsum = np.insert(np.cumsum(ioncount_each), 0, 0)
	# protcount, peptcount = len(peptcount_each), len(ioncount_each)
	protcount = len(peptcount_each)

	it_all, tracking_it_all = [], []
	eachprot_timeprofile = []
	eachprot_sum, eachprot_mean, eachprot_std = [], [], []

	Wmat = Wmat.tocsr()
	Vmat = sparse.csr_matrix(Vmat)

	if sparse.issparse(Hmat):
		refresh_H = Hmat.copy().todense()
	else:
		refresh_H = Hmat.copy()

	convergence = False
	# dimension arrangement lkij
	# print ' ***constrain by dividing pept'
	for it in xrange(max_iter):	
		it_all.append(it)
		numer = Wmat.dot(refresh_H) + eps
		numer = Wmat.transpose().dot(Vmat/numer)
		refresh_H_copy = refresh_H.copy()
		refresh_H = np.multiply(refresh_H, numer)
		del(numer)
		eachprot_sum_it = []
		# for each in range(protcount):
		for each in range(protcount):
			from_k = peptcount_cumsum[each]
			to_k = peptcount_cumsum[each+1]
			refresh_Hg = refresh_H[from_k:to_k,slice(None)].copy()
			if weight != 0:
				deno = 1 + ((weight*peptcount_eachsqrt[each])/(np.sum(abs(refresh_Hg))+eps)) # log l1 group 
				refresh_H[from_k:to_k,:] = refresh_Hg/deno
				del(refresh_Hg)
				del(deno)
				eachprot_sum_it.append(np.sum(refresh_H[from_k:to_k,:]))
				# eachprot_sumax_it.append(np.sum(refresh_H[from_k:to_k,:],axis=1))
				if each == protcount-1:
					cost_it = kl_calc(Vmat, Wmat, refresh_H, Vmat_mean, peptcount_eachsqrt, eachprot_sum_it, weight, eps) # calc last before clean
					tracking_it_all.append(cost_it)
					convergence = tol_test(it, convergence, tracking_it_all)
		if convergence:
			break

	for each in range(protcount):
		from_k = peptcount_cumsum[each]
		to_k = peptcount_cumsum[each+1]
		refresh_Hg = refresh_H[from_k:to_k,slice(None)].copy()
		refresh_Hg_process = (refresh_H[from_k:to_k,slice(None)]*Vmat_mean)
		eachprot_sum.append(np.sum(refresh_Hg_process))
		eachprot_mean.append(np.mean(refresh_Hg_process))
		eachprot_timeprofile.append(np.ravel(np.sum(refresh_Hg_process,axis=0)))
	eachprot_timeprofile.append(np.ravel(np.sum(refresh_H[-noise_number:,slice(None)],axis=0))*Vmat_mean) #append noise

	print ' process report',' report time: ', timeit.default_timer()
	it_all = pd.DataFrame(np.column_stack([np.arange(len(tracking_it_all)),tracking_it_all]), columns=['it', 'tracking'])
	nmf_identification(outfile, Vmat*Vmat_mean, Wmat, sparse.csr_matrix(refresh_H*Vmat_mean), finalMS1_df,\
						globalparam_list, prot_peptcount, [Wmat.tocsc()[slice(None),-noise_number:].todense()], noise_mean, weight,\
						it_all, initRT, eachprot_sum, eachprot_mean, eachprot_timeprofile,\
						iso_maxnumber, cos_mat, peak_report)
	print ' saving matrices with renormalization by Vmat_mean:', Vmat_mean
	nmf_savemat(Vmat*Vmat_mean, Wmat, sparse.csr_matrix(refresh_H*Vmat_mean), weight, outfile)