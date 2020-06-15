from __future__ import division
import numpy as np
import pandas as pd
from pyteomics import achrom

def rt_prediction(calib_file, finalMS1_df, RT_equation):	
	print ' from:', calib_file
	calib_df = pd.read_excel(calib_file)
	initRT_tuple = []
	if RT_equation != False:		
		finalMS1_df['rtpeak_calib'] = finalMS1_df['rtpeak'].apply(RT_equation)
		# print ' adjust rt with minus 30 s (1 dynamic exclusion)'
		for idx, row in finalMS1_df.iterrows():
			pept = str(row['pept'])
			rt = row['rtpeak_calib']
			initRT_tuple.append((idx,pept,rt))
	else:
		for idx, row in finalMS1_df.iterrows():
			pept = str(row['pept'])
			rt = row['rtpeak']
			initRT_tuple.append((idx,pept,rt))
		 
	try:
		reten_start = calib_df['PepRtimeStart']
		reten_end = calib_df['PepRtimeEnd']
		initRT_width = np.mean((reten_end - reten_start).values)
	except:
		initRT_width = 1
		print ' warnings: using default initRT_width:', initRT_width
	
	return initRT_tuple, initRT_width



		