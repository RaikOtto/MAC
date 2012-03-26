import functions as f


def Write_MatlabFile(c, d):
			
	text = write_opener(c, d)
	text = write_concentrations(c,text)			
	text = write_reactions(c, d, text)
	text = assign_stoichiometry(c,text)
	
	return text
	
		
def write_opener(c, d):
	
	text  = '%%% '+c['Networkname']+' %%%\r\n'
	text += '\r\n'
	text += 'function dSdt = '+c['Networkname']+'(concentration_variable,t)\r\n'
	text += '\r\n'
	text += 'global para'	+c['Networkname']+'\r\n'
	text += 'VM = para'		+c['Networkname']+'.VM;\r\n'
	
	types = [reac_type for reac_type in [d[key]['type'] for key in d.keys() ] ]
	if 'MM' in types or 'MMIR' in types:	text += 'KM = para'	+c['Networkname']+'.KM;\r\n'
		
	text += 'S_0 = para'	+c['Networkname']+'.S0;\r\n'
	text += 'Keq = para'	+c['Networkname']+'.Keq;\r\n'
	text += 'KREG = para'	+c['Networkname']+'.KREG;\r\n'
	text += 'KREG_nh = para'+c['Networkname']+'.KREG_nh;\r\n\r\n'
	text += 'if (length(t)>1); VM=ones(1,length(VM)); V0 = t;end;\r\n'
	text += '\r\n'
	return text
	

def write_concentrations(c,text):
	
    text += '% assign concentration variables\r\n'
    ext_temp_text = '\r\n% external concentrations (constant metabolites)\r\n' #because external concentrations come later
    for met in c['met_names']: 
		if '_x' in met:	ext_temp_text += met + ' = S_0(' + str(c['met_names'].index(met)+1-c['nr_dyn']) + ');\r\n'
		else: text += met + ' = concentration_variable(' + str(c['met_names'].index(met)+1) + ');\r\n'
    text += ext_temp_text
    return text	
	
def write_reactions(c, d, text):
	"""
	This function writes all .m File reactions.
	"""	
	
	KREG_counter, KREG_nh_counter = 0, 0
	
	for reac in sorted([int(key) for key in d.keys()]):
		
		reac     = str(reac)		
		reaction = eval('f.'+d[reac]['type']+".replace('$',reac)") 
		
		# here the predefined reactions get used, taken from function.py
		r_type   = d[reac]['type'] # reaction type, e.g MMIR		
		text  += '\r\n%reaction nr.: '+ reac + ' (' + d[reac]['reaction'] + ')\r\n'
		
		# substrates
		text += 'S = [ '+" ".join( d[reac]['subs']) + ' ]; NS = [ ' + " ".join( d[reac]['subs_stoc'] ) + ' ]; '
		if r_type != 'MA' and r_type != 'MAIR': 
			text += 'KS = [ ' + ' '.join( d[reac]['subs_ks'] ) + ' ];'
		text += '\r\n'
		
		# products
		if 'IR' not in r_type:
			text += 'P = [ '+" ".join( d[reac]['prods'] )+' ]; NP = [ ' + " ".join( d[reac]['prods_stoc'] ) + ' ]; '
			if r_type != 'MA': 
				text += 'KP = [ '+ ' '.join( d[reac]['prods_ks'] ) +' ];'
			text += '\r\n'
	
		text = add_modifiers( text, c, d, reac, KREG_counter, KREG_nh_counter, reaction  )
		
	return text
	
	
def add_modifiers( text, c, d, reac, KREG_counter, KREG_nh_counter, reaction):
	

	## modifiers ##

	# case only inhibitors 
	if d[reac]['activators'] == [] and d[reac]['inhibitors'] != []:	

		text += 'I = [ '+(' ').join( d[reac]['inhibitors'] )+' ]; KIX =[ '
		for item in d[reac]['inhibitors']: KREG_counter += 1; text += 'KREG('+str(KREG_counter)+') '
		text += ']; nh = [ ' 
		for item in d[reac]['inhibitors_nh']: KREG_nh_counter += 1; text += 'KREG_nh('+ str(KREG_nh_counter) +') '
		text += '];\r\n'
		text += f.Rate_Inhibit
		text += reaction +' * REGI;\r\n' 
		
	# Case only activators	
	elif d[reac]['activators'] != [] and d[reac]['inhibitors'] == []:	

		text += 'A = [ '+(' ').join( d[reac]['activators'] )+' ]; KAX =[ '
		for item in d[reac]['activators']: KREG_counter += 1; text += 'KREG('+str(KREG_counter)+') '
		text += ']; nh = [ ' 
		for item in d[reac]['activators_nh']: KREG_nh_counter += 1; text += 'KREG_nh('+ str(KREG_nh_counter) +') '
		text += '];\r\n'
		text += f.Rate_Activate
		text += reaction +' * REGA;\r\n' 	
		
	# Case activators and inhibitors		
	elif d[reac]['activators'] != [] and d[reac]['inhibitors'] != []: 

		text += 'A = [ '+(' ').join( d[reac]['activators'] )+' ]; KAX =[ '
		for item in d[reac]['activators']: KREG_counter += 1; text += 'KREG('+str(KREG_counter)+') '
		text += ']; nh = [ ' 
		for item in d[reac]['activators_nh']: KREG_nh_counter += 1; text += 'KREG_nh('+ str(KREG_nh_counter) +') '
		text += '];\r\n'
		
		text += 'I = [ '+(' ').join( d[reac]['inhibitors'] )+' ]; KIX =[ '
		for item in d[reac]['inhibitors']: KREG_counter += 1; text += 'KREG('+str(KREG_counter)+') '
		text += ']; nh = [ ' 
		for item in d[reac]['inhibitors_nh']: KREG_nh_counter += 1; text += 'KREG_nh('+ str(KREG_nh_counter) +') '
		text += '];\r\n'
		
		text += f.Rate_Act_Inh
		text += reaction + ' * REGA * REGI; \r\n'

	# Case no activators and inhibitors	
	else: text += reaction + ';\r\n'

	text = text + '\r\n'
			
	return text	
	
def assign_stoichiometry(c, text):
	
	text += '\r\n% ASSIGN STOICHIOMETRY \r\n'
	mets_dyn = len([met for met in c['met_names'] if '_x' not in met])
	
	for i in range(1,mets_dyn+1):

		text += 'dSdt('+str(i)+') = '
		NTline = [entry for entry in c['NT'][i-1]]
		for n in range(0,len(NTline)):
			if int(NTline[n]) != 0:
				text += NTline[n] + ' * v(' + str(n+1) + ') '
		text += ';\r\n'
	
	text += '\r\n'
	
	text += 'if (length(t)>1);\r\n'
	text += 'ix = find(V0==0);\r\n'
	text += 'dSdt = V0./v;\r\n'
	text += 'dSdt(ix) = -1; % for reactions in equilibrium VM is assigned differently\r\n'
	text += 'end;\r\n'
	
	return text
