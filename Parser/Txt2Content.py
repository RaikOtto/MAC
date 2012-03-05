import re, Defaults, sys

def CreateParaFile(raw_data, c, d):
	"""
	This function creates the Matlab Parameter File. All default 
	definitions are to be found in the default.py file.
	"""

	Default_Reac_List = Defaults.Default_Reac_List
	data = trimm_data(raw_data, Default_Reac_List) # responsible for default values,
		#e.g. MM for all non-default reactions !!! Default values are in Defaults.py
	d = load_dictionary(d,data) # loads trimmed data into dictionary
	d = parse_reaction(d) # determins the reactions
	c = parse_names(c,d) # loads the sub and prod names
	d = parse_modifiers(c,d, Default_Reac_List) # responsible for activators and inhibitors, as well as nh-values
	d = load_default_KM(d) # checks as well for external metabolite -> irreversible reaction type
	c = build_NT(c,d) # finds KNR and KNS as well, KNR = diff. kin. coefficient from stoch.
	c = find_KMR_KMM(c,d) # KM values for all functions
	c = find_KREGs(c,d) # which reactions are regulated
	d = calculate_stoch_ks(c,d)
	
	return c, d
	
def trimm_data(raw_data, Default_Reac_List):
	"""
	Cuts the strings read from the txt into proper chunks of clean pieces.
	"""
	data = []
	
	Default_Reac_Type = Defaults.ReacType
	Default_Keq_Value = Defaults.Keq

	for line in raw_data:
		
		rtype_found = False
		keq_found   = False
		km_found	= False
		split_mods_line = []
		line = line.strip()
		# ignore empty and comment lines
		if 	len(line) == 0 or line == '\r\n' or line.startswith('%') or line == '\r': continue
		# checks if reaction has a name
		elif len(line.split(':')) == 1: line = 'no label: ' + line
		# split line up into description and modifiers	
		splitline = line.split(';')
		
		if len(splitline) > 1: # if len == 0 => no modifiers
			split_mods_line = splitline[1].split(' ') # split modifiers
			
			for mod in split_mods_line:
				
				mod = mod.strip()
				if mod == '' or mod.startswith('C='): continue
				
				elif mod.startswith('-n') or mod.strip().startswith('+n'): 
					split_mods_line[split_mods_line.index(mod)] = ' '+mod.strip()
				
				elif mod.startswith('-'): 
					split_mods_line[split_mods_line.index(mod)] = ' -n4'+mod.strip()[1:]
				
				elif mod.startswith('+'): 
					split_mods_line[split_mods_line.index(mod)] = ' +n4'+mod.strip()[1:]
					
				elif mod in Default_Reac_List: rtype_found = True
				
				elif mod.startswith('KEQ='): keq_found = True
				
				elif mod.startswith('KM='): km_found = True
				
				elif mod.startswith('VM'): continue
				
				elif mod[0] == ('c'): continue
				
				else: print 'Error, unknown modifier: '+ mod + '!\r\n'
		
		if not keq_found:   split_mods_line.append(Default_Keq_Value)
		if not rtype_found: split_mods_line.append(Default_Reac_Type)
		
		splitline = splitline[0] + '; ' # entry zero contains the reaction and the name, the ; separates the mods
		
		for entry in split_mods_line:
			
			if entry == ' ' or entry == '' : continue
			splitline += ' '  + entry
		
		data.append(splitline)
		
	return data
	
def load_dictionary(d, data):
	
	i = 1
		
	for line in data:
		d[str(i)] = {}
		split = line.split(':')
		d[str(i)]['label'] = split[0]
		content = split[1].split(';')
		d[str(i)]['number'] = str(i)
		d[str(i)]['reaction'] = content[0].strip()
		
		if len(content) >1:
			d[str(i)]['modifiers'] = content[1]
		else:
			d[str(i)]['modifiers'] = None
		i += 1
		
	return d
	
def parse_reaction(d):
	"""
	This function changes the arrows in the reaction definitions, as well as it initilizes the reactions in the d dictionary
	"""

	for key in d.keys():
		
		initilize_d(d, key)
		
		met_external = False

		if '=>' in d[key]['reaction']: 
		
			d[key]['reaction'] = d[key]['reaction'].replace('=>','->')
			
		if '->' in d[key]['reaction']: content = d[key]['reaction'].strip().split('->')
		else: content = d[key]['reaction'].strip().split('=')
		
		for elem in content[0].strip().split(' '):
			
			try:
				int(elem.replace("'",''))
			except:
				if elem != '+' and elem != '-':
					if len(elem) > 0:
						if d[key].has_key('subs'): d[key]['subs'].append(elem)
						else: d[key]['subs'] = [elem]
						
		
		for elem in content[1].strip().split(' '):
			
			try:
				int(elem.replace("'",''))
			except:
				if elem != '+' and elem != '-':
					if len(elem) > 0:
						if d[key].has_key('prods'): d[key]['prods'].append(elem)
						else: d[key]['prods'] = [elem]
		
	return d
	
def parse_modifiers(c, d, Default_Reac_List):
	"""
	This function searches for all know types of modificators loads the information into the dictionary
	"""
	
	Default_Conc = Defaults.Conc
	
	for key in d.keys():
		
		if d[key]['modifiers'] != None:
			
			for mod in d[key]['modifiers'].split(' '):
				
				mod = mod.strip()
				
				if len(mod) == 0: continue
				
				if mod.startswith('-'): 
				
					inhibitorname = find_name(c['met_names'],mod)
					if not	d[key]['KREG_KM'].has_key(inhibitorname): 
							d[key]['KREG_KM'][inhibitorname] = {}
					
					if '_' in mod: 
					
						KM = re.search("_\S*",mod).group(0).replace('_','')
						c['KREG_KM'].append(KM)
						mod = mod.replace('_'+KM,'')
						d[key]['KREG_KM'][inhibitorname] = KM
						c['KREG_KM'][inhibitorname].append(KM)
		
					else: 
						d[key]['KREG_KM'][inhibitorname] = Defaults.KREG_KM
						c['KREG_KM'].append(Defaults.KREG_KM)					
					
					d[key]['inhibitors'].append(inhibitorname)
					d[key]['inh_reg'].append(str(c['met_names'].index(inhibitorname)+1))
					d[key]['inhibitors_nh'].append(mod.replace(inhibitorname,'')[2:])

				elif mod.startswith('+'): 
				
					activatorname = find_name(c['met_names'],mod)
					if not	d[key]['KREG_KM'].has_key(activatorname): 
							d[key]['KREG_KM'][activatorname] = {}					
					
					if '_' in mod: 
					
						KM = re.search("_\S*",mod).group(0).replace('_','')
						mod = mod.replace('_'+KM,'')
						d[key]['KREG_KM'][activatorname] = KM
						c['KREG_KM'].append(KM)
						
					else: 
						d[key]['KREG_KM'][activatorname] = Defaults.KREG_KM
						c['KREG_KM'].append(Defaults.KREG_KM)
						
					d[key]['activators'].append(activatorname)
					d[key]['act_reg'].append(str(c['met_names'].index(activatorname)+1))
					d[key]['activators_nh'].append(mod.replace(activatorname,'')[2:])

				elif mod.startswith('KEQ='): d[key]['KEQ'] = mod.replace('KEQ=','')
				
				elif mod.startswith('VM='): d[key]['VM'] = mod.replace('VM=','')
				
				elif mod.startswith('KM='): 
					string 	= mod.replace('KM=','') # met and Km number
					Km_dict = d[key]['Km']

					name, value = fill_Km( string, Km_dict)
					d[key]['Km'][name] = value
										
				elif mod in Default_Reac_List: d[key]['type'] = mod
				
				elif mod.startswith('C='):
					
					met  = mod.replace('C=','')
					name = find_name(c['met_names'],met)
					
					if c['Conc_Dict'].has_key(name) is not True:
						c['Conc_Dict'][name] = met.replace( name+'_', '')
					else: print 'Warning! Double definition for reaction nr.: ' + str(key) + ' metabolite ' + name + '!\r\n'
			d = Set_Defaults( d, key)
	
	for entry in c['met_names']:
		if c['Conc_Dict'].has_key(entry) is not True: c['Conc_Dict'][entry] = Default_Conc
	
	return d
	
def fill_Km( string, Km_dict):
	
	name  = re.search( '\w*_', string).group(0)[:-1]
	value = re.search( '_\w*', string).group(0)[1:]
	Km_dict[name] = value
	
	return name, value
	
def load_default_KM(d):
	
	for key in d.keys():
		
		met_external = False
		
		for met in d[key]['subs']:	
			if not d[key]['Km'].has_key(met): d[key]['Km'][met] = Defaults.Km
			if '_x' in met: met_external = True; met_found = met
			
		for met in d[key]['prods']: 
			if '_x' in met: met_external = True; met_found = met
		
		#if met_external and 'IR' not in d[key]['type']:
		#	print 'External metabolite ('+met_found+') found in reaction nr. '+key+', changing reaction type ' +\
		#		d[key]['type'] + ' to ' + d[key]['type'] + 'IR!'
		#	d[key]['type'] = d[key]['type']+'IR'
		#	d[key]['reaction'] = d[key]['reaction'].replace('=','->')
		
	return d
	
def initilize_d(d, key):
	
	d[key]['inhibitors'] 	= []
	d[key]['inh_reg'] 		= []
	d[key]['inhibitors_nh'] = []
	d[key]['activators'] 	= []
	d[key]['act_reg'] 		= []
	d[key]['activators_nh'] = []
	d[key]['Km']  			= {}
	d[key]['KREG_KM']		= {}
	d[key]['VM']			= Defaults.VM
	
	return d
	
def Set_Defaults( d, key):
	
	Default_Reac = Defaults.ReacType
	
	if not d[key].has_key('inhibitors'): d[key]['inhibitors'] = None
	if not d[key].has_key('activators'): d[key]['activators'] = None
	if d[key].has_key('type') is not True: 
		if '->' not in d[key]['reaction'] : d[key]['type'] = Default_Reac
	
	extern = False
	for sub in d[key]['subs']:	
		if '_x' in sub: extern = True
	
	if not 'IR' in d[key]['type']:
		
		if '->' in d[key]['reaction'] or extern: 
			d[key]['type'] = d[key]['type'].strip() + 'IR'
			print 'Warning, changing reaction Nr.'+str(key)+ ' to irreversible!'	
			
	del d[key]['modifiers']
	
	return d
	
		
def find_KREGs(c,d):
	"""
	Identifies what reactions are regulated by activators and inhibitors
	"""
	
	KREG_R, KREG_M, KREG_nh, KREG_T = [], [], [], []
	
	for i in range(1,int(c['nr_reactions'])+1):
		
		if d[str(i)]['activators'] != None:
			
			for met in d[str(i)]['activators']:
				
				activator = find_name(c['met_names'],met)
				KREG_R.append(str(i))
				KREG_M.append(str(c['met_names'].index(activator)+1))
				KREG_nh.append(str(d[str(i)]['activators_nh']\
				[d[str(i)]['activators'].index(activator)]))
				KREG_T.append('+1')
				
		if d[str(i)]['inhibitors'] != None:
			
			for met in d[str(i)]['inhibitors']:
				
				inhibitor = find_name(c['met_names'],met)
				KREG_R.append(str(i))
				KREG_M.append(str(c['met_names'].index(inhibitor)+1))
				KREG_nh.append(str(d[str(i)]['inhibitors_nh']\
				[d[str(i)]['inhibitors'].index(inhibitor)]))
				KREG_T.append('-1')
				
	c['KREG_R']  = KREG_R
	c['KREG_M']  = KREG_M
	c['KREG_nh'] = KREG_nh
	c['KREG_T']  = KREG_T
	
	return c
	
	
def find_name(namelist, met):
	
	best_match = ''
	
	for name in namelist:
		match = re.search(name, met)
		
		if '_sre.SRE_Match' in str(type(match)):
			if len(match.group(0)) > len(best_match):	
				best_match = match.group(0)
				
	return best_match
	
	
def find_KMR_KMM(c,d):
	"""
	Finds KMR and KMM values
	
	It will be neccessary to modify this step once different types of
	reactions different from MM and MA are introduced
	"""
	
	KMR, KMM = [], []

	for i in range(1,int(c['nr_reactions'])+1):
		
		for t in range(len(d[str(i)]['subs'])):
			
			KMR.append(str(i))		
			KMM.append(c['met_names'].index(d[str(i)]['subs'][t])+1)
			
		if 'IR' not in d[str(i)]['type']:
			
			for n in range(len(d[str(i)]['prods'])):
				
				KMR.append(str(i))
				KMM.append(c['met_names'].index(d[str(i)]['prods'][n])+1)
		
	c['KMR'], c['KMM'] = KMR, KMM
	
	return c
	
	
def build_NT(c,d):
	"""
	Here the stoichiometry matrix is created.
	"""
	
	unsortKNR = []
	unsortKNS = []
	c['KNR_found'] = False
	
	line = [ ' 0' for entry in range(int(c['nr_reactions'] ))]
	NT 	 = [line[:] for i in range(c['nr_dyn'] + c['nr_ext'] )]
	
	for elem in c['met_names']:

		for key in d.keys():
			
			if elem in d[key]['reaction']:
				
				content = d[key]['reaction'].split(' ')
				
				for i in range(0,len(content)):
					
					if content[i] == elem:		
					
						try:
							
							if "'" in content[i-1]:
								
								c['KNR_found'] = True
								unsortKNR.append(key)
								unsortKNS.append(str(c['met_names'].index(elem)+1))
								content[i-1] = content[i-1].replace("'",'')
								
							coeff = str(int(content[i-1]))
							
						except:	coeff = '1'			
										
						column = int(d[key]['number']) - 1
						row = c['met_names'].index(elem)
						
						if elem in d[key]['subs']: NT[row][column] = ' -'+coeff
						else: NT[row][column] = ' +'+coeff
						
	c['NT'] = NT
	c['KNR'], c['KNS'] = Find_KNR_KNS(unsortKNR,unsortKNS)
	
	return c
	
	
def Find_KNR_KNS(unsortKNR,unsortKNS):

	if unsortKNR != []: # at least one reaction with dif. stoch from kin. coeff found
	
		sortKNR  , sortKNS 	 = [unsortKNR[0]], [unsortKNS[0]]
		unsortKNR, unsortKNS =  unsortKNR[1:],  unsortKNS[1:]
		
		for i in range(0,len(unsortKNR)):
			
			KNR, KNS = unsortKNR[i], unsortKNS[i]
			for t in range(len(sortKNR)-1,-1,-1):
				
				if int(KNR) < int(sortKNR[t]):
					tmpKNR, tmpKNS = sortKNR[t], sortKNS[t]
					sortKNR[t], sortKNS[t] = KNR, KNS
					
					try:
						sortKNR[t+1], sortKNS[t+1] = tmpKNR, tmpKNS
					except:
						sortKNR.append(tmpKNR)
						sortKNS.append(tmpKNS)
						
				else:
					if t+1 == len(sortKNR):
						sortKNR.append(KNR)
						sortKNS.append(KNS)
						break
						
	# no reaction with dif. stoch from kin. coeff found
	else: sortKNR, sortKNS = [], []
		
	return sortKNR, sortKNS
	
	
def parse_names(c,d):
	"""
	Out of a string this function finds the best match for a metabolite in a random string out of the list of the known metabolites.
	"""
	
	sortarray  = []
	
	for key in d.keys():
		for sub in d[key]['subs']: 
			if sub not in sortarray: sortarray.append(sub)
		for prod in d[key]['prods']: 
			if prod not in sortarray: sortarray.append(prod)

	met_names = sorted(sortarray)
	dyn_list, ext_list = [], []
	
	for elem in met_names:
		if '_x' in elem: 
			ext_list.append(elem)
		else: dyn_list.append(elem)
	
	c['met_names'] 	= dyn_list + ext_list
	c['nr_dyn']		= len(dyn_list)
	c['nr_ext']		= len(ext_list)
	c['nr_reactions'] = max([int(key) for key in d.keys() ])
	
	return c
	
	
def calculate_stoch_ks(c,d):
	
	km_counter = 0
	
	knr, kns = c['KNR'], c['KNS']
	knr_dict = {}
	i = 0
	
	while i < len(knr):
		
		if not knr_dict.has_key(knr[i]): knr_dict[knr[i]] = [kns[i]]
		else: knr_dict[knr[i]].append(kns[i])
		i += 1
		
	for reac in sorted([int(key) for key in d.keys()]):
		
		reac = str(reac)
		d[reac]['subs_stoc'], d[reac]['prods_stoc'] = [], []
		d[reac]['subs_ks']  , d[reac]['prods_ks']   = [], [] 
		if not knr_dict.has_key(reac) : knr_dict[reac] = []
		
		for sub in d[reac]['subs']:
			
			sub_index = str(c['met_names'].index(sub) + 1)
			if sub_index not in knr_dict[reac]:	
			
				nt_entry = c['NT'][int(sub_index)-1][int(reac)-1]
				nt_entry = str(abs(int(nt_entry)))
				d[reac]['subs_stoc'].append( nt_entry )
				
			else: d[reac]['subs_stoc'].append('1')
				
			km_counter += 1	
			d[reac]['subs_ks'].append('KM('+str(km_counter)+')')
			
		for prod in d[reac]['prods']:
			
			prod_index = str(c['met_names'].index(prod) +1)
			if prod_index not in knr_dict[reac]: 
				
				nt_entry = c['NT'][int(prod_index)-1][int(reac)-1]
				nt_entry = str(abs(int(nt_entry)))
				d[reac]['prods_stoc'].append( nt_entry )
				
			else: d[reac]['prods_stoc'].append('1')

			if 'IR' not in d[reac]['type']:
				
				km_counter += 1
				d[reac]['prods_ks'].append('KM('+str(km_counter)+')')
	
	c['KM_count'] = km_counter
	
	return d
