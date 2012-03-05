import re, MatlabParaFile as MP

def Create_Octave_file(c,d):
	
	text = MP.Create_Matlab_Para_File(c,d)
	
	text = write_sim_header(c,d,text)
	
	return text
		
def write_sim_header(c,d, text):
	
	title =  ' Simulation of '+c['Networkname'] + ' '
	length_title = len(title)
	spacer = '%%%' + ''.join(['%' for i in range(length_title)]) + '%%%\r\n'
	text += '\r\n'
	text += spacer + '%%%' + title + '%%%\r\n' + spacer
	
	text += 'source("../Scripts/functions.m")\r\n'
	text += '\r\n'

	text  = Add_default_values(text,c,d)
	text += '\r\n'
		

	text += 'T = [0:0.1:50];\r\n'
	text += '[S] = lsode("'+c['Networkname']+'", X0, T);\r\n'
	
	text += 'plot(T,[S]);\r\n'
	text += '%axis([0,25,-2,2]);\r\n'
	text += 'title("'+c['Networkname']+'");\r\n'
	text += 'xlabel("time");\r\n'
	text += 'ylabel("concentration");\r\n'
	
	text += 'legend( '
	i = 0
	while i < c['nr_dyn']: text += '"'+c['met_names'][i]+'", '; i += 1
	text = text[:-2] + ' );\r\n'
	#text = text + '[J,EE] = '+c['Networkname']+'_Jacobi(N,X0);\r\n'
	
	return text
	
def Add_default_values(text,c,d):
	
	# SX and Keq part not final
	
	nr_reacs = c['nr_reactions']
	
	text += '\r\n% Default values %\r\n\r\n'
	
	keq_list = [kq for kq in [ d[str(entry)]['KEQ'] for entry in range(1,nr_reacs+1) ]]
	text += 'Keq = [ '+ ' '.join(keq_list) +'];\r\n'
	
	# filter KM if not needed
	types = [reac_type for reac_type in [d[key]['type'] for key in d.keys() ] ]
	if 'MM' in types or 'MMIR' in types:
		Kms = []
		for item in range(0,int(c['KM_count'])): Kms.append('1')
		
		KM = []
		for item in range(0, c['KM_count'] ): KM.extend(['1'])
	
		for reac in d.keys():
			
			for sub in d[reac]['subs']:
				
				if d[key]['Km'].has_key(sub):
					
					index = int(re.search('\(\w*',d[reac]['subs_ks'][0]).group(0)[1:]) - 1
					KM[index] = d[key]['Km'][sub]
					
			for prod in d[reac]['prods']:
				
				if d[key]['Km'].has_key(prod) and 'IR' not in d[reac]['type']:
					index = int(re.search('\(\w*',d[reac]['prods_ks'][0]).group(0)[1:]) - 1
					KM[index] = d[key]['Km'][prod]

		KM_vec = ' '
		for item in KM: KM_vec = KM_vec + ' ' + str(item)
		text += 'KM = ['+ KM_vec[2:] + '];\r\n'
	
	# X0 and S0 defaults
	X0 = 'X0 = [ '
	S0 = 'S0 = [ '
	# VM values
	text += 'V0  = [ '
	reacs = [ str(key) for key in sorted( [ int(key) for key in d.keys() ] ) ]
	for reac in reacs: text += d[reac]['VM'] + ' '
	text += '];\r\n\r\n'
	
	for met in c['met_names']:
		
		text += met + ' = ' + str(c['Conc_Dict'][met]) + ';\r\n'
		if '_x' not in met:	X0 += met + ' '
		else: S0 += met + ' '
		
	text += '\r\n' + X0 + '];\r\n'
	text += S0 + '];\r\n\r\n'
	
	text += 'if (min([X0,S0]) <0), printf("ERROR, negative concentration: %f\\n",min([X0,S0]));end;\r\n'
	
	text += '\r\n'
	text += 'global para'+c['Networkname']+'\r\n'
	text += 'para'+c['Networkname']+'.Keq = Keq;\r\n'
	text += 'para'+c['Networkname']+'.KREG = KREG;\r\n'
	text += 'para'+c['Networkname']+'.KREG_nh = KREG_nh;\r\n'
	text += 'para'+c['Networkname']+'.N  = N;\r\n'
	text += 'para'+c['Networkname']+'.VM  = V0;\r\n'
	text += 'para'+c['Networkname']+'.S0  = S0;\r\n'
	types = [reac_type for reac_type in [d[key]['type'] for key in d.keys() ] ]
	if 'MM' in types or 'MMIR' in types:	text += 'para'+c['Networkname']+'.KM = KM;\r\n'	
	text += '\r\n'
	
	
	# VM values
	text += 'VM  = ' + c['Networkname'] + '(X0,V0); % normalize VM parameter\r\n\r\n'
	
	# VM checks
	text = write_VM_checks(text)
	
	text += 'para'+c['Networkname']+'.VM  = VM;\r\n'	
	text += '\r\n% End default values %\r\n'
	
	return text
	
def write_VM_checks(text):
	
	text += '% check for VM < 1, these have to be provided individually\r\n'
	text += 'if (find(VM==-2)), printf(" ERROR, The metabolic state is inconcistent");end;\r\n'
	text += 'if (find(VM==-1)), printf(" The flux vector contains equilibrium reactions: Specify Vmax manually");end;\r\n'
	
	return text
