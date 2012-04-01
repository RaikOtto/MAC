import re, MatlabParaFile as MP

def Create_Octave_file(c,d):

	text = write_sim_file_header(c, d)	
	text += MP.Create_Matlab_Para_File(c,d)+ '\r\n\r\n'
	text = move_KREG( text, True, False, c )
	
	text += Add_default_values( c, d )	
	text = write_sim_without_normalisation(c, d, text)
	text = write_sim_with_normalisation(c, d, text)
	
	text += create_title(' Monte-Carlo Simulations with scaled Parameters ')
	text += '\r\n[J,EE] = '+c['Networkname']+'_Jacobi(N,X0); \r\n'
	text += create_title('% End of simulations %')
	
	return text
		
def write_sim_file_header(c, d):
	
	# just creates the comment titles in a convenient way
	text = create_title(' Network information '+c['Networkname'] + ' ')
	text += 'source("../Scripts/functions.m")\r\n'
	text += 'global para'+c['Networkname']+'\r\n\r\n'
	
	return text
	
def Add_default_values( c, d ):
	
	# SX and Keq part not final
	
	nr_reacs = c['nr_reactions']
	text = ''
	
	text += '% Default values %\r\n\r\n'
	
	text = move_KREG( text, False, True, c )	
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
	text += 'V0  = [ ' #abs(K);\r\n'#
	reacs = [ str(key) for key in sorted( [ int(key) for key in d.keys() ] ) ]
	for reac in reacs: text += d[reac]['VM'] + ' '
	text += '];\r\n\r\n'
	
	text += 'V0 = K; \r\n' ### ADDED !!! ###
	
	for met in c['met_names']:
		
		text += met + ' = ' + str(c['Conc_Dict'][met]) + ';\r\n'
		if '_x' not in met:	X0 += met + ' '
		else: S0 += met + ' '
	
	text += '\r\n% End default values %\r\n'
		
	text += '\r\n' + X0 + '];\r\n'
	text += S0 + '];\r\n\r\n'
	
	text += Add_Para_links( c )
	
	types = [reac_type for reac_type in [d[key]['type'] for key in d.keys() ] ]
	if 'MM' in types or 'MMIR' in types:	text += 'para'+c['Networkname']+'.KM = KM;\r\n'		
	text += '\r\n'

	return text	
	
def Add_Para_links( c ):
	
	text = 'para'+c['Networkname']+'.Keq = Keq;\r\n'
	text += 'para'+c['Networkname']+'.KREG = KREG;\r\n'
	text += 'para'+c['Networkname']+'.KREG_nh = KREG_nh;\r\n'
	text += 'para'+c['Networkname']+'.VM  = V0;\r\n'
	text += 'para'+c['Networkname']+'.S0  = S0;\r\n'
	text += 'para'+c['Networkname']+'.X0  = X0;\r\n'

	return text

def write_sim_without_normalisation(c,d, text):
	
	text += create_title('% Simulation %')

	text += 'subplot(2,1,1)\r\n'
	text += 'T = [0:0.1:50];\r\n\r\n'
	text += '[S] = lsode("'+c['Networkname']+'", X0, T);\r\n\r\n'
	text += 'plot(T,[S]);\r\n'
	text += 'title("'+c['Networkname']+' without normalised parameters");\r\n'
	text += 'xlabel("time");\r\n'
	text += 'ylabel("concentration");\r\n'
	
	text += 'legend( '
	i = 0
	while i < c['nr_dyn']: text += '"'+c['met_names'][i]+'", '; i += 1
	text = text[:-2] + ' );\r\n\r\n\r\n'
	
	return text
	
	
def write_sim_with_normalisation(c,d, text):

	text += create_title('% Simulation with scaled parameters %')
	
	text += filter_scaled_parameters( Add_default_values( c, d ) )
	
	# Normalise VM values
	text += 'para'+c['Networkname'] + '.VM = '+c['Networkname'] + '(X0,V0); % normalize VM parameter\r\n\r\n'

	# VM checks
	text = write_VM_checks(text, c)	
	
	text += 'subplot(2,1,2)\r\n'
	text += 'T = [0:0.1:50];\r\n\r\n'
	text += '[S] = lsode("'+c['Networkname']+'", 0.5 * X0, T);\r\n\r\n'
	text += 'plot(T,[S]);\r\n'
	text += 'title("'+c['Networkname']+' with normalised parameters");\r\n'
	text += 'xlabel("time");\r\n'
	text += 'ylabel("concentration");\r\n'

	text += 'legend( '
	i = 0
	while i < c['nr_dyn']: text += '"'+c['met_names'][i]+'", '; i += 1
	text = text[:-2] + ' );\r\n\r\n\r\n'
	text += 'print("jpg","plot.jpg");\r\n'
	return text	
	
		
def create_title(title):
	
	length_title = len(title)
	spacer = '%%%' + ''.join(['%' for i in range(length_title)]) + '%%%\r\n'
	
	return spacer + '%%%' + title + '%%%\r\n' + spacer + '\r\n'
	
	
def write_VM_checks( text, c):
	
	text += '% check for VM < 1, these have to be provided individually\r\n'
	text += 'if (find(para'+c['Networkname']+'.VM==-2)), printf(" ERROR, The metabolic state is inconcistent\\r\\n");end;\r\n'
	text += 'if (find(para'+c['Networkname']+'.VM==-1)), printf(" The flux vector contains equilibrium reactions: Specify Vmax manually\\r\\n");end;\r\n'
	text += '%if ( length(find( round( abs(N*V0\') * 10^15 ) ==0)) != size(N)(1) ), printf("Warning: V0 vector not in the kernel of the matrix\\r\\n");end;\r\n\r\n'
	
	return text
	
	
def move_KREG( text, delete, move, c ):
	"""
	This helper function moves the KREG vectors within the simulation file since it is needed at a different place than in the matlab file.
	"""
	
	relocate = MP.build_KREG( text, c )
	
	if delete : text = text.replace(relocate,'')
	if move: text = relocate	
	
	return text
	
def filter_scaled_parameters( text ):
	"""
	KREG_nh, Keq and their links to global para are not needed and are filtered
	"""

	#text = re.sub('Keq.*', '' , text)
	text = re.sub('.*Keq.*\r\n', '' , text)
	text = re.sub('.*KREG_nh.*\r\n', '' , text)
	
	return text
