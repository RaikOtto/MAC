import re

def Create_Matlab_Para_File(c,d):
	
	text = ''
	text = Met_Names_Dyn_ext(text, c, d)
	text = Reaction_types(text, c, d)
	text = NT_matrix(text, c , d)

	text += '\r\n'
	text += 'N = NT(1:length(metabolite_dyn),:);\r\n'
	text += 'K = null(N);\r\n'
	text += '\r\n'
	text = Remaining_Information(text,c,d)
	
	return text
	
def Remaining_Information(text,c,d):
	"""
	Fills information for KNR, KNS, KMR, KMM, KREG_R, KREG_M, KREG_nh,
	 KREG_T, SX and Keq. SQ and Keq not final
	"""
	
	if c['KNR_found']:
		text += 'KNR = [ '
		for item in c['KNR']: text += str(item) + ' '
		text += ']; % Which reac. has different kin to stoch. coefficient\r\n'
		text += 'KNS = [ '
		for item in c['KNS']: text += str(item) + ' '
		text += ']; % What met. has different kin to stoch. coefficient\r\n'
		text += '\r\n'
	
	
	text += 'KREG_R = [ '
	for item in c['KREG_R']: text += str(item) + ' '
	text += ']; % What reaction is modified\r\n'
	text += 'KREG_M = [ '
	for item in c['KREG_M']: text += str(item) + ' '
	text += ']; % What metabolite serves as a modifier \r\n'
	text += 'KREG_nh = [ '
	for item in c['KREG_nh']: text += str(item) + ' '
	text += ']; % strength of modification \r\n'
	text += 'KREG_T = [ '
	for item in c['KREG_T']: text += str(item) + ' '
	text += ']; % Modifier activator or inhibitor \r\n'
	
	return text
	
def build_KREG(text, c):
	"""
	Separated from the other KREG vector since it is used in the simulation file
	"""
	
	text += 'KREG = [ '
	for item in c['KREG_KM']: text += str(item) + ' '
	text += ']; % KM value for modification \r\n'	
	
	return text
	
def NT_matrix(text, c, d):
	"""
	Function only works for st. coeff up to size 9
	"""
	text += '\r\n%       '
	for t in range(c['nr_reactions']): text += str(t+1) + ' '
	text += '\r\n'
	text += 'NT = [ '

	for i in range( len(c['NT']) ):
		line = '       '
		if i == 0:	line = ''
		for item in c['NT'][i]: 
			spacer = ''
			for r in range(2-len(str(item))): spacer += ' '
			line += spacer + str(item)
		if i != len(c['NT'])-1:	text += line + ' ; %'+ c['met_names'][i] +'\r\n'
		else: text += line + ']; %'+ c['met_names'][i] +'\r\n'
	return text
	
def Reaction_types(text,c,d):

	text += '\r\n'
	text += 'REACTION_TYPES = ['
	for i in range(1,len(d.keys())+1):
		text += " '" + d[str(i)]['type'] +"'"
	text = text[:-1] # get rid of the last " ;"
	text += "' ];\r\n"
	
	return text
	
def Met_Names_Dyn_ext(text,c,d):
	
	text += 'metabolite_names = ['
	
	for name in c['met_names']:
		text += " '" + name + "'"
		
	text += ' ]; \r\n'
	text += 'metabolite_dyn = ['
	
	for i in range(1,c['nr_dyn']+1):
		text += ' '+str(i)
	text += ' ];\r\n'
	
	text += 'metabolite_ext = ['
	
	for i in range(c['nr_ext']):
		text += ' '+str(c['nr_dyn']+1+i)
	text += ' ];\r\n'
	
	return text
