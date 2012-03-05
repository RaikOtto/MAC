import functions as f, MatlabFile as MF

def Write_Jacobi_MatlabFile(c, d):
	"""	
	This script writes the Jacobi_Matrix file. It iterates over the reactions and creates their derivates according to the metabolites. Metabolites are either included in the reaction itself or appear as inhibition or activiation factors. Products in irreversible reactions are not included, under the condition that they do not appear as activators in inhibitors.
	"""
			
        text = write_opener(c, d)
        text = write_concentrations(c,text)
        text = write_JM_initiator(c,text)

        for reac in sorted([int(key) for key in d.keys()]): # the idea is to use the matlabfile method 
            reac = str(reac)				    			# to avoid having two definitions of write_reactions
            text = text + MF.write_reactions(c, {reac:d[reac]}, '').replace("* REGA",'').replace("* REGI",'')
				# reference to MatlabFile.py -> writes the part which is identical
            text = write_JM_entries(c,d[reac],text)
    
        text += '\r\n'
        text += 'J = N*EE;\r\n'
	
	return text
		
def write_opener(c,d ):
	
	text  = '%%% '+c['Networkname']+'_Jacobi %.%%\r\n'
	text += '\r\n'
	text += 'function [J,EE] = '+c['Networkname']+'_Jacobi(N,concentration_variable)\r\n'
	text += '\r\n'
	text += 'global para'	+c['Networkname']+'\r\n'
	text += 'VM = para'		+c['Networkname']+'.VM;\r\n'
	
	types = [reac_type for reac_type in [d[key]['type'] for key in d.keys() ] ]
	if 'MM' in types or 'MMIR' in types:	text += 'KM = para'		+c['Networkname']+'.KM;\r\n'
	
	text += 'S_0 = para'	+c['Networkname']+'.S0;\r\n'
	text += 'Keq = para'	+c['Networkname']+'.Keq;\r\n'
	text += 'KREG = para'	+c['Networkname']+'.KREG;\r\n'
	text += 'KREG_nh = para'+c['Networkname']+'.KREG_nh;\r\n'

	text += '\r\n'
	
	return text

def write_concentrations(c,text):
	
    text += '% assign concentration variables\r\n'
    ext_temp_text = '\r\n% external concentrations (constant metabolites)\r\n' #because external concentrations come later
    
    for met in c['met_names']: 
		if '_x' in met:	ext_temp_text += met + ' = S_0(' + str(c['met_names'].index(met)+1-c['nr_dyn']) + ');\r\n'
		else: text += met + ' = concentration_variable(' + str(c['met_names'].index(met)+1) + ');\r\n'
    text += ext_temp_text
    
    text += '\r\n'
    
    return text	
    
def write_JM_initiator(c,text):

	text += '% Initiate Jacobian matrix\r\n'
	text += 'EE = zeros('+str(c['nr_reactions'])+','+ str(c['nr_dyn']) +');\r\n'
	
	return text    
	
def write_JM_entries(c,d,text):
	"""
	This function creates the jacobi matlab commands according to whether a metabolite is a product/ substrate and/ or 
	activator/ inhibitor. All combinations are possible: e.g. no product, substrate, but inhibitor.
	"""
	
	inhibs 	  = d['inhibitors']
	activs    = d['activators']
	subs      = d['subs']
	prods     = d['prods']
	r_type    = d['type']
	
	mets      = subs   + prods
	mods      = inhibs + activs
	
	for mod in mods:
		if mod not in mets: mets += [mod] # as to avoid double entries
	
	for met in mets: 
		
		if '_x' in met: continue 
		# no external mets in jacobi, not even as inhibitor or activator
		
		index_met = str(mets.index(met) + 1 )
		cord 	  = str(c['met_names'].index(met) + 1 )
		reac_num  = d['number']
		
		if int(index_met) > len(subs):
			
			index_met = str( (int(index_met) - len(subs) * -1 ) )
			
			
		### Function of the metabolite in the reaction ###		
			
		# substrate or product in a reversible reaction
		if met in subs or (met in prods and 'IR' not in r_type):
			
			if met in subs: index_sub_or_prod = str(subs.index(met) + 1)
			else: 			index_sub_or_prod = str( -1*(prods.index(met) + 1) )
			text = prep_text(text, met, inhibs, activs, reac_num, r_type, index_sub_or_prod)
			
			# case not activator and not inhibitor but product or substrate in a reversible reaction
			if met not in inhibs and met not in activs:
				text = sub_or_prod_and_not_active_and_not_inhibit(text, met, inhibs, activs, reac_num, cord)
			
			# case activator		
			elif met not in inhibs and met in activs:
				text = sub_or_prod_and_active(text, met, inhibs, activs, reac_num, cord)
				
			# case inhibitor
			elif met in inhibs and met not in activs:
				text = sub_or_prod_and_inhib(text, met, inhibs, activs, reac_num, cord)
				
			# case activator and inhibitor
			else:
				text = sub_or_prod_and_active_and_inhibit(text, met, inhibs, activs, reac_num, cord)
		
		# just inhibitor and/or met product and reaction irreversible or metabolite not in subs and prods
		elif met in inhibs and met not in activs:
			text = only_inhib_or_active_and_met_inhib(text, met, inhibs, activs, reac_num, cord)
		
		# just activator and/or met product and reaction irreversible or metabolite not in subs and prods
		elif met not in inhibs and met in activs:
			text = only_inhib_or_active_and_met_active(text, met, inhibs, activs, reac_num, cord)
			
		# activator and inhibitor and/or met product and reaction irreversible or metabolite not in subs and prods
		elif met in inhibs and met in activs:
			text = only_inhib_and_active_and_met_is_both(text, met, inhibs, activs, reac_num, cord)
		
		# case metabolite products in IR reaction ant not inhibitor and not activator	
		elif met not in inhibs and met not in activs: pass
						
		else: print 'Error! Metabolite ' + met + ' cant be attributed to any role in the reaction '+ reac_num +'!\r\n'
	
	return text
	
def prep_text(text, met, inhibs, activs, reac_num, r_type, index_sub_or_prod):
		"""
		if   inhibs == [] and activs == []:
			text += eval('f.'+r_type+'_deriv').replace('$',reac_num).replace("&",index_sub_or_prod)+';\r\n'
		elif inhibs != [] and activs == []:
			text += eval('f.'+r_type+'_deriv').replace('$',reac_num).replace("&",index_sub_or_prod)+' * REGI;\r\n'
		elif inhibs == [] and activs != []:
			text += eval('f.'+r_type+'_deriv').replace('$',reac_num).replace("&",index_sub_or_prod)+' * REGA;\r\n'
		else:
			text += eval('f.'+r_type+'_deriv').replace('$',reac_num).replace("&",index_sub_or_prod)+' * REGA * REGI;\r\n'
		
		return text
		"""
		return text + eval('f.'+r_type+'_deriv').replace('$',reac_num).replace("&",index_sub_or_prod)+';\r\n'

def sub_or_prod_and_not_active_and_not_inhibit(text, met, inhibs, activs, reac_num, cord):
	
	if inhibs == [] and activs == []:
		text += 'EE('+ reac_num +','+ cord +') = dvds;%'+met+'\r\n'
		
	elif inhibs != [] and activs == []:
		text += 'EE('+ reac_num +','+ cord +') = dvds * REGI;%'+met+'\r\n'
		
	elif inhibs == [] and activs != []:
		
		text += 'EE('+ reac_num +','+ cord +') = dvds * REGA;%'+met+'\r\n'
		
	else:
		text += 'EE('+ reac_num +','+ cord +') = dvds * REGI * REGA;%'+met+'\r\n'
		
	return text
	
def sub_or_prod_and_active(text, met, inhibs, activs, reac_num, cord):
	
	index_activs = str( activs.index(met) + 1 )
	text += f.Activate.replace("v($) = VM($) *","").replace("&",index_activs)+';\r\n'
	
	if inhibs == []:
		text += 'EE('+ reac_num +','+cord+') = dvds * REGA + v(' + reac_num + ') * dfda;%'+met+'\r\n'
		
	else:
		text += 'EE('+ reac_num +','+cord+') = dvds * REGA * REGI + v(' + reac_num + ') * dfda * REGI;%'+met+'\r\n'
	
	return text
	
def sub_or_prod_and_inhib(text, met, inhibs, activs, reac_num, cord):
	
	index_inhibs = str( inhibs.index(met) + 1 )
	text += f.Inhibit.replace("v($) = VM($) *","").replace("&",index_inhibs)+';\r\n'	
	
	if activs == []:	
		text += 'EE('+ reac_num +','+cord+') = dvds * REGI + v(' + reac_num + ') * dfdi;%'+met+'\r\n'	
			
	else:
		text += 'EE('+ reac_num +','+cord+') = dvds * REGA * REGI + v(' + reac_num + ') * dfdi * REGA;%'+met+'\r\n'
		
	return text
	
def sub_or_prod_and_active_and_inhibit(text, met, inhibs, activs, reac_num, cord):
	
	index_activs = str( activs.index(met) + 1 )
	index_inhibs = str( inhibs.index(met) + 1 )
	text += f.Inhibit.replace("v($) = VM($) *","").replace("&",index_inhibs)+';\r\n'		
	text += 'EE('+ reac_num +','+cord+') = dvds * REGI + v(' + reac_num + ') * dfdi;%'+met+'\r\n'	
	text += f.Activate.replace("v($) = VM($) *","").replace("&",index_activs)+';\r\n'
	text += 'EE('+ reac_num +','+cord+') = dvds * REGA + v(' + reac_num + ') * dfda;%'+met+'\r\n'
	text += 'EE('+ reac_num +','+cord+') = '+\
		'dvds * REGI * REGA + '+\
		'v(' + reac_num + ') * dfdi + '+\
		'v(' + reac_num + ') * dfda;%'+met+'\r\n'
	
	return text	
	
def only_inhib_or_active_and_met_active(text, met, inhibs, activs, reac_num, cord):
	
	index_activs = str( activs.index(met) + 1 )
	text += f.Activate.replace("v($) = VM($) *","").replace("&",index_activs)+';\r\n'
	text += 'EE('+ reac_num +','+cord+') = v(' + reac_num + ') * dfda;%'+met+'\r\n'	
	
	if inhibs != []: 
		text = text[:-5]
		text = text[:-1] + text[-1:].replace(';','')
		text = text +  ' * REGI;%'+met+'\r\n'	
			
	return text
	
def only_inhib_or_active_and_met_inhib(text, met, inhibs, activs, reac_num, cord):
	
	index_inhibs = str( inhibs.index(met) + 1 )
	text += f.Inhibit.replace("v($) = VM($) *","").replace("&",index_inhibs)+';\r\n'
	text += 'EE('+ reac_num +','+cord+') = v(' + reac_num + ') * dfdi;%'+met+'\r\n'
	
	if activs != []: 
		text = text[:-5]
		text = text[:-1] + text[-1:].replace(';','')
		text = text +  ' * REGA;%'+met+'\r\n'	
		
	return text
	
def only_inhib_and_active_and_met_is_both(text, met, inhibs, activs, reac_num, cord):
	
	index_inhibs = str( inhibs.index(met) + 1 ); index_activs = str( activs.index(met) + 1 )
	text += f.Activate.replace("v($) = VM($) *","").replace("&",index_activs)+';\r\n'
	text += f.Inhibit.replace( "v($) = VM($) *","").replace("&",index_inhibs)+';\r\n'
	text += 'EE('+ reac_num +','+cord+') = v(' + reac_num + ') * dfdi * REGA + v(' + reac_num + ') * dfda * REGI ;%'+met+'\r\n'
		
	return text
