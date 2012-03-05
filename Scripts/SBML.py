from libsbml import *
import logging

def Create_SBML_file(c, d):
	"""
	This script creates the SBML file using LibSBML. It follows the pattern of creating the
	1. Network
	2. Compartment
	3. Reactions
	4. Metabolites
	5. Modifiers
	"""

	mets = c['met_names']
	
	model, sbmlDoc  = initialise(c['Networkname'])
	
	load_species(model, mets)
	load_reactions(model, d, mets, sbmlDoc)	
	load_initials(sbmlDoc, model, c)
		
	result_string = sbmlDoc.toSBML()
	result_string = '<?xml version="1.0" encoding="UTF-8"?>\r\n' + result_string	
		
	return result_string
	
def initialise(Networkname):
	"""
	Create the SBML document, model, compartment and set basic initial values.
	"""

	sbmlDoc = SBMLDocument(2,4)
	model = sbmlDoc.createModel()
	model.setId('Metabolic_model')
	model.setName(Networkname)

	comp = model.createCompartment()
	logging.info("Compartment: %s" % comp.getName())
	comp.initDefaults()
	comp.setId('compartment')
	comp.setSize(1)
	comp.setName('cytosol')
	model.addCompartment(comp)
	
	model.setVolumeUnits('litre')
	model.setSubstanceUnits('mol')
	
	return model, sbmlDoc
	
def load_species(model, mets):
	
	for met in mets:

		sp = model.createSpecies()
		sp.setId('species_'+str(int(mets.index(met))+1))
		sp.setName(met)
		if '_x' in sp.getName(): 
			sp.setConstant(1)
			sp.setBoundaryCondition(1)
		model.addSpecies(sp)
		
def load_reactions(model, d, mets, sbmlDoc):
	"""
	Responsible for creating the reactions.
	"""
	
	for key in sorted([ int(key) for key in d.keys()]): 

		r = model.createReaction()
		r.setId('reaction_nr_'+d[str(key)]['number'])
		r.setName(d[str(key)]['label'])
		r.setCompartment('compartment')
		
		create_kinetic_law(r,d[str(key)], mets)
		
		if 'IR' not in d[str(key)]['type']: r.setReversible(1)
		else:	r.setReversible(0)
		
		subs = d[str(key)]['subs']
		for sub in subs:
			
			# procuring general species parameter
			stoc = d[str(key)]['subs_stoc'][subs.index(sub)]
			species = 'species_' + str(int(mets.index(sub))+1)
			
			# loading parameters
			r.createReactant().setSpecies(species)
			re = r.getReactant(species)
			re.setStoichiometry(int(stoc))
			re.setName(sub)
			
			if '_x' in sub: 
				r.setReversible(0)
				
		prods = d[str(key)]['prods']			
		for prod in d[str(key)]['prods']:
			
			# procuring general species parameter
			stoc = d[str(key)]['prods_stoc'][prods.index(prod)]
			species = 'species_' + str(int(mets.index(prod))+1)
			
			# loading parameters
			r.createProduct().setSpecies(species)
			re = r.getProduct(species)
			re.setStoichiometry(int(stoc))
			re.setName(prod)
			
			if '_x' in prod: 
				r.setReversible(0)				

			
def create_kinetic_law(r,data, mets):
	"""
	Firstly loads the corresponding rate law and modifys it thereafter.
	"""
	
	l = r.createKineticLaw()
	filename = sys.path[0]+'/Scripts/Rate_Laws/'+data['type']
	Rate_Law_Definition = open( filename ).read()
	
	Rate_Law_Definition = InsertSubs(data, mets, Rate_Law_Definition)

	if r.isSetReversible != False:		
		Rate_Law_Definition = InsertProds(data, mets, Rate_Law_Definition)
		
	Rate_Law_Definition = load_modifiers(Rate_Law_Definition, r,data, mets, l)
	
	if data['type'] == 'MA':	AddKeq( l, data )
		
	elif data['type'] == 'MM':
		
		AddKs( l, data, find_Ks( r, data, mets ) )
		AddKeq( l, data )
		AddFunction( l, data )
				  
	if   data['type'] == 'MMIR':	AddKs( l, data, find_Ks( r, data, mets ) )
	
	Add_VM(data, l, float(data['VM'])) 	# Add VM Parameter to rate law in libSBML
	
	l.setMath(parseFormula(Rate_Law_Definition))

	r.setKineticLaw(l)
		
def InsertSubs( data, mets, Rate_Law_Definition ):

	sub_string = ''
	
	for sub in data['subs']: 
			index = data['subs'].index(sub)
			
			species_name = ' * species_' + str(int(mets.index(sub))+1)
			
			if int(data['subs_stoc'][index]) > 1: power = '^'+str(data['subs_stoc'][index])
			else: power = ''
			sub_string += species_name + power
			
	Rate_Law_Definition = Rate_Law_Definition.replace('* substrates',sub_string)
	Rate_Law_Definition = Rate_Law_Definition.replace('substrates',sub_string[2:])	
	
	return Rate_Law_Definition	

	
def InsertProds( data, mets, Rate_Law_Definition ):

	prod_string = ''
	for prod in data['prods']: 
			
			species_name = ' * species_' + str(int(mets.index(prod))+1)

			index = data['prods'].index(prod)
			if int(data['prods_stoc'][index]) > 1: power = '^'+str(data['prods_stoc'][index])
			else: power = ''
			prod_string += species_name + power			
			
	Rate_Law_Definition = Rate_Law_Definition.replace('* products',prod_string)
	Rate_Law_Definition = Rate_Law_Definition.replace('products', prod_string[2:])
	
	return Rate_Law_Definition
	
		
def AddFunction(RateLaw, data):
	"""
	This function returns the value for the denominator function. 
	It's provisorical yet.
	"""
	
	f = RateLaw.createParameter()
	f.setName('f_function')
	f.setId('f_function')
	f.setValue(1) # provisorical
	RateLaw.addParameter(f)
	
			
def AddKeq( RateLaw, data ):
	
	Keq = RateLaw.createParameter()
	Keq.setName('Keq')
	Keq.setId('Keq')
	Keq.setValue( float(data['KEQ']) )
	Keq.setUnits("dimensionless")
	RateLaw.addParameter(Keq)
	
def AddKs(RateLaw, data, value):

	Ks = RateLaw.createParameter()
	Ks.setName('Ks')
	Ks.setId('Ks')
	Ks.setValue( value )
	Ks.setUnits("dimensionless")
	RateLaw.addParameter(Ks)
		
		
def find_Ks(r, data, mets):
	
	value = 1.0 # for multiplication
	for sub in data['subs']: value = value * float(data['Km'][sub])
	
	return value

def load_initials(sbmlDoc, model, c):
	
	for species in sbmlDoc.getModel().getListOfSpecies():
		species.setInitialConcentration(float(c['Conc_Dict'][species.getName()]))	
		species.setCompartment('compartment')		

def Add_VM( data, l, value):		
	
	VM = l.createParameter()
	VM.setName('VM')
	VM.setId('VM')		
	VM.setValue(value)		
	VM.setUnits("dimensionless")
	l.addParameter(VM)		

		
def load_modifiers(formula, r, data, mets, l):
	
	act_string = ''
	
	for mod_name in data['activators']:
		
		index_sbml = str( int(mets.index(mod_name)) + 1 ) 
		r.createModifier().setSpecies('species_'+index_sbml)
		mod = r.getModifier('species_'+index_sbml)
		mod.setName(mod_name)	
		
		Sbml_Name = 'species_'+index_sbml
		index_reaction = data['activators'].index(mod_name) 
		nh_value = data['activators_nh'][index_reaction]
		Power = '('+Sbml_Name+' / Kax_'+ mod_name +' )^nh_'+ mod_name
		act_string += '( '+ Power +' / ( 1 + '+ Power +' ) ) *'
		
		add_nh( nh_value, mod_name, l )

		add_Mod( mod_name, l, 'Kax_', float(data['KREG_KM'][mod_name]) )
		
	formula = formula.replace('Kact *', act_string)[:-1] # get rid of * char

	inh_string = ''
	
	for mod_name in data['inhibitors']:
		
		index_sbml = str( int(mets.index(mod_name)) + 1 ) 
		r.createModifier().setSpecies('species_'+index_sbml)
		mod = r.getModifier('species_'+index_sbml)
		mod.setName(mod_name)	
		
		Sbml_Name = 'species_'+index_sbml
		index_reaction = data['inhibitors'].index(mod_name) 
		nh_value = data['inhibitors_nh'][index_reaction]
		Power = '(' + Sbml_Name + '/ Kix_'+ mod_name +' )^nh_'+ mod_name
		inh_string += '(1 / ( 1 + ' + Power + ' ) ) *'
		
		add_nh( nh_value, mod_name, l )
		add_Mod( mod_name, l, 'Kix_', float(data['KREG_KM'][mod_name]))
	
	formula = formula.replace('Kinh *', inh_string)

	return formula
	
def add_nh( nh_value, mod_name, RateLaw):
	
	nh = RateLaw.createParameter()
	nh.setName('nh_'+mod_name)
	nh.setId('nh_'+mod_name)
	nh.setValue( float(nh_value) )
	nh.setUnits("dimensionless")
	RateLaw.addParameter(nh)	

def add_Mod( mod_name, RateLaw, ActInh, value ):
	
	Mod = RateLaw.createParameter()
	Mod.setName( ActInh+mod_name )
	Mod.setId( ActInh+mod_name)
	Mod.setValue( value )
	Mod.setUnits( "dimensionless" )
	RateLaw.addParameter( Mod )
	
