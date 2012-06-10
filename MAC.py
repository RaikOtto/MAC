import sys, Parser, re ,os, Scripts, traceback

def CreateMatlabFile(filename):
	"""
	This is the main script file. It coordinates all processes. All 
	five text files are created by their own script files in the script 
	directory. MAC.py calls these scripts.
	"""
	
	### preparation ###
	networkname = re.sub("\.\S*",'', filename)
	osfilename = sys.path[0] + '/' + networkname	
	d = {}			 													# Reaction dictionary, only for reactions
	c = {'Networkname': networkname, 'KREG_KM': [], 'Conc_Dict' : {}}	# Content dictionary, for background information
	raw_data = [line for line in open(osfilename+'.txt').readlines()]
	
	try:	content = Parser.Txt2Content.CreateParaFile( raw_data, c ,d)
	except: print traceback.print_exc(10); sys.exit('Error parsing the txt! Execution aborted!')
		
	#content = Scripts.Check_Content.Check_Content(content[0],content[1]) yet 2 do
	c 		= content[0]
	d 		= content[1]
	
	# Create the files
	#Para_file   = Scripts.MatlabParaFile.Create_Matlab_Para_File(c, d)
	Matlab_file = Scripts.MatlabFile.Write_MatlabFile(c, d)
	Jacobi_file = Scripts.Jacobi_Matrix.Write_Jacobi_MatlabFile(c, d)
	Octave_file = Scripts.OctaveFile.Create_Octave_file(c,d)
		
	#Create directories if they don't exist
	if os.path.exists(osfilename): os.system('rm -rf '+osfilename+'/*')
	else: os.makedirs(osfilename)
	
	# Move files
	open(osfilename+'/'+networkname+'.m','w').write(Matlab_file)
	open(osfilename+'/'+networkname+'_Jacobi.m','w').write(Jacobi_file)
	open(osfilename+'/'+networkname+'_Sim.m','w').write(Octave_file)
	
	try: # check if there is LibSBML installed:
	
		import libsbml
		SBML_file   = Scripts.SBML.Create_SBML_file(c, d)
		open(osfilename+'/'+networkname+'.xml','w').write(SBML_file)
		#Scripts.Copasi_Jacobi.Show_Copasi_Jacobi(c)
		
	except:
		
		print 'LibSBML not found, no SBML file created.'
		error_message = traceback.format_exc(10)		
		open('Error_information.txt','w').write( error_message )
	
if __name__ == "__main__":
	
	if len(sys.argv[:]) == 1: raise Exception('Need txt file')
	CreateMatlabFile(sys.argv[1]) 
