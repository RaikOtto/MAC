import os

model_files_path = "/home/carpathy/Dropbox/ITB/MAC/ModelFiles"
mac_path = "/home/carpathy/Dropbox/ITB/MAC"


def clean_up(model, mac_path):
	
	if os.path.isfile(model): os.system('rm '+model)
	
	model = model.replace('.txt','')
	if os.path.isfile(model+'_octave'): os.system('rm '+model+'_octave')
	
	if os.path.exists(model):
		
		os.system('rm '+model+'/*')
		os.system('rmdir '+model)	
		
def octave_test(model):
	
	model = model.replace('.txt','')
	
	text = 'cd '+model+' \r\n'
	text += model + '_Sim'+'\r\n'
	text += 'cd ../'+' \r\n'
	text += 'exit'+' \r\n'
	
	fileh = open(model+'_octave','w').write(text)
	
	os.system('octave '+model+'_octave')
		

for model in os.listdir(model_files_path):
	
	print '- '+model+' -'
	clean_up(model, mac_path)
	os.system('cp ModelFiles/'+model +' '+ model)
	
	os.system('python2.6 MAC.py '+model)
	#octave_test(model)
	
	clean_up(model, mac_path)
	
	
