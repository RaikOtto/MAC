
def Check_Content(c, d):
	
	d = check_for_MM(d)
	return (c,d)
	
def check_for_MM(d):

	for key in sorted([int(entry) for entry in d.keys()]):
		key = str(key)
		error = False
		if d[key]['type'] == 'MM' or d[key]['type'] == 'MMIR':
			if len(d[key]['subs']) != 1 or len(d[key]['prods']) != 1 and d[key]['type'] != 'MMIR': 
				error = True
			else:
				for stoc in d[key]['subs_stoc'] + d[key]['prods_stoc']:
					if abs(int(stoc)) != 1: error = True
		if error == True and d[key]['type'] == 'MM':
			print 'Error: Changing Reaction Nr.'+key+' ('+d[key]['reaction']+') from type MM to MA!'
			d[key]['type'] = 'MA' 
		elif error == True and d[key]['type'] == 'MMIR':
			print 'Error: Changing Reaction Nr.'+key+' ('+d[key]['reaction']+') from type MMIR to MAIR!'
			d[key]['type'] = 'MAIR' 
	return d
