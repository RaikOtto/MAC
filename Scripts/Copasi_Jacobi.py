from COPASI import *
import sys

def Show_Copasi_Jacobi( c ):
	
	Networkname = c['Networkname']
	dataModel = CCopasiRootContainer.addDatamodel()
	result=True
	assert CCopasiRootContainer.getRoot() != None
	assert dataModel != None
	assert CCopasiRootContainer.getDatamodelList().size() == 1
	
	result = dataModel.importSBML(Networkname+"/"+Networkname+".xml")
	
	CCopasiMessage.clearDeque()
	result=True


	if result != True and  mostSevere < CCopasiMessage.ERROR:
		print >> sys.stderr, "Sorry. Model could not be imported."
		return

	model = dataModel.getModel()
	assert model != None
	
	if model != None:
		model.applyInitialValues()
		jacobian=FloatMatrix()
		model.calculateJacobian(jacobian, 1e-12, 1.0)
		stateTemplate = model.getStateTemplate()
		userOrder = stateTemplate.getUserOrder()

		nameVector=[]
		entity = None
		status=-1

		for i in range(0,userOrder.size()):
			entity = stateTemplate.getEntity(userOrder.get(i))
			assert entity != None

			status = entity.getStatus()

			if (status == CModelEntity.ODE or (status == CModelEntity.REACTIONS and entity.isUsed())):
				nameVector.append(entity.getObjectName())

		assert len(nameVector) == jacobian.numRows()

		print "Jacobian Matrix:"
		print ""
		print "%7s" % (" "),

		for i in range(0,len(nameVector)):
			print "%7s" % (i+1),

		print ""

		for i in range(0,len(nameVector)):
			print "%7s" % (nameVector[i]),

			for j in range(0,len(nameVector)):
				print "%7.3f" % (jacobian.get(i,j)),

			print ""
