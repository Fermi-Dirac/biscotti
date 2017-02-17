from classes import qecalc
newcontrol  = {'title' : 'demo',
               'calculation' : 'relax',
               'prefix' : 'demo',
               'max_seconds' : 250500,
               'tstress': True}

newcalc = qecalc.QECalcIn(control= newcontrol)

folder = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Learning\test'
file = 'test.in'
for key in newcalc.system:
    print (key)
    #print("key: " + key + "Val " + str(newcalc.namelistdict['CONTROL'][key]) + " type " + str(type(newcalc.control[key])))

