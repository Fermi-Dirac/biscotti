#This function just duplicates the LabVIEW recursive file list function. Later i'll make it better.
def RecursiveFileList(filepath = "", matchpattern='*') :
    allfiles = [ y for x in os.walk(filepath) for y in glob(os.path.join(x[0],matchpattern)) ]
    return allfiles
