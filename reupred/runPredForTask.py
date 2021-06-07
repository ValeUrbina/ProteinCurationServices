import os
name= sys.argv[1] 
wdir= sys.argv[2]  
with open(wdir +  name, 'r') as f:
    templatesList = f.read().split("\n")
for pdbid in templatesList:
    chain = pdbid[4]
    id = pdbid[:4]
    os.system('./run.sh '+id+' '+wdir +' '+chain+' chain False False')
