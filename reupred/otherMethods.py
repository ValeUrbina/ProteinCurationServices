import os
import shutil,re
from Bio.PDB import *
from alignmentManager import *
import subprocess, io
import threading
import gzip
import urllib
from threading import Thread
from threading import Thread
import os,Bio.SeqIO, Bio.AlignIO,Bio.PDB.DSSP
import shutil
from Bio.SubsMat import MatrixInfo
import seaborn as sns
import logging
import time
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Process, Manager
import numpy as np
import matplotlib
import matplotlib.pyplot
import matplotlib.pyplot as plt
import numpy as np
plt.switch_backend('agg')
matplotlib.use('Agg')

matrix = MatrixInfo.blosum62
re_info = re.compile('Aligned length=\s*(?P<len>[0-9]+), RMSD=\s*(?P<rmsd>[0-9\.]+), Seq_ID=.+=\s*(?P<id>[0-9\.]+)')
re_score = re.compile('TM-score=\s*(?P<score>[0-9\.]+)\s+\(if normalized by length of Chain_2\)')
ResiduesCodes = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER': '*',
                 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'XAA': 'X', 'UNK': 'X'}


def charge_align(sc, target, templatesList, PDBdir, ALIdir, SRULdir, DSSPdirectory,Dsspexe,TMalignexe):
    """
        This function calculates and loads the alignment objects, containing the results of aligning the target against SRUL
              :param sc: class description (III/IV/V)
              :param target: name of the target
              :param templatesList: list of SRUL templates
              :param PDBdir: path to the pdb of the target
              :param ALIdir: path where the alignment files will be stored
              :param SRULdir: path to the SRUL library
              :param DSSPdirectory: path where the secondary structure files will be stored
              :type sc: string
              :type target: string
              :type templatesList: list of strings
              :type PDBdir: string
              :type ALIdir: string
              :type SRULdir: string
              :type DSSPdirectory: string
              :return templatesTMdata: good alignment objects
              :type templatesTMdata: list of alignment objects that cover the thresholds
              :return: nali: number of alignments inside templatesTMdata
              :rtype: nali: integer
              :return: avgTM: average Tmscore value
              :rtype: avgTM: float
    """
    try:
        templatesTMdata,templatesTMdataName = [],[]
        avgTM, nali, avg = 0, 0, 0

        if target:

            alignmentProtvsUnits = ProteinAlignmentSet(target, DSSPdirectory, PDBdir,Dsspexe,TMalignexe)
            #alignmentProtvsUnits.calculate_dssp(DSSPdirectory, PDBdir)

            t = Thread(target=alignmentProtvsUnits.load_alignments, args=(templatesList, PDBdir, ALIdir, SRULdir, sc,TMalignexe))
            t.start()
            while t.is_alive():
                pass
            templatesTMdata = alignmentProtvsUnits.get_alignments()
            templatesTMdataName = alignmentProtvsUnits.get_alignments_name()
            value = alignmentProtvsUnits.get_min_rmsd()
            valueTM, avgTM = alignmentProtvsUnits.calc_max_tmscore()
            nali = len(templatesTMdata)
            return value, templatesTMdata, nali, avgTM, templatesTMdataName
        else:
            return [], templatesTMdata, nali, avgTM, templatesTMdataName
    except Exception, e:
        logging.error("Problem occurred while loading or calculating the alignments agains SRUL "+ str(e))


def parse_tmalign(output):
    """
    Parse TM-Align output
    'score' ... TM-Align score (float)
    'rmsd'  ... RMSD calculated by TM-Align (float)
    'len'   ... length (in residues) of the structure alignment (float)
    'id'    ... sequence identity based on the structure alignment (float)
    """
    if output.count('\n') < 7:
        logging.error('no TM-Align result')
    r = {}
    try:
        r = re_info.search(output).groupdict()
        r.update( re_score.search(output).groupdict())
        r = { k : float(v) for k, v in r.items() }
    except AttributeError, why:
        logging.error('Could not find info in TMAlign output')
    return r


def get_sequences_tmalign(output):
    sequences = []
    lines = output.split("\n")
    for ln in lines:
        if re.search('^[A-Z-]+$', ln):
            sequences.append(ln)
    return sequences


def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]


def score_pairwise(seq1, seq2, matrixtemp):
    # counts similar pairs
    score = 0
    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        if '-' not in pair:
            scoreal = score_match(pair, matrixtemp)
            if scoreal > 0:
                score += 1
    return float(score)/min([len(seq1), len(seq2)])

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def draw_matrix(title, alignment_matrix):
    # DataFrame
    sorted_labels = alignment_matrix.keys()
    sorted_labels.sort(key=natural_keys)
    scores = [[score_pairwise(alignment_matrix[unit1][unit2]['sequences'][0], alignment_matrix[unit1][unit2]['sequences'][1], matrix) for unit2 in sorted_labels] for unit1 in sorted_labels]
    df = pd.DataFrame(data=scores, columns=sorted_labels)
    sns.set(style="white")
    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(11, 9))
    # Generate a mask for the upper triangle
    mask = np.zeros_like(scores, dtype=np.bool)
    mask[np.tril_indices_from(mask, k=-1)] = True
    myannotation = np.array([[str(round(scores[rindex][cindex], 2)) if mask[rindex][cindex] == False else '' for cindex in range(len(mask[rindex]))] for rindex in range(len(mask))])
    sns.set(font_scale=1)
    # Draw the heatmap with the mask and correct aspect ratio
    g = sns.heatmap(df, annot=myannotation, fmt='', center=0.5, cmap='Reds', cbar=False,
                    square=True, linewidths=.5, cbar_kws={"shrink": .5, "orientation": "horizontal"},
                    yticklabels=sorted_labels, xticklabels=sorted_labels)
    if re.search("\.", title):
        plt.savefig(title[:re.search(".", title).start()-1], transparent=False, bbox_inches='tight', pad_inches=0)
    else:
        plt.savefig(title, transparent=False, bbox_inches='tight', pad_inches=0)
    plt.close()


def calc_alignments_matrix(Listunit, destination, target,UnitsDir,TMalignexe):
    """
        This function executes TMalign of all vs all the units of the region and calls draw_matrix procedure that creates the matrix
              :param Listunit : List of the units to evaluate
              :param destination : path where the TMalign results are going to be stored and where the units fragments are stored
              :param target: name of the target, used to name the matrix
              :type  Listunit : list of strings
              :type destination : string
              :type target : string
    """
    try:
        matrix,matrixM = {},{}
        report = open(destination + target + '_PairwiseReport', 'w')
        cytoreport = open(destination + target + '_NetworkReport', 'w')
        cytoreport.write('unitA unitB sequence_A TMscore \n')
        report.write('unitA unitB sequence_A sequence_B TMscore SeqId rmsd \n')
        for unit in Listunit:
            unit_name = unit[re.search('^[A-Za-z0-9]*_._', unit).end():re.search('_[A-Za-z0-9]*$', unit).start()]
            matrixM[unit_name] = {}

            matrix[unit]  = {}
            for otherunit in Listunit:
                tmOutput = subprocess.check_output(TMalignexe +' ' + UnitsDir+unit + '.pdb ' + UnitsDir + otherunit + '.pdb -a  T ', shell=True)


                other_unit_name = otherunit[
                                  re.search('^[A-Za-z0-9]*_._', otherunit).end():re.search('_[A-Za-z0-9]*$', otherunit).start()]
                TM_align_data = parse_tmalign(tmOutput)
                TM_align_data['sequences'] = get_sequences_tmalign(tmOutput)

                matrixM[unit_name][other_unit_name] = TM_align_data
                matrix[unit][otherunit] = TM_align_data
                report.write(unit + ' ' + otherunit + ' ' + matrix[unit][otherunit]['sequences'][0] + ' ' + matrix[unit][otherunit]['sequences'][1] + ' ' + str(matrix[unit][otherunit]['score']) + ' ' + str(matrix[unit][otherunit]['id']) + ' ' + str(matrix[unit][otherunit]['rmsd']) + '\n')
                cytoreport.write(unit + ' ' + otherunit + ' ' + str(matrix[unit][otherunit]['score'])+'\n')
        draw_matrix(destination + target + '_matrix', matrixM)
        report.close()
        cytoreport.close()
    except Exception, e:
        logging.error("error while creating PairwiseReport for the creation of the matrix" + str(e))


def aligned_dssp(filename, dic,dsspexe):
    """
        This function executes TMalign of all vs all the units of the region and calls draw_matrix procedure that creates the matrix
              :param filename : name of the unit fragment, including the complete path
              :type  filename : strings
              :return secstr : string
              :type secstr : string
    """
    try:
        secstr = ''
        p = PDBParser()
        structure = p.get_structure('unit', filename + '.pdb')
        model = structure[0]
        dssp = DSSP(model, filename + '.pdb',dssp=dsspexe)
        for dsspvals in list(dssp.keys()):
            secstr = secstr + dssp[dsspvals][2]
        if dic ==True:
            return dssp
        else:
            return secstr
    except Exception, e:
        logging.error("Error while trying to obtain the dssp of the aligned units | "+ filename + ".pdb" +  str(e))
        return ''


def create_aux_unit_mus(min_, max_, origen, destination, target, tempName, chain):
    """
        This function that creates the units fragments
              :param  : name of the unit fragment, including the complete path
              :param  : min_ : start pdb index for the fragment unit
              :param  : max_ : end pdb index for the fragment unit
              :param  : origen : path to the pdb from the fragment will be created
              :param  : destination : path where the new unit is going to be stored
              :param  : target : name of the protein used to create the fragment
              :param  : tempName : name for the new unit fragment
              :param  : chain : chain of the protein
              :type  : min_ : string
              :type  : max_ : string
              :type  : origen : string
              :type  : destination : string
              :type  : target : string
              :type  : tempName : string
              :type  : chain : string
              :return  : the numbers of residues in the unit fragment
              :type  : int
    """
    try:

        fpdbUnit = open(destination + tempName + ".pdb", 'w')
        if os.path.isfile(origen + target + ".pdb"):
            pdbFile = open(origen + target + ".pdb", 'r')
        else:
            pdbFile = open(origen + target, 'r')
        parser = PDBParser()
        k, countAux = 1, 0

        structureAux = parser.get_structure('temp', pdbFile)
        for res in structureAux[0][chain].get_residues():
            if res.get_id()[0] == ' ':
                for a in res:
                    if len(a.get_name()) == 4:
                        spc = ''
                    elif len(a.get_name()) == 1:
                        spc = '   '
                    elif len(a.get_name()) == 2:
                        spc = '  '
                    else:
                        spc = ' '
                    spp = '     '
                    if int(res.get_id()[1]) in range(int(min_), int(max_)+1):
                        if k >= 10000:
                                fpdbUnit.write('ATOM  '+'%5s'%str(k)+'  '+a.get_name() + spc + '%3s'%res.get_resname() + '  ' + target[4:5] + '%4s'%str(res.get_id()[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0])+'%8s'%str('%.3f'%a.get_coord()[1])+'%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor())+spp + a.get_id()[0] + '\n')
                        else:
                                fpdbUnit.write('ATOM   '+'%4s'%str(k)+'  '+a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id()[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0])+'%8s'%str('%.3f'%a.get_coord()[1])+'%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        k += 1
                        if a.get_name() == 'C':
                            countAux += 1
        fpdbUnit.close()
        return countAux
    except Exception, e:
        logging.error("Error while creating  the unit fragments" +  str(e))

def evaluating_temp_mustang(musfile):

    if os.path.isfile(musfile):
        alignment = Bio.AlignIO.read(open(musfile), "fasta")
        insertionlist = {}
        perfectAligned = 0
        alilen = len(alignment[0])

        for i in range(0, alilen):
            count = 0
            listid = []
            alignment[pair].seq = []
            for pair in range(0, len(alignment)):
                if alignment[pair].seq[i] == '-':
                    count += 1
                else:
                    listid.append(pair)
            if float(float(count) / float(len(alignment))) >= 0.65:
                for id in listid:
                    if id in insertionlist:
                        insertionlist[id] = insertionlist[id] + ',' + str(i)
                    else:
                        insertionlist[id] = str(i)
            elif float(float(len(alignment) - count) / float(len(alignment))) > 0.50:
                perfectAligned += 1
        if float(float(perfectAligned) / float(len(alignment[0]))) <= float(0.22):
            # os.remove(dbfilename)
            logging.debug("Low quality result" + str(float(float(perfectAligned) / float(len(alignment[0])))))

            return False
    return True

def calcinsertions(dbfilename, musfile, UNITSdir):
    """
        This function obtains the alignments together with the corresponding pdb indexed of the insertions and calls getinsertionpdbpositions to write them
              :param  : dbfilename : name of the db file in where the insertions will be added
              :param  : musfile : name of mustang resulting of the units alignments (.afasta) file
              :param  : UNITSdir : path to the unit fragments  used for mustang
              :type  : dbfilename : string
              :type  : musfile : string
              :type  : UNITSdir : string
    """
    try:
        logging.debug("adding to db file the insertions " + dbfilename)
        InsertionList=[]
        name,alignment=[],[]
        if os.path.isfile(musfile):
            k=0
            p=0
            flag= False


            lineList = [line.rstrip('\n') for line in open(musfile)]
            lineList = list(filter(None, lineList))  # fastest

            for line in lineList:

                if not flag :
                    name.append(line[1:])
                    flag =True
                else:
                    alignment.append(line)
                    flag=False


            #alignment = Bio.AlignIO.read(open(musfile), "fasta")
            insertionlist = {}
            perfectAligned = 0
            alilen = len(alignment[0])

            for i in range(0,alilen): #position in line
                count=0
                listid=[]

                for pair in range(0,len(alignment)): #lines

                    if alignment[pair][i] == '-':
                        count += 1
                    else:
                        listid.append(pair) #id linea

                if float(float(count)/float(len(alignment))) >= 0.65:
                    for id in listid:
                        if id in insertionlist:
                            insertionlist[id] = insertionlist[id] + ',' + str(i)
                        else:
                            insertionlist[id] = str(i)
                elif float(float(len(alignment)-count)/float(len(alignment))) > 0.50:
                    perfectAligned += 1

            if float(float(perfectAligned) / float(len(alignment[0]))) <= float(0.22):
                #os.remove(dbfilename)
                logging.INFO("Low quality result"+str(float(float(perfectAligned) / float(len(alignment[0])))))

                return False,True,InsertionList
            else:

                subseqList = {}

                for pair in insertionlist: #para cada linea
                    values = insertionlist[pair].split(',')
                    inirange = (name[pair]).split('_')[1]
                    insstart = values[0]
                    ident = name[pair]
                    insend = 99999
                    for j in range(1,len(values)):
                        if int(values[j]) - int(values[j-1]) == 1:
                            insend = values[j]
                        else:
                            if insend != 99999 and (int(insend)-int(insstart)) >= 2:

                                if ident not in subseqList:

                                    subseqList[ident] = str(((alignment[pair])[int(insstart):int(insend)+1]).translate(None,'-'))
                                else:

                                    subseqList[ident] = subseqList[ident] + ',' + str(((alignment[pair])[int(insstart):int(insend)+1]).translate(None,'-'))
                            insstart = values[j]

                    if pair and (int(insend) - int(insstart)) >= 2:
                        if ident not in subseqList:
                            subseqList[ident] = str(((alignment[pair])[int(insstart):int(insend) + 1]).translate(None,'-'))
                        else:
                            subseqList[ident] = subseqList[ident] + ',' + str(
                                ((alignment[pair])[int(insstart):int(insend) + 1]).translate(None,'-'))
                    insend = 0

                InsertionList=getinsertionpdbpositions(subseqList, dbfilename, UNITSdir)
        return True,False,InsertionList
    except Exception, e:
        logging.error("A problem occurred while calculating the insertions" + str(e))

def getinsertionpdbpositions(subseqList, dbfilename, UNITSdir):
    """
        This function obtains the insertions and writes them in the db file
              :param  : subseqList : list of units to evaluate
              :param  : dbfilename : name of the db file in where the insertions will be added
              :param  : UNITSdir : path to the unit fragments  used for mustang
              :type  : subseqList : list of unit names
              :type  : dbfilename : string
              :type  : UNITSdir : string
    """
    try:
        InsertionList=[]
        if subseqList != {}:
            logging.debug("adding to file getinsertionpdbpositions " + dbfilename)
            filedb = open(dbfilename, 'a')
            for id in subseqList:
                seq = ''
                pdbFile = UNITSdir+id
                p = PDBParser()
                structure = p.get_structure(pdbFile, pdbFile)
                for model in structure:
                    for chain in model:
                        seqpos = list()
                        for residue in chain:
                            if is_aa(residue.get_resname(), standard=True):
                                seq=seq+ResiduesCodes[residue.get_resname()]
                                seqpos.append(residue.get_id()[1])
                            else:
                                seq = seq + ResiduesCodes[residue.get_resname()]
                                seqpos.append(-1)
                insertions = subseqList[id].split(',')
                for subseq in insertions:
                    pos = seq.find(subseq)
                    if int(seqpos[pos + len(subseq) - 1]) - int(seqpos[pos]) >= 3:
                        InsertionList.append('INS\t' + str(seqpos[pos]) + ' ' + str(seqpos[pos+len(subseq)-1]) + '\n')
            return InsertionList
    except Exception, e:
        logging.error("An error occurred while writting the insertions" + str(e))


def deleteunitsofregion(musUnits,fastpListAA):
    for line in musUnits:
        fastpListAA.remove(line)
    return fastpListAA

def evaluatemustangResult(fastpListAA, tempdir,OriginalTargetname, chain, Mustangexe):
    valResp, changed =False,False
    prevLast=-999
    musUnits=[]
    musfile = open(tempdir+'tryMus' + '.mus', 'w')
    musfile.write('>' + tempdir + '\n')

    musUnitsfinal=[]
    for unittext in fastpListAA:


            unit=unittext.split('\t')

            if  prevLast==-999 or int(unit[0]) - int(prevLast) ==1 :
                musfile.write('+' + 'Repeat' + '_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + '.pdb\n')
                create_aux_unit_mus(unit[0], unit[1], tempdir,tempdir,
                                    OriginalTargetname + '.pdb',
                                    'Repeat' + '_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' , chain)
                musUnits.append(unittext)
            else:


                musfile.close()
                if os.path.isfile(tempdir+'tryMus' + '.mus'):
                    command = Mustangexe+' -f ' + tempdir+'tryMus' + '.mus' + ' -F fasta  -r ON  -o ' + tempdir+'tryMus' + '.mustang'
                    logging.debug(command)
                    os.system(command)
                    if os.path.isfile(tempdir+'tryMus' + '.mustang.afasta') and os.path.isfile(tempdir+'tryMus' + '.mustang.pdb'):
                        valid = True#evaluating_temp_mustang(tempdir+'tryMus' + '.mustang.afasta')
                        os.system('rm '+ tempdir+'tryMus' + '.mustang.afasta')
                        os.system('rm ' + tempdir + 'tryMus' + '.mustang.pdb')
                        if not valid:
                            musUnitsfinal = musUnitsfinal + musUnits
                            musUnits.append(unittext)
                            changed = True
                        valResp = valResp or valid
                        musUnits=[]

                    else:
                        valResp = valResp or False
                        musUnitsfinal = musUnitsfinal + musUnits
                        changed = True
                        musUnits = []
                        musUnits.append(unittext)
                musfile = open(tempdir + 'tryMus' + '.mus', 'w')
                musfile.write('>' + tempdir + '\n')
            prevLast=unit[1]
    musfile.close()
    if os.path.isfile(tempdir + 'tryMus' + '.mus'):

        os.system(
            Mustangexe+' -f ' + tempdir + 'tryMus' + '.mus' + ' -F fasta  -r ON  -o ' + tempdir + 'tryMus' + '.mustang')
        if os.path.isfile(tempdir + 'tryMus' + '.mustang.afasta') and os.path.isfile(
                                tempdir + 'tryMus' + '.mustang.pdb'):
            valid = True#evaluating_temp_mustang(tempdir + 'tryMus' + '.mustang.afasta')
            if not valid:
                musUnitsfinal = musUnitsfinal + musUnits
                changed = True
            valResp = valResp or valid
        else:
            valResp = valResp or False
            musUnitsfinal = musUnitsfinal+musUnits
            changed = True
    if changed:
        fastpListAA = deleteunitsofregion(musUnitsfinal, fastpListAA)


    return valResp, changed, fastpListAA


def create_pdb_description_file(wdirOriginal ,pdb_header,name):

    info=open(wdirOriginal+name+'.info','w')
    for keys, values in pdb_header.items():
        if keys =='source' or keys =='compound':
            info.write(keys+'\n')
            #for keys2, values2 in values:
             #   info.write(str(keys2 + '\t' + str(values2) + '\n')
        else:
            info.write(keys +'\t'+ str(values)+'\n')
    info.close()


def intersect(a, b):
    result=[]
    valT=False
    for i in b:
        if isinstance(i,list):
            result.append(intersect(a,i))
        else:
            if i in a:
                 result.append(i)
                 valT=True
    return valT


def reevaluate_regions(folder, target, OriginalTargetname, regID, Raphaelperiod, chain,Mustangexe,TMalignexe,dsspexe,wdirOriginal):
    """
        This function sets all the results in the right folder structure
              :param  : folder : path to the chain folder in where the region folder will be created
              :param  : target : name of the target pdb
              :param  : OriginalTargetname : name of the original input file
              :param  : regID : number ot the region that is been evaluated
              :param  : Raphaelperiod : Raphael's value of the target
              :param  : chain :  chain of the target pdb
              :type  : folder : string
              :type  : target : string
              :type  : OriginalTargetname : string
              :type  : regID : number
              :type  : Raphaelperiod : float
              :type  : chain : string
    """
    try :

        regionList, unitList, insList,regions, inslistStart, inslistEnd, selectedTemplateList, types, folds = [], [], [], [], [], [], [], [], []
        idreg=0
        isValid=False
        regionsListdata = {}
       # print 'mkdir ' + folder + 'Region' + str(idreg + 1)
        for countregions in range(0, regID + 1):

            logging.debug('mkdir ' + folder + 'Region' + str(idreg + 1))

            if os.path.isfile(folder + 'region' + str(idreg)):
                with open(folder + 'region' + str(idreg), 'r') as f:
                    info = f.read().split("\n")
                counterunit=0
                for line in info:
                    row = line.split("\t")
                    if row[0] == 'MASTER':
                        selectedTemplateList.append(row)
                    if row[0] == "REG":

                        infreg = row[1].split(" ")
                        regionsListdata['Region' + str(idreg + 1)] = range(int(infreg[0]), int(infreg[1]) + 1)
                        regions.append([range(int(infreg[0]), int(infreg[1]) + 1)])
                        if infreg[3] == "X":
                            regionList.append([int(infreg[0]), int(infreg[1]), infreg[2], 0])
                        else:
                            regionList.append([int(infreg[0]), int(infreg[1]), infreg[2], int(infreg[3])])
                            if len(infreg)>6:
                                if infreg[6]=='-' and  infreg[7] =='-':
                                    types.append([int(infreg[0]), int(infreg[1]) , 0, 0])
                                elif  infreg[7] =='-':
                                    types.append([int(infreg[0]), int(infreg[1]), infreg[6], 0 ])
                                elif infreg[6]=='-' and  infreg[7] =='-':
                                    types.append([int(infreg[0]), int(infreg[1]) , 0, int(infreg[7])])
                                else:
                                    types.append([int(infreg[0]), int(infreg[1]), infreg[6], int(infreg[7])])
                            if  len(infreg) > 4:
                                if infreg[4]== '-':
                                    if  infreg[5]== '-':
                                        folds.append([int(infreg[0]), int(infreg[1]), 0, 0])
                                    else:
                                        folds.append([int(infreg[0]), int(infreg[1]), 0,  int(infreg[5])])
                                else:
                                    if  infreg[5]== '-':
                                        folds.append([int(infreg[0]), int(infreg[1]),int(infreg[4]), 0])
                                    else:
                                        folds.append([int(infreg[0]), int(infreg[1]), int(infreg[4]),  int(infreg[5])])
                    if row[0] == "TYPE":

                        inftype = row[1].split(" ")
                        types.append([int(inftype[0]), int(inftype[1]), inftype[2], int(inftype[3])])
                    if row[0] == "FOLD":

                        inffold= row[1].split(" ")
                        folds.append([int(inffold[0]), int(inffold[1]), inffold[2], int(inffold[3])])
                    elif row[0] == "UNIT":
                        counterunit=counterunit+1
                        infunit = row[1].split(" ")
                        unitList.append([int(infunit[0]), int(infunit[1])])
                    elif row[0] == "INS":
                        infins = row[1].split(" ")
                        insList.append([int(infins[0]), int(infins[1])])
                        inslistEnd.append(int(infins[1]))
                        inslistStart.append(int(infins[0]))

                if counterunit>2:
                    os.system('mkdir ' + folder + 'Region' + str(idreg + 1))

                else:
                    idreg=idreg+1
            if os.path.isfile(folder + 'region' + str(countregions)):
                os.remove(folder + 'region' + str(countregions))
            idreg = idreg+1



        regioncount = 1

        if insList!=[]:
            logging.warning("Some insertions had to be excluded from the units fragments to be able to align them")

        valuesRegErase=[]
        logging.debug('these are the valid regions'+ str(regionList))
        #print "pass 0"
        for reg in regionList:
            logging.debug("evaluating region" + str(reg))
            Listunit, dsspList = [], {}
            valuesinregion = range(int(reg[0]), int(reg[1]))
            musfile = open(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mus', 'w')
            os.system('mkdir '+folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units')
            musfile.write('>' + folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units/\n')
           # print "pass 1"
            for unit in unitList:
                if int(unit[0]) in valuesinregion:
                    musfile.write('+' + OriginalTargetname + '_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + str(regioncount) + '.pdb\n')
                    create_aux_unit_mus(unit[0], unit[1], folder, folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units/',  OriginalTargetname + '.pdb', OriginalTargetname + '_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + str(regioncount),chain)
                    Listunit.append(OriginalTargetname + '_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + str(regioncount))
                    dsspList[OriginalTargetname+'_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + str(regioncount) + '.pdb'] = (aligned_dssp(folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units/' + OriginalTargetname + '_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + str(regioncount),False,dsspexe))
            musfile.close()
           # print "pass 2", folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mus'
            if os.path.isfile(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mus'):
                #mustOutput=''
                os.system(Mustangexe+' -f ' + folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mus' + ' -F fasta  -r ON  -o ' + folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang')
                #print "inside", mustOutput
                #mustOutput = subprocess.check_output(Mustangexe+' -f ' + folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mus' + ' -F fasta  -r ON  -o ' + folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang', shell=True)
                if os.path.isfile(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang.pdb'):
                    files = [wdirOriginal + "header",
                             folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(
                                 regioncount) + '.mustang.pdb']
                    concat = ''.join([open(f).read() for f in files])
                    auxfile = open(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang.pdb', 'w')
                    auxfile.write(concat)
                    auxfile.close()
                    calc_alignments_matrix(Listunit, folder + 'Region' + str(regioncount) + '/', OriginalTargetname + 'region' + str(regioncount), folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units/',TMalignexe)
            if os.path.isfile(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang.pdb') and os.path.isfile(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang.afasta'):
                with open(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang.afasta') as f:
                    mustanginfo = f.read().split("\n")
                alignmus = {}
                fdssp = open(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.aligned_dssp', 'w')

                seq = ''
                print "mustang info" , mustanginfo
                for line in mustanginfo:
                    if line != '':
                        if line[0] == '>':
                            name = line[1:]
                        else:
                            seq += line
                    else:
                        alignmus[name] = seq

                        dssp, counter = '', 0
                        for let in alignmus[name]:
                            if let == '-':
                                dssp += '-'
                            elif counter < len(dsspList[name]):
                                dssp += dsspList[name][counter]

                                counter += 1
                            else:
                                dssp += '-'

                        fdssp.write('>' + name + '\n' + dssp + '\n')
                        seq = ''
                fdssp.close()
               # print "alignus",alignmus
              #  print sorted(alignmus.items(), key=itemgetter(0))
                # for key in alignmus:
                #
                #     dssp, counter = '', 0
                #     for let in alignmus[key]:
                #         if let == '-':
                #             dssp += '-'
                #         elif counter < len(dsspList[key]):
                #             dssp += dsspList[key][counter]
                #
                #             counter += 1
                #         else:
                #             dssp += '-'
                #     fdssp.write('>'+key + '\n' + dssp + '\n')


                isValidparc,eraseReg,InsertionList=calcinsertions(folder + OriginalTargetname + ".db", folder + 'Region' + str( regioncount) + '/' + OriginalTargetname + 'region' + str(  regioncount) + '.mustang.afasta', folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units/')
               # print "CalcINSertions", isValidparc,eraseReg,InsertionList
                if eraseReg:

                    os.system('rm -r '+folder +'Region' + str( regioncount))
                    valuesRegErase = valuesRegErase + regionsListdata['Region' + str( regioncount)]
                else:
                    shutil.make_archive(
                        folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(
                            regioncount) + '_units', 'zip',
                        folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(
                            regioncount) + '_units/')
                    shutil.rmtree(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(
                         regioncount) + '_units/')

                isValid= isValidparc or isValid
            else:
                os.system('rm -r ' + folder + 'Region' + str(regioncount))
                valuesRegErase = valuesRegErase + regionsListdata['Region' + str(regioncount)]

                logging.warning("It was not possible to calculate the units alignment")
                isValid = False or isValid
                eraseReg=False
                InsertionList=[]

            regioncount += 1

        InsertionList=[]

        if isValid:
            logging.info("adding to file reevaluate_regions " + folder + OriginalTargetname + ".db")
            outf = open(folder + OriginalTargetname + ".db", "w")

            outf.write("SOURCE\tReUPred\n")
            outf.write("DATE\t" + (time.strftime("%d/%m/%Y")) + "\n")
            outf.write("PDB\t" + OriginalTargetname[:len(OriginalTargetname)-2] + "\n")
            outf.write("CHAIN\t" + chain + "\n")
            if Raphaelperiod and chain in Raphaelperiod:
                outf.write("RAPHAEL\t" + str(Raphaelperiod[target[:5]]) + "\n")
            else:
                outf.write("RAPHAEL\t*\n")
            outf.write("REPEATED\tYes\n")
            for r in regionList:

                if not intersect(valuesRegErase,range(int(r[0]),int(r[1]))):

                    outf.write("REG\t" + str(r[0]) + " " + str(r[1]) + " " + str(r[2]) + " " + str(r[3]) + "\n")

            for selectedTemplate in  selectedTemplateList:
                valuesmaster = selectedTemplate[1].split(' ')
                if not intersect(valuesRegErase, range(int(valuesmaster[0]),int(valuesmaster[1]))):
                    outf.write("MASTER\t" + str(valuesmaster[0]) + " " + str(valuesmaster[1]) + ' ' + str(valuesmaster[2]) + " " + str(valuesmaster[3]) + ' ' + str(valuesmaster[4]) +"\n")

            for unittype in types:
                if not intersect(valuesRegErase, range(int(unittype[0]),int(unittype[1]))):
                    outf.write("TYPE\t" + str(unittype[0]) + " " + str(unittype[1]) + ' ' + str(unittype[2]) + " " + str(unittype[3]) + "\n")

            for folds in  folds:
                if not intersect(valuesRegErase,range(int(folds[0]),int(folds[1]))):
                    outf.write("FOLD\t" + str(folds[0]) + " " + str(folds[1]) + ' ' + str(folds[2]) + " " + str(folds[3]) +"\n")

            for uni in unitList:
                if not intersect(valuesRegErase, range(int(uni[0]),int(uni[1]))):
                    outf.write("UNIT\t" + str(uni[0]) + " " + str(uni[1]) + "\n")

            for r in insList:
                if r[1] - r[0] > 3:
                    if not intersect(valuesRegErase, range(int(r[0]),int(r[1]))):
                        outf.write("INS\t" + str(r[0]) + ' ' + str(r[1]) + "\n")
            #print InsertionList
            for r in InsertionList:
                if r!=[] or r!='':
                    outf.write(r)
            outf.close()

        if os.path.isfile(folder + OriginalTargetname + ".db"):
            os.system('cp '+folder + OriginalTargetname + ".db " + folder + OriginalTargetname + ".reupred" )
        if os.path.isfile(folder + 'temp' + chain + '.db'):
            os.remove(folder + 'temp' + chain + '.db')
        if os.path.isfile(folder + 'temp_' + chain + '.mapping'):
            os.system('cp ' + folder + 'temp_' + chain + '.mapping ' + folder + OriginalTargetname + ".mapping")
            os.system('rm ' + folder + 'temp_' + chain + '.mapping')
        if os.path.isfile(folder + 'temp_' + chain + '.pdb'):
            os.system('cp ' + folder + 'temp_' + chain + '.pdb ' + folder + OriginalTargetname + ".pdb")
            os.system('rm '+folder + 'temp_' + chain + '.pdb')
        files = [wdirOriginal + "header", folder + OriginalTargetname + ".pdb"]
        concat = ''.join([open(f).read() for f in files])
        auxfile = open(folder + OriginalTargetname + ".pdb", 'w')
        auxfile.write(concat)
        auxfile.close()
        return isValid
    except Exception, e:
        logging.error("An error occurred while setting all the resulting data in place" + str(e))


def reevaluate_regions_validation(folder, OriginalTargetname, chain,Mustangexe,TMalignexe,dsspexe, wdirOriginal):
    """
        This function sets all the results in the right folder structure
              :param  : folder : path to the chain folder in where the region folder will be created
              :param  : target : name of the target pdb
              :param  : OriginalTargetname : name of the original input file
              :param  : regID : number ot the region that is been evaluated
              :param  : Raphaelperiod : Raphael's value of the target
              :param  : chain :  chain of the target pdb
              :type  : folder : string
              :type  : target : string
              :type  : OriginalTargetname : string
              :type  : regID : number
              :type  : Raphaelperiod : float
              :type  : chain : string
    """
    try :
        regionList, unitList, insList,regions, inslistStart, inslistEnd, selectedTemplateList, types, folds = [], [], [], [], [], [], [], [], []
        countregions=1
        if os.path.isfile(folder + OriginalTargetname + '.db'):
            with open(folder + OriginalTargetname + '.db', 'r') as f:
                info = f.read().split("\n")
            try:
                for line in info:
                    row = line.split("\t")
                    if row[0] == 'MASTER':
                        selectedTemplateList.append(row)
                    elif row[0] == "REG":
                        os.system('mkdir ' + folder + 'Region' + str(countregions  ))
                        logging.debug('mkdir ' + folder + 'Region' + str(countregions ))
                        countregions=countregions+1
                        infreg = row[1].split(" ")
                        regions.append([range(int(infreg[0]), int(infreg[1]) + 1)])
                        if infreg[3] == "X":
                            regionList.append([int(infreg[0]), int(infreg[1]), infreg[2], 0])

                        else:
                            regionList.append([int(infreg[0]), int(infreg[1]), infreg[2], int(infreg[3])])

                    elif row[0] == "UNIT":
                        infunit = row[1].split(" ")
                        unitList.append([int(infunit[0]), int(infunit[1])])
                    elif row[0] == "INS":
                        infins = row[1].split(" ")
                        insList.append([int(infins[0]), int(infins[1])])
                        inslistEnd.append(int(infins[1]))
                        inslistStart.append(int(infins[0]))
            except Exception,e :
                logging.error("there are errors in the db file format"+folder + OriginalTargetname + '.db')
                logging.warning("there are errors in the db file format"+folder + OriginalTargetname + '.db')

        regioncount = 1
        for reg in regionList:
            Listunit, dsspList = [], {}
            valuesinregion = range(int(reg[0]), int(reg[1]))
            musfile = open(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mus', 'w')
            os.system('mkdir '+folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units')
            musfile.write('>' + folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units/\n')
            for unit in unitList:
                if int(unit[0]) in valuesinregion:
                    musfile.write('+' + OriginalTargetname + '_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + str(regioncount) + '.pdb\n')
                    create_aux_unit_mus(unit[0], unit[1], folder, folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units/',  OriginalTargetname + '.pdb', OriginalTargetname + '_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + str(regioncount),chain)
                    Listunit.append(OriginalTargetname + '_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + str(regioncount))
                    dsspList[OriginalTargetname+'_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + str(regioncount) + '.pdb'] = (aligned_dssp(folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units/' + OriginalTargetname + '_' + str(unit[0]) + '_' + str(unit[1]) + '_reg' + str(regioncount),False,dsspexe))
            musfile.close()
            if os.path.isfile(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mus'):
                os.system(Mustangexe+' -f ' + folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mus' + ' -F fasta  -r ON  -o ' + folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang')
                if os.path.isfile(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang.pdb'):
                    files = [wdirOriginal + "header",
                             folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(
                                 regioncount) + '.mustang.pdb']
                    concat = ''.join([open(f).read() for f in files])
                    auxfile = open(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang.pdb', 'w')
                    auxfile.write(concat)
                    auxfile.close()

                logging.debug('mustang -f ' + folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mus' + ' -F fasta  -r ON  -o ' + folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang')
            calc_alignments_matrix(Listunit, folder + 'Region' + str(regioncount) + '/', OriginalTargetname + 'region' + str(regioncount), folder + 'Region' + str(regioncount) + '/'+OriginalTargetname+'region' + str(regioncount)+'_units/',TMalignexe)
            if os.path.isfile(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang.afasta'):
                with open(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + '.mustang.afasta') as f:
                    mustanginfo = f.read().split("\n")
                alignmus = {}
                seq = ''
                for line in mustanginfo:
                    if line != '':
                        if line[0] == '>':
                            name = line[1:]
                        else:
                            seq += line
                    else:
                        alignmus[name] = seq
                        seq = ''
                fdssp = open(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(regioncount) + ' .aligned_dssp', 'w')
                for key in alignmus:
                    dssp, counter = '', 0
                    for let in alignmus[key]:
                        if let == '-':
                            dssp += '-'
                        else:
                            dssp += dsspList[key][counter]
                            #dssp += dsspList[key][counter]
                            counter += 1
                    fdssp.write('>'+key + '\n' + dssp + '\n')
                fdssp.close()
                shutil.make_archive(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(
                    regioncount) + '_units', 'zip',
                                    folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(
                                        regioncount) + '_units/')
                shutil.rmtree(folder + 'Region' + str(regioncount) + '/' + OriginalTargetname + 'region' + str(
                                        regioncount) + '_units/')
            regioncount += 1
        if os.path.isfile(folder + 'temp' + chain + '.db'):
            os.remove(folder + 'temp' + chain + '.db')
        if os.path.isfile(folder + 'temp_' + chain + '.mapping'):
            os.system('cp ' + folder + 'temp_' + chain + '.mapping ' + folder + OriginalTargetname + ".mapping")
            os.system('rm ' + folder + 'temp_' + chain + '.mapping')
        if os.path.isfile(folder + 'temp_' + chain + '.pdb'):
            os.system('cp ' + folder + 'temp_' + chain + '.pdb ' + folder + OriginalTargetname + ".pdb")
            os.system('rm '+folder + 'temp_' + chain + '.pdb')
        files = [wdirOriginal + "header", folder + OriginalTargetname + ".pdb"]
        concat = ''.join([open(f).read() for f in files])
        auxfile = open(folder + OriginalTargetname + ".pdb", 'w')
        auxfile.write(concat)
        auxfile.close()
    except Exception, e:
        logging.error("An error occurred while setting all the resulting data in place" + str(e))


def identify_order_predictor(target, PDBdir, ALIdir, SRULdir, DSSPdirectory, templatesListIII, templatesListIV, templatesListV,Dsspexe,TMalignexe):
    """
        This function identifies the order of the classes for the prediction process using reduced SRUL
              :param  : target : name of the target
              :param  : PDBdir : path to the padb file
              :param  : SRULdir : path to the SRUL library
              :param  : DSSPdirectory : path to the dssp directory
              :param  : templatesListIII : Reduced list of SRUL unit for class III
              :param  : templatesListIV : Reduced list of SRUL unit for class IV
              :param  : templatesListV : Reduced list of SRUL unit for class V
              :type  : target : string
              :type  : PDBdir : string
              :type  : SRULdir : string
              :type  : DSSPdirectory : string
              :type  : templatesListIII : list of strings
              :type  : templatesListIV : list of strings
              :type  : templatesListV : list of strings
              :return  : order : list of class index 0 for class III, 1 for class IV and 2 for class V
              :return  : existaligns : flag to identif if good alignments exist
              :type  : order : list of integer
              :type  : existaligns :boolean
    """
    try:
        existaligns = True
        eliminate=False
        templatesTMdataNameIII, templatesTMdataNameIV, templatesTMdataNameV = [], [], []
        try:
            templatesTMdataIII, templatesTMdataIV, templatesTMdataV = [], [], []
            rmsd_valIII, templatesTMdataIII, naliIII, tmIII, templatesTMdataNameIII = charge_align('III', target, templatesListIII, PDBdir, ALIdir + '/', SRULdir, DSSPdirectory, Dsspexe,TMalignexe)
            rmsd_valIV, templatesTMdataIV, naliIV, tmIV, templatesTMdataNameIV = charge_align('IV', target, templatesListIV, PDBdir, ALIdir + '/', SRULdir, DSSPdirectory,Dsspexe,TMalignexe)
            rmsd_valV, templatesTMdataV, naliV, tmV, templatesTMdataNameV = charge_align('V', target, templatesListV, PDBdir, ALIdir + '/', SRULdir, DSSPdirectory,Dsspexe,TMalignexe)

        except  Exception, e:
            logging.error("Unable to align target against SRUL"+ str(e))
        try:
            #print rmsd_valIII, rmsd_valIV, rmsd_valV
           # print tmIII,tmIV, tmV
            if rmsd_valIII == 999 or rmsd_valV == 999 or rmsd_valIV == 999:
                eliminate =True
                if  ((rmsd_valIII == 999 or rmsd_valV == 999 or rmsd_valIV == 999 ) and ( rmsd_valIII == rmsd_valIV or rmsd_valIII ==  rmsd_valV or rmsd_valV == rmsd_valIV )) or (rmsd_valIII == 0.0 or rmsd_valV == 0.0 or rmsd_valIV == 0.0):
                    order = identify_class(rmsd_valIII, rmsd_valIV, rmsd_valV, False, eliminate)
                elif abs(float(tmIII) - float(tmIV)) < float(0.04) or abs(float(tmIII) - float(tmV)) < float(0.04) or abs(float(tmV) - float(tmIV)) < float(0.04) :
                    order = identify_class(rmsd_valIII, rmsd_valIV, rmsd_valV, False, eliminate)
                else:
                    order = identify_class(float(tmIII) * naliIV, float(tmIV) * naliIV, float(tmV) * naliIV,False, eliminate)
               # print "RMSD"
            elif abs(float(tmIII) - float(tmIV)) < float(0.04) or abs(float(tmIII) - float(tmV)) < float(0.04) or abs(float( tmV) - float(tmIV)) < float(0.04):
                order = identify_class(rmsd_valIII, rmsd_valIV, rmsd_valV, False, eliminate)
                #print "rmsd no neliminate"
            else:
                order = identify_class(float(tmIII) * naliIV, float(tmIV) * naliIV, float(tmV) * naliIV, True, False)
        except Exception, e:
            logging.error("Unable to calculate the order for the prediction process"+ str(e))
        if templatesTMdataV == templatesTMdataIV and templatesTMdataIII == templatesTMdataIV and templatesTMdataIII == []:
            logging.info("There are no alignments against the SRUL")
            existaligns = False
        return order, existaligns, templatesTMdataNameIII ,templatesTMdataNameIV, templatesTMdataNameV
    except Exception, e:
        logging.error("An error occurred while trying to obtain the order of the classes for the prediction process"+ str(e))


def load_templates_list(directory, templateList, templatesTMdataNameIII ,templatesTMdataNameIV, templatesTMdataNameV):
    templatesListV, templateListV1, templatesListV2, templateListV3, templatesListV4, templateListV5 = [], [], [], [], [], []
    templatesListIII, templateListIII1, templateListIII2, templateListIII3, templateListIII4, templateListIII5, templateListIII6 = [], [], [], [], [], [], []
    templatesListIV, templateListIV1, templatesListIV2, templateListIV3, templatesListIV4, templateListIV5, templatesListIV6, templateListIV7, templatesListIV8, templateListIV9, templatesListIV10 = [], [], [], [], [], [], [], [], [], [], []
    data_directory=directory + 'Data/'
    with open(data_directory + templateList + 'V1', 'r') as f:
        templateListV1 = f.read().split("\n")
    with open(data_directory + templateList + 'V2', 'r') as f:
        templateListV2 = f.read().split("\n")
    with open(data_directory + templateList + 'V3', 'r') as f:
        templateListV3 = f.read().split("\n")
    with open(data_directory + templateList + 'V4', 'r') as f:
        templateListV4 = f.read().split("\n")
    with open(data_directory + templateList + 'V5', 'r') as f:
        templateListV5 = f.read().split("\n")
    with open(data_directory + templateList + 'IV1', 'r') as f:
        templateListIV1 = f.read().split("\n")
    with open(data_directory + 'clusterLipocalin', 'r') as f:
        LIPOCALINlike = f.read().split("\n")
    with open(data_directory + 'clusterOstac', 'r') as f:
        OSTAClike = f.read().split("\n")
    with open(data_directory + templateList + 'IV3', 'r') as f:
        templateListIV3 = f.read().split("\n")
    with open(data_directory + templateList + 'IV4', 'r') as f:
        templateListIV4 = f.read().split("\n")
    with open(data_directory + templateList + 'IV5', 'r') as f:
        templateListIV5 = f.read().split("\n")
    with open(data_directory + templateList + 'IV6', 'r') as f:
        templateListIV6 = f.read().split("\n")
    with open(data_directory + templateList + 'IV7', 'r') as f:
        templateListIV7 = f.read().split("\n")
    with open(data_directory + templateList + 'IV8', 'r') as f:
        templateListIV8 = f.read().split("\n")
    with open(data_directory + templateList + 'IV9', 'r') as f:
        templateListIV9 = f.read().split("\n")
    with open(data_directory + templateList + 'IV10', 'r') as f:
        templateListIV10 = f.read().split("\n")
    with open(data_directory + templateList + 'III1', 'r') as f:
        templateListIII1 = f.read().split("\n")
    with open(data_directory + templateList + 'III2', 'r') as f:
        templateListIII2 = f.read().split("\n")
    with open(data_directory + 'clusterTpr', 'r') as f:
        templateListTPRlike = f.read().split("\n")
    with open(data_directory + 'clusterAnkyrin', 'r') as f:
        templateListANKlike = f.read().split("\n")
    with open(data_directory +'clusterArmadillo', 'r') as f:
        templateListARMlike = f.read().split("\n")
    with open(data_directory + 'clusterTal', 'r') as f:
        templateListTALlike = f.read().split("\n")
    with open(data_directory + 'clusterPumilio', 'r') as f:
        templateListPUMlike = f.read().split("\n")
    with open(data_directory + templateList + 'III4', 'r') as f:
        templateListIII4 = f.read().split("\n")
    with open(data_directory + templateList + 'III5', 'r') as f:
        templateListIII5 = f.read().split("\n")
    with open(data_directory + templateList + 'III6', 'r') as f:
        templateListIII6 = f.read().split("\n")
    for name in templatesTMdataNameIII:
        if name in templateListIII1:
            templatesListIII= templatesListIII+ templateListIII1
        elif name in templateListIII2:
            templatesListIII= templatesListIII+ templateListIII2
        elif name in templateListTPRlike:
            templatesListIII= templatesListIII+ templateListTPRlike
        elif name in templateListANKlike:
            templatesListIII= templatesListIII+ templateListANKlike
        elif name in templateListARMlike:
            templatesListIII= templatesListIII+ templateListARMlike
        elif name in templateListTALlike:
            templatesListIII= templatesListIII+ templateListTALlike
        elif name in templateListPUMlike:
            templatesListIII= templatesListIII+ templateListPUMlike
        elif name in templateListIII4:
            templatesListIII= templatesListIII+ templateListIII4
        elif name in templateListIII5:
            templatesListIII= templatesListIII+ templateListIII5
        elif name in templateListIII6:
            templatesListIII= templatesListIII+ templateListIII6
    for name in templatesTMdataNameIV:
        if name in templateListIV1:
            templatesListIV= templatesListIV+ templateListIV1
        elif name in LIPOCALINlike:
            templatesListIV= templatesListIV+ LIPOCALINlike
        elif name in OSTAClike:
            templatesListIV= templatesListIV+ OSTAClike
        elif name in templateListIV3:
            templatesListIV = templatesListIV + templateListIV3
        elif name in templateListIV4:
            templatesListIV = templatesListIV + templateListIV4
        elif name in templateListIV5:
            templatesListIV = templatesListIV + templateListIV5
        elif name in templateListIV6:
            templatesListIV = templatesListIV + templateListIV6
        elif name in templateListIV7:
            templatesListIV = templatesListIV + templateListIV7
        elif name in templateListIV8:
            templatesListIV = templatesListIV + templateListIV8
        elif name in templateListIV9:
            templatesListIV = templatesListIV + templateListIV9
        elif name in templateListIV10:
            templatesListIV = templatesListIV + templateListIV10
    for name in templatesTMdataNameV:
        if name in templateListV1:
            templatesListV = templatesListV + templateListV1
        elif name in templateListV2:
            templatesListV = templatesListV + templateListV2
        elif name in templateListV3:
            templatesListV = templatesListV + templateListV3
        elif name in templateListV4:
            templatesListV = templatesListV + templateListV4
        elif name in templateListV5:
            templatesListV= templatesListV+ templateListV5
    templatesListIII=list(set(templatesListIII))
    templatesListIV=list(set(templatesListIV))
    templatesListV=list(set(templatesListV))
    logging.info("the predictor did " + str(len(templatesListIII)) + " alignments in class III")
    logging.info("the predictor did " + str(len(templatesListIV))+ " alignments in class IV")
    logging.info("the predictor did " + str(len(templatesListV)) + " alignments in class V")
    #exclude code
    # with open(data_directory  + "ExcludeListAllSRUL60", 'r') as f:
    #     excludedFiles = f.read().split("\n")
    # tempLIII=[]
    # print " "
    # for line in templatesListIII:
    #     if line not in excludedFiles:
    #         tempLIII.append(line)
    # templatesListIII=list(tempLIII)
    # tempLIV=[]
    # for line in templatesListIV:
    #     if line not in excludedFiles:
    #         tempLIV.append(line)
    # templatesListIV=list(tempLIV)
    # tempLV=[]
    # for line in templatesListV:
    #     if line not in excludedFiles:
    #         tempLV.append(line)
    # templatesListV=list(tempLV)

    #exclude code
    logging.info("the predictor did " + str(len(templatesListIII)) + " alignments in class III")
    logging.info("the predictor did " + str(len(templatesListIV))+ " alignments in class IV")
    logging.info("the predictor did " + str(len(templatesListV)) + " alignments in class V")
    return templatesListIII, templatesListIV, templatesListV


def found_order_predictor(target, PDBdir, ALIdir, SRULdir, DSSPdirectory, listorder, templatesTMdataIII, templatesTMdataIV, templatesTMdataV, templatesTMdataNameIII ,templatesTMdataNameIV, templatesTMdataNameV, directory, templateList,Dsspexe,TMalignexe):
    """
        This function identifies the order of the classes for the prediction process using the whole SRUL on the possible classes
              :param : target : name of the target
              :param : PDBdir : path to the padb file
              :param : SRULdir : path to the SRUL library
              :param : DSSPdirectory : path to the dssp directory
              :param : templatesListIII : Reduced list of SRUL unit for class III
              :param : templatesListIV : Reduced list of SRUL unit for class IV
              :param : templatesListV : Reduced list of SRUL unit for class V
              :param : listorder : list of classes to evaluate
              :type : target : string
              :type : PDBdir : string
              :type : SRULdir : string
              :type : DSSPdirectory : string
              :type : templatesListIII : list of strings
              :type : templatesListIV : list of strings
              :type : templatesListV : list of strings
              :type : listorder : list of integer
              :return : order : list of class index 0 for class III, 1 for class IV and 2 for class V
              :return : existaligns : flag to identif if good alignments exist
              :return : templatesTMdataIII : list of SRUL units that aligned good in the class III
              :return : templatesTMdataIV : list of SRUL units that aligned good in the class IV
              :return : templatesTMdataV : list of SRUL units that aligned good in the class V
              :type : order : list of integer
              :type : templatesTMdataIII : list of string
              :type : templatesTMdataIV : list of string
              :type : templatesTMdataV : list of string
              :type : existaligns :boolean
    """
    try:
        existaligns = True
        rmsd_valIII = 999
        rmsd_valV = 999
        rmsd_valIV = 999
        templatesListIII, templatesListIV, templatesListV = load_templates_list(directory, templateList, templatesTMdataNameIII ,templatesTMdataNameIV, templatesTMdataNameV)
        try:

            if 0 in listorder:
                rmsd_valIII, templatesTMdataIII, naliIII, tmIII, templatesTMdataNameIII = charge_align('III', target, templatesListIII, PDBdir, ALIdir + '/', SRULdir, DSSPdirectory,Dsspexe,TMalignexe)
            if 1 in listorder:

                rmsd_valIV, templatesTMdataIV, naliIV, tmIV, templatesTMdataNameIV = charge_align('IV', target, templatesListIV, PDBdir, ALIdir + '/', SRULdir, DSSPdirectory,Dsspexe,TMalignexe)
            if 2 in listorder:
                rmsd_valV, templatesTMdataV, naliV, tmV, templatesTMdataNameV = charge_align('V', target, templatesListV, PDBdir, ALIdir + '/', SRULdir, DSSPdirectory,Dsspexe,TMalignexe)

            if rmsd_valIII == 0.0 or rmsd_valV == 0.0 or rmsd_valIV == 0.0 or rmsd_valIII == 999 or rmsd_valV == 999 or rmsd_valIV == 999:
                order = identify_class(rmsd_valIII, rmsd_valIV, rmsd_valV, False, True)
            else:
                order = identify_class(float(tmIII) * naliIV, float(tmIV) * naliIV, float(tmV) * naliIV, True, False)

        except  Exception, e:
            logging.error("Unable to calculate the order for the prediction process"+ str(e))

        if templatesTMdataV == templatesTMdataIV and templatesTMdataIII == templatesTMdataIV and templatesTMdataIII == []:
            logging.info("There are no alignments against the SRUL")
            existaligns = False
        return order, templatesTMdataIII, templatesTMdataIV, templatesTMdataV, existaligns
    except Exception, e:
        logging.error("An error occurred while trying to align target against the SRUL"+ str(e))


def delete_region_folders(folder):
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        if os.path.isdir(file_path):
            shutil.rmtree(file_path)

def create_output_non_predicted(folder, target, Raphaelperiod):
    """
        This function writes an empty db file in case the target if not predicted as a repeat
              :param  : folder : path  in where the db file will be stored
              :param  : target : name of the target
              :param  : Raphaelperiod : Raphael's value for the target
              :type  : folder : string
              :type  : target : string
              :type  : Raphaelperiod : float
    """
    try:
        logging.debug("adding to file create_output_non_predicted " + folder + target + ".db")
        delete_region_folders(folder)
        outf = open(folder + target + ".db", "w")
        outf.write("#Prediction made with ReUPred, using RepeatsDB Lite\n")
        outf.write("SOURCE\tReUPred\n")
        outf.write("PDB\t" + target[:4] + "\n")
        outf.write("CHAIN\t" + target[5:6] + "\n")
        if Raphaelperiod and target[:5] in Raphaelperiod:
            outf.write("RAPHAEL\t" + str(Raphaelperiod[target[:5]]) + "\n")
        else:
            outf.write("RAPHAEL\t*\n")
        outf.write("REPEATED\tNo\n")
        outf.close()
    except Exception, e:
        logging.error("An error occurred while trying to create the empty db file")


def obtain_pdb_input(AllOriginalTargets, AllOriginalTargetschains, quantity, targetname, chain, PDBdir, wdirOriginal):
    try:
        if quantity == 'chain':
            AllOriginalTargets.append(targetname + chain)
            AllOriginalTargetschains.append(chain)
            os.system("mkdir " + wdirOriginal + '_' + chain)
            logging.debug("mkdir " + wdirOriginal + '_' + chain)
            logging.debug('Protein chain to evaluate: ' + targetname + chain)
        elif quantity == 'all':
            if os.path.isfile(PDBdir + targetname):
                pdb = MMCIFParser().get_structure(targetname, PDBdir + targetname)
            elif os.path.isfile(PDBdir + targetname + '.pdb'):
                pdb = MMCIFParser().get_structure(targetname, PDBdir + targetname + '.pdb')
            elif os.path.isfile(PDBdir + 'pdb' + targetname + '.cif'):
                pdb = MMCIFParser().get_structure(targetname, PDBdir + 'pdb' + targetname + '.cif')
            logging.debug('All chains to evaluate in protein ' + targetname)
            for chainN in pdb.get_chains():
                if targetname + chainN.get_id() not in AllOriginalTargets:
                    AllOriginalTargets.append(targetname + chainN.get_id())
                    AllOriginalTargetschains.append(chainN.get_id())
                    chain = targetname + chainN.get_id()
                    os.system("mkdir " + wdirOriginal + '_' + chainN.get_id())
                    logging.debug("mkdir " + wdirOriginal + '_' + chainN.get_id())
                    logging.debug('Protein chain to evaluate: ' + chain)
        elif quantity == 'first':
            if os.path.isfile(PDBdir + targetname):
                pdb = MMCIFParser().get_structure(targetname, PDBdir + targetname)
            elif os.path.isfile(PDBdir + targetname + '.pdb'):
                pdb = MMCIFParser().get_structure(targetname, PDBdir + targetname + '.pdb')
            elif os.path.isfile(PDBdir + 'pdb' + targetname + '.cif'):
                pdb = MMCIFParser().get_structure(targetname, PDBdir + 'pdb' + targetname + '.cif')
            first = True
            for chainN in pdb.get_chains():
                if targetname + chainN.get_id() not in AllOriginalTargets and first:
                    AllOriginalTargets.append(targetname + chainN.get_id())
                    AllOriginalTargetschains.append(chainN.get_id())
                    chain = targetname + chainN.get_id()
                    os.system("mkdir " + wdirOriginal + '_' + chainN.get_id())
                    logging.debug("mkdir " + wdirOriginal + '_' + chainN.get_id())
                    first = False
                    logging.debug('Protein first chain to evaluate: ' + chain)
    except Exception, e:
        logging.error("An error occurred while trying to retrieve the pdb chains |" + targetname)
        logging.warning("An error occurred while trying to retrieve the pdb chains" + targetname)

def create_directory(directory):
    try:
        if not os.path.isdir(directory):
            os.makedirs(directory)
    except  Exception, e:
            logging.error( "was not able no find the directory | " + directory)


def create_all_directories(PDBregionDir, ALIdir, OUTdir):
    create_directory(PDBregionDir)
    create_directory(ALIdir + '/')
    create_directory(OUTdir)
    create_directory(OUTdir + 'III')
    create_directory(OUTdir + 'IV')
    create_directory(OUTdir + 'V')


def create_folders_tree(musdir, wdir, TEMPLATESround, TEMPLATESroundIII, TEMPLATESroundIV, TEMPLATESroundV, DSSPdirectory, OUTdirDB, tempName, PDBregionDir, ALIdir, OUTdir):
    create_all_directories(PDBregionDir, ALIdir, OUTdir)
    create_directory(wdir+musdir)
    create_directory(OUTdirDB)
    create_directory(DSSPdirectory)
    create_directory(wdir + "dsspRegion/")
    for i in range(4):
        create_templates_folders(wdir + tempName + str(i + 1), TEMPLATESround)
        create_templates_folders(wdir + tempName + str(i + 1), TEMPLATESroundIII)
        create_templates_folders(wdir + tempName + str(i + 1), TEMPLATESroundIV)
        create_templates_folders(wdir + tempName + str(i + 1), TEMPLATESroundV)


def get_pdb_information(pdbfile):
    handle = open(pdbfile, 'r')
    header_dict = parse_pdb_header(handle)
    handle.close()
    return header_dict

def obtain_header_sec_str(pdbfile, folder):

    nlineT = subprocess.check_output(
        'grep \"^ATOM\" -n -m 1 '  + folder+pdbfile + ' | cut -d : -f 1', shell=True)
    nline=int(nlineT)-1

    pdbhead = subprocess.check_output(
        'sed -n 1,'+str(int(nline))+ 'p ' + folder + pdbfile, shell=True)

    headefile = open(folder+'header','w')
    for lines in pdbhead:
        headefile.write(lines)
    headefile.close()


def download_pdb_chains(entry, pdbFolder, exist, quantity, chainUsr,pdbFolderDB,wdirOriginal):
    try:

        io = PDBIO()
        gotit ,found = False,False
        pdb_header = {}

        if exist:
            pdb = PDBParser().get_structure(entry.upper(), pdbFolder + entry)
            obtain_header_sec_str(entry, pdbFolder)
            pdb_header= get_pdb_information(pdbFolder + entry)
            print pdb_header
            if quantity == 'first':

                for chain in pdb.get_chains():

                    pdb_header['chain '+chain.get_id()] = 'non requested'
                    if not gotit:
                        chainUsr = chain.get_id()
                        gotit = True
                        pdb_header['chain ' + chainUsr] = 'requested'

            if quantity == 'all':


                for chainO in pdb.get_chains():
                    pdb_header['chain ' + chainO.get_id()] = 'requested'
                    chain = entry + chainO.get_id()
                    io.set_structure(chainO)
                    os.system('mkdir ' + pdbFolder + '_' + chainO.get_id())
                    io.save(pdbFolder + chain+ '.pdb')
                    io.save(pdbFolder +'_'+chainO.get_id()+'/'+  entry +'_'+ chainO.get_id() + '.pdb')

                    create_mapping(pdbFolder + chain + '.pdb', pdbFolder +'_'+chainO.get_id()+'/'+  entry +'_'+ chainO.get_id() + '.mapping', chain)
                    files = [wdirOriginal + "header",pdbFolder +'_'+chainO.get_id()+'/'+  entry +'_'+ chainO.get_id() + '.pdb']
                    concat = ''.join([open(f).read() for f in files])
                    auxfile = open(pdbFolder +'_'+chainO.get_id()+'/'+  entry +'_'+ chainO.get_id() + '.pdb', 'w')
                    auxfile.write(concat)
                    auxfile.close()

            elif quantity == 'chain' or quantity == 'first':
                chain = entry + chainUsr
                flag=True

                for chainO in pdb.get_chains():
                        pdb_header['chain ' + chainO.get_id()] = 'non requested'
                        if  flag and quantity == 'first':
                            io.set_structure(chainO)
                            flag = False
                        if quantity == 'chain' and chainO.get_id() == chainUsr:
                            io.set_structure(chainO)

                os.system('mkdir ' + pdbFolder + '_' + chainO.get_id())
                io.save(pdbFolder+'_'+chainO.get_id()+'/'+  entry +'_'+ chainO.get_id() + '.pdb')
                io.save(pdbFolder+ chain + '.pdb')
                create_mapping(pdbFolder + chain+ '.pdb', pdbFolder +'_'+chainO.get_id()+'/'+  entry +'_'+ chainO.get_id()+ '.mapping', chain)
                files = [wdirOriginal + "header",pdbFolder+'_'+chainO.get_id()+'/'+  entry +'_'+ chainO.get_id() + '.pdb']
                concat = ''.join([open(f).read() for f in files])
                auxfile = open(pdbFolder+'_'+chainO.get_id()+'/'+  entry +'_'+ chainO.get_id() + '.pdb', 'w')
                auxfile.write(concat)
                auxfile.close()

        else:
            if entry != '' or os.stat(pdbFolder + entry + '.pdb').st_size < 1:
                k = 1
                flag = True
                pl = PDBList()

                #print pdbFolderDB + entry[1:3] + '/' + 'pdb' + entry[:4] + '.ent.gz'
                if not os.path.isfile(pdbFolderDB +entry[1:3]+'/' + 'pdb' + entry[:4] + '.cif'):
                    # try:
                    #     print entry[:4]

                    # except Exception, e:

                        logging.info("there is no ent file | " +pdbFolderDB +entry[1:3]+'/' + 'pdb' + entry[:4] + '.ent' + ' trying other options to retrieve')
                        url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId="
                        pdbid = url + entry
                        open(pdbFolder+'pdb' + entry[:4] + '.cif' , "w").write(urllib.urlopen(pdbid).read())
                        open( pdbFolder+ entry[:4] + '.pdb', "w").write(urllib.urlopen(pdbid).read())
                        pdb_header = get_pdb_information(pdbFolder+'pdb' + entry[:4] + '.cif')
                        pdbname = 'pdb' + entry[:4] + '.cif'
                        pdb = PDBParser().get_structure(entry[:4].upper(), pdbFolder + 'pdb' + entry[:4] + '.cif')
                        obtain_header_sec_str( 'pdb' + entry[:4] + '.cif', pdbFolder)
                        os.system('cp ' + pdbFolder + 'pdb' + entry + '.cif ' + pdbFolder + 'reupred_' + entry)
                else:
                    os.system('cp ' + pdbFolderDB +entry[1:3]+'/' + 'pdb' + entry[:4] + '.cif '+pdbFolder+ 'pdb' + entry + 'cif ' )
                    os.system('gzip -d '+pdbFolder+ 'pdb' + entry + '.cif '+pdbFolder+ 'pdb' + entry + '.cif')

                    pdb = PDBParser().get_structure(entry.upper(), pdbFolder + 'pdb' + entry + '.cif')
                    os.system('cp '+pdbFolder + 'pdb' + entry + '.cif '+pdbFolder+'reupred_'+entry)


                    pdb_header = get_pdb_information(pdbFolder + 'pdb' + entry + '.cif')
                    obtain_header_sec_str('pdb' + entry + '.cif', pdbFolder)
                    pdbname = 'pdb' + entry + '.cif'
                if quantity == 'first':
                    for chain in pdb.get_chains():
                        pdb_header['chain ' + chain.get_id()] = 'non requested'
                        if not gotit:
                            chainUsr = chain.get_id()
                            gotit = True
                if quantity == 'all':


                    for chainO in pdb.get_chains():
                        pdb_header['chain ' + chainO.get_id()] = 'non requested'
                        chain = entry + chainO.get_id()
                        io.set_structure(chainO)

                        os.system('mkdir ' + pdbFolder + '_' + chainO.get_id())
                        io.save(pdbFolder + '_' + chainO.get_id() + '/' +  entry +'_'+ chainO.get_id() + '.pdb')
                        io.save(pdbFolder +  chain + '.pdb')
                        create_mapping(pdbFolder + chain + '.pdb', pdbFolder +'_'+chainO.get_id()+'/'+  entry +'_'+ chainO.get_id()+ '.mapping', chain)
                        files = [wdirOriginal + "header", pdbFolder + '_' + chainO.get_id() + '/' +  entry +'_'+ chainO.get_id() + '.pdb']
                        concat = ''.join([open(f).read() for f in files])
                        auxfile = open(pdbFolder + '_' + chainO.get_id() + '/' +  entry +'_'+ chainO.get_id() + '.pdb', 'w')
                        auxfile.write(concat)
                        auxfile.close()
                elif quantity == 'chain' or quantity == 'first':
                    flag = True
                    for chainO in pdb.get_chains():
                        pdb_header['chain ' + chainO.get_id()] = 'non requested'
                        if  flag and quantity == 'first':
                            io.set_structure(chainO)
                            flag = False
                            selected = chainO.get_id()
                        if quantity == 'chain' and chainO.get_id() == chainUsr:
                            io.set_structure(chainO)
                            found= True
                            selected = chainO.get_id()
                    if not found:
                        logging.warning("The chain does not exist in the structure protein"  )
                    chain = entry + chainUsr

                    os.system('mkdir ' + pdbFolder + '_' + selected)
                    io.save(pdbFolder + '_' + selected + '/' +  entry +'_'+ selected+ '.pdb')
                    io.save(pdbFolder + chain + '.pdb')
                    create_mapping(pdbFolder + chain + '.pdb', pdbFolder +'_'+selected+'/'+  entry +'_'+ selected + '.mapping', chain)
                    files = [wdirOriginal + "header", pdbFolder + '_' + selected + '/' +  entry +'_'+ selected+ '.pdb']
                    concat = ''.join([open(f).read() for f in files])
                    auxfile = open(pdbFolder + '_' + selected + '/' +  entry +'_'+ selected+ '.pdb', 'w')
                    auxfile.write(concat)
                    auxfile.close()
            else:
                logging.warning("there is no input file or the selected pdb is tooo short")

        successpdb = True

    except  Exception, e:
        logging.error("Error retrieving PDB file" + str(e))
        successpdb = False
    return chainUsr, successpdb, pdb_header


def create_mapping(pdbfilename, mappfilename, pdbId):
    if os.path.isfile(pdbfilename):
        pdbFile = open(pdbfilename,'r')
        mappfile = open(mappfilename,'w')
        parser = PDBParser()
        structureAux = parser.get_structure(pdbId[:4], pdbFile)
        model = structureAux[0]
        secstr = ''
        unable = False
        try:
            dssp = DSSP(model, pdbfilename)
            for dsspvals in list(dssp.keys()):
                secstr = secstr + dssp[dsspvals][2]
        except :
            logging.warning("Unable to obtain secondary structure |"+pdbfilename)
            unable = True
        dsspunit = ''
        count = 0
        for residue in structureAux[0][pdbId[4:5]].get_residues():
            if residue.get_resname() in ResiduesCodes and  (secstr=='' or count<len(secstr)):
                if secstr!='' and not unable:
                    mappfile.write(str(count) + ' ' + str(residue.get_id()[1]) + ' ' + ResiduesCodes[residue.get_resname()] + ' ' + secstr[count] + '\n')
                else:
                    mappfile.write(
                        str(count) + ' ' + str(residue.get_id()[1]) + ' ' + ResiduesCodes[residue.get_resname()] + ' ' +
                         '- \n')
                count += 1
            elif unable and residue.get_resname() in ResiduesCodes :
                mappfile.write(
                    str(count) + ' ' + str(residue.get_id()[1]) + ' ' + ResiduesCodes[residue.get_resname()] + ' ' +
                    '.' + '\n')



def create_folders_target(directory, templatesWork , target, PDBdir):
    if not os.path.isdir(directory + templatesWork):
        os.makedirs(directory + templatesWork)
    if not os.path.isdir(directory + templatesWork + target[:5]):
        os.makedirs(directory + templatesWork + target[:5])
    if os.path.isdir(directory + templatesWork + target[:5]):
        os.system('cp ' + PDBdir + target[:5] + '.pdb ' + directory + templatesWork + target[:5] + '/')
        logging.debug('cp ' + PDBdir + target[:5] + '.pdb ' + directory + templatesWork + target[:5] + '/')
    if not os.path.isdir(directory + templatesWork + 'units'):
        os.makedirs(directory + templatesWork + 'units')


def move_files(origenUnitPath, target, destinationPath, origenPath):
    # moves all the predicted files to auxiliar folders
    os.system('mv '+origenUnitPath + target[:5] + '*repeat* ' + destinationPath)
    os.system('mv '+origenPath + target[:5]+'.len ' + destinationPath)
    os.system('mv '+origenPath + target[:5]+'.fastp ' + destinationPath)
    os.system('mv '+origenPath + target[:5]+'.mus ' + destinationPath)
    logging.debug('moving files' + ' mv '+origenUnitPath + target[:5] + '*repeat* ' + destinationPath)

def identify_classification_class_v(target, directory, templatesTMdata, inSecondRound, strUnits, region, returnVal, strins):
    sc = "0 - - 0\n"
    descsc = "\tUnknown "
    data_directory = directory + 'Data/'
    with open(data_directory + "templateListV1", 'r') as f:
        AlphaTemplates = f.read().split("\n")
    with open(data_directory + "templateListV2", 'r') as f:
        BetaTemplates = f.read().split("\n")
    with open(data_directory + "templateListV3", 'r') as f:
        AlphaBetaTemplates = f.read().split("\n")
    with open(data_directory + "templateListV4", 'r') as f:
        BetaSandwichTemplates = f.read().split("\n")
    with open(data_directory + "templateListV5", 'r') as f:
        AlphaBetaSandwichTemplates = f.read().split("\n")
    if templatesTMdata[inSecondRound][0] in AlphaTemplates:
        descsc = "\tPredicted subclass Alpha Beads "
        sc = "1 - - - -\n"
    if templatesTMdata[inSecondRound][0] in BetaTemplates:
        descsc = "\tPredicted subclass Beta Beads "
        sc = "2 - - - -\n"
    if templatesTMdata[inSecondRound][0] in AlphaBetaTemplates:
        descsc = "\tPredicted subclass Alpha/Beta Beads "
        sc = "3 - - - -\n"
    if templatesTMdata[inSecondRound][0] in BetaSandwichTemplates:
        descsc = "\tPredicted subclass Beta Sandwich Beads "
        sc = "4 - - - -\n"
    if templatesTMdata[inSecondRound][0] in AlphaBetaSandwichTemplates:
        descsc = "\tPredicted subclass Alpha Beta Sandwich Beads "
        sc = "5 - - - -\n"
    if returnVal:
        return target[:4] + "\t" + target[4:5] + "\t" + str(region.count(';') + 1) + "\t" + region + "\t" + strUnits[1:] + "\t" + strins + "\t" + "V" + "\t" + sc, "V " + sc
    else:
        return "", "V " + sc


def identify_classification_class_iii(target, directory, templatesTMdata, inSecondRound, strUnits, region, returnVal, strins):
    sc = "0 - - 0\n"
    unitype = ""
    descsc = "\tUnknown "
    data_directory=directory+'Data/'
    with open(data_directory + "templateListIII1", 'r') as f:
        BetaTemplates = f.read().split("\n")
    with open(data_directory + "templateListIII2", 'r') as f:
        AlphaBetaTemplates = f.read().split("\n")
    with open(data_directory + 'clusterTpr', 'r') as f:
        templateListTPRlike = f.read().split("\n")
    with open(data_directory + 'clusterAnkyrin', 'r') as f:
        templateListANKlike = f.read().split("\n")
    with open(data_directory +'clusterArmadillo', 'r') as f:
        templateListARMlike = f.read().split("\n")
    with open(data_directory + 'clusterTal', 'r') as f:
        templateListTALlike = f.read().split("\n")
    with open(data_directory + 'clusterPumilio', 'r') as f:
        templateListPUMlike = f.read().split("\n")

    with open(data_directory + "templateListIII4", 'r') as f:
        spiralTemplates = f.read().split("\n")
    with open(data_directory + "templateListIII5", 'r') as f:
        singleLayerTemplates = f.read().split("\n")

    with open(data_directory + "templateListIII6", 'r') as f:
        Box = f.read().split("\n")
    if templatesTMdata[inSecondRound][0] in BetaTemplates:
        descsc = "\tPredicted subclass Beta solenoid "
        sc = "1 - - - -\n"
    if templatesTMdata[inSecondRound][0] in AlphaBetaTemplates:
        descsc = "\tPredicted subclass Alpha/Beta solenoid "
        sc = "2 - - - -\n"
    if templatesTMdata[inSecondRound][0] in templateListTPRlike or templatesTMdata[inSecondRound][0] in templateListANKlike or templatesTMdata[inSecondRound][0] in templateListARMlike  or templatesTMdata[inSecondRound][0] in templateListPUMlike  or templatesTMdata[inSecondRound][0] in templateListTALlike:
        descsc = "\tPredicted subclass Alpha solenoid "
        sc = "3"
        with open(data_directory + "templateAnkyrin", 'r') as f:
            ankyrin = f.read().split("\n")
        with open(data_directory + "templateSel1", 'r') as f:
            sel1 = f.read().split("\n")
        with open(data_directory + "templatePpta", 'r') as f:
            ppta = f.read().split("\n")
        with open(data_directory + "templateTpr", 'r') as f:
            tpr = f.read().split("\n")
        with open(data_directory + "templateHeat", 'r') as f:
            heat = f.read().split("\n")
        with open(data_directory + "templateArmadillo", 'r') as f:
            armadillo = f.read().split("\n")
        with open(data_directory + "templateImportin", 'r') as f:
            importin = f.read().split("\n")
        with open(data_directory + "templateTal", 'r') as f:
            tal = f.read().split("\n")
        with open(data_directory + "templatePumilio", 'r') as f:
            pumilio = f.read().split("\n")
        with open(data_directory + "templateAnapc", 'r') as f:
            anapc = f.read().split("\n")
        with open(data_directory + "templateCrm1c", 'r') as f:
            crm1c = f.read().split("\n")
        with open(data_directory + "templateAdaptinN", 'r') as f:
            adaptinn = f.read().split("\n")
        with open(data_directory + "templateXpo1", 'r') as f:
            xpo1 = f.read().split("\n")
        if templatesTMdata[inSecondRound][0] in ankyrin:
            descsc = descsc + " ankyrin "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in sel1:
            descsc = descsc + " sel1 "
            sc = sc + " - - -"
            unitype = unitype + " -\n"

        elif templatesTMdata[inSecondRound][0] in heat:
            descsc = descsc + " heat "
            sc = sc + " - - -"
            unitype = unitype + " -\n"

        elif templatesTMdata[inSecondRound][0] in ppta:
            descsc = descsc + " ppta "
            sc = sc + " - - -"
            unitype = unitype + " -\n"

        elif templatesTMdata[inSecondRound][0] in tpr:
            descsc = descsc + " tpr"
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in anapc:
            descsc = descsc + " anapc "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in adaptinn:
            descsc = descsc + " adaptinN "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in tpr:
            descsc = descsc + " Xpo1 "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in crm1c:
            descsc = descsc + " crm1c "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in tal:
            descsc = descsc + " tal "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in armadillo:
            descsc = descsc + " armadillo "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in importin:
            descsc = descsc + " importin "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in pumilio:
            descsc = descsc + " pumilio "
            sc = sc + " -"
            unitype = unitype + " -\n"
        else:
            descsc = descsc + " other "
            sc = sc + " 0"
            unitype = unitype + " 0\n"


    if templatesTMdata[inSecondRound][0] in spiralTemplates:
        descsc = "\tPredicted subclass trimer beta spirals  "
        sc = "4\n"
        with open(data_directory + "templateOstA", 'r') as f:
            osta = f.read().split("\n")
        with open(data_directory + "templateGlycoshydro", 'r') as f:
            glycohydro = f.read().split("\n")
        with open(data_directory + "templateCwBinding", 'r') as f:
            cwbinding = f.read().split("\n")
        with open(data_directory + "templateLactamase", 'r') as f:
            lactamase = f.read().split("\n")
        if templatesTMdata[inSecondRound][0] in osta:
            descsc = descsc + " OstA "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in glycohydro:
            descsc = descsc + " Glycohydro "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in lactamase:
            descsc = descsc + " Lactamase "
            sc = sc + " - - -"
            unitype = unitype + " -\n"

        elif templatesTMdata[inSecondRound][0] in cwbinding:
            descsc = descsc + " CW_binding "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        else:
            descsc = descsc + " other "
            sc = sc + " 0"
            unitype = unitype + " 0\n"


    if templatesTMdata[inSecondRound][0] in singleLayerTemplates :
        descsc = "\tPredicted subclass single layer beta "
        sc = "5\n"
        with open(data_directory + "templateLipoprotein1", 'r') as f:
            lipoprotein1 = f.read().split("\n")
        with open(data_directory + "templateTcp10c", 'r') as f:
            tcp10c = f.read().split("\n")
        if templatesTMdata[inSecondRound][0] in tcp10c:
            descsc = descsc + " Tcp10_c "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in lipoprotein1:
            descsc = descsc + " Lipoprotein_1 "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in osta:
            descsc = descsc + " OstA "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        else:
            descsc = descsc + " other "
            sc = sc + " -"
            unitype = unitype + " -\n"

    if templatesTMdata[inSecondRound][0] in Box:
        descsc = "\tPredicted subclass Box"
        sc = "6 - - - -\n"
    if returnVal:
        return target[:4] + "\t" + target[4:5] + "\t" + str(region.count(';') + 1) + "\t" + region + "\t" + strUnits[1:] + "\t" + strins + "\t" + "III" + "\t" + sc, "III." + sc
    else:
        return "", "III " + sc + unitype



def identify_classification_class_iv(target, directory, templatesTMdata, inSecondRound, strUnits, region, returnVal, strins):
    sc = "0 - - - 0\n"
    unitype = ""
    data_directory = directory+'Data/'
    descsc = "\tUnknown "
    with open(data_directory + "templateListIV1", 'r') as f:
        TimBarrelTemplates = f.read().split("\n")
    with open(data_directory + 'clusterLipocalin', 'r') as f:
        LIPOCALINlike = f.read().split("\n")
    with open(data_directory + 'clusterOstac', 'r') as f:
        OSTAClike = f.read().split("\n")
    with open(data_directory + "templateListIV3", 'r') as f:
        BetaTrefoilTemplates = f.read().split("\n")
    with open(data_directory + "templateListIV4", 'r') as f:
        BetaPropellerTemplates = f.read().split("\n")
    with open(data_directory + "templateListIV5", 'r') as f:
        AlphaBetaPrismTemplates = f.read().split("\n")
    with open(data_directory + "templateListIV6", 'r') as f:
        AlphaBarrelTemplates = f.read().split("\n")
    with open(data_directory + "templateListIV7", 'r') as f:
        AlphaBetaBarrel = f.read().split("\n")
    with open(data_directory + "templateListIV8", 'r') as f:
        AlphaBetaPropeller = f.read().split("\n")
    with open(data_directory + "templateListIV9", 'r') as f:
        AlphaBetaTrefoilTemplates = f.read().split("\n")
    with open(data_directory + "templateListIV10", 'r') as f:
        AlignPrism = f.read().split("\n")

    if templatesTMdata[inSecondRound][0] in TimBarrelTemplates:
        descsc = "\tPredicted subclass Tim Barrel "
        sc = "1 - - - -\n"
    if templatesTMdata[inSecondRound][0] in LIPOCALINlike or  templatesTMdata[inSecondRound][0] in OSTAClike:
        descsc = "\tPredicted subclass Beta Barrel "
        sc = "2\n"
        with open(data_directory + "clusterOstac", 'r') as f:
            ostac = f.read().split("\n")
        with open(data_directory + "clusterLipocalin", 'r') as f:
            lipocalin = f.read().split("\n")
        with open(data_directory + "templatetonBdeprec", 'r') as f:
            tonb = f.read().split("\n")
        with open(data_directory + "templatePorin", 'r') as f:
            porin = f.read().split("\n")
        if templatesTMdata[inSecondRound][0] in lipocalin:
            descsc = descsc + " Lipocalin "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in porin:
            descsc = descsc + " Porin "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        elif templatesTMdata[inSecondRound][0] in tonb:
            descsc = descsc + " TonBdeprec "
            sc = sc + " - - -"
            unitype = unitype + " -\n"
        else:
            descsc = descsc + " other "
            sc = sc + " 0"
            unitype = unitype + " 0\n"
    if templatesTMdata[inSecondRound][0] in BetaTrefoilTemplates:
        descsc = "\tPredicted subclass Beta Trefoil "
        sc = "3 - - - -\n"
    if templatesTMdata[inSecondRound][0] in BetaPropellerTemplates:
        descsc = "\tPredicted subclass Beta Propeller "
        sc = "4 - - - -\n"
    if templatesTMdata[inSecondRound][0] in AlphaBetaPrismTemplates:
        descsc = "\tPredicted subclass Alpha Beta Prism "
        sc = "5 - - - -\n"
    if templatesTMdata[inSecondRound][0] in AlphaBarrelTemplates:
        descsc = "\tPredicted subclass Alpha Barrel "
        sc = "6 - - - -\n"
    if templatesTMdata[inSecondRound][0] in AlphaBetaBarrel:
        descsc = "\tPredicted subclass Alpha Beta Barrel "
        sc = "7 - - - -\n"
    if templatesTMdata[inSecondRound][0] in AlphaBetaPropeller:
        descsc = "\tPredicted subclass Alpha Beta Propeller "
        sc = "8 - - - -\n"
    if templatesTMdata[inSecondRound][0] in AlphaBetaTrefoilTemplates:
        descsc = "\tPredicted subclass Alpha Beta Trefoil "
        sc = "9 - - - -\n"
    if templatesTMdata[inSecondRound][0] in AlignPrism:
        descsc = "\tPredicted subclass Aligned Prism "
        sc = "10 - - - -\n"
    if returnVal:
        return target[:4] + "\t" + target[4:5] + "\t" + str(region.count(';') + 1) + "\t" + region + "\t" + strUnits[1:] + "\t" + strins + "\t" + "IV" + "\t" + sc, "IV." + sc
    else:
        return "", "IV " + sc


def create_final_output(target, OUTdirDB, region, fastpListAA, Raphaelperiod, strins, sc, difreg, selectedTemplate):

    logging.debug("adding to file create_final_output " +OUTdirDB + target[:5] + ".db")
    if not difreg:
        outf = open(OUTdirDB + target[:5] + ".db", "a")
    else:
        outf = open(OUTdirDB + target[:5] + ".db", "w")
        outf.write("SOURCE\tReUPred\n")
        outf.write("PDB\t" + target[:4] + "\n")
        outf.write("CHAIN\t" + target[4:5] + "\n")
        if Raphaelperiod and target[:5] in Raphaelperiod:
            outf.write("RAPHAEL\t" + str(Raphaelperiod[target[:5]]) + "\n")
        else:
            outf.write("RAPHAEL\t*\n")
    reg = region.split(";")
    for r in reg:
        outf.write("REG\t" + r + " " + sc +"\n")
        valuesMaster = selectedTemplate.split('_')
        outf.write("MASTER\t" + r + " " + valuesMaster[0] + ' ' + valuesMaster[1] + ' ' + valuesMaster[2] + "\n")

    for r in fastpListAA:
        uni = r.split("\t")
        outf.write("UNIT\t" + uni[0] + " " + uni[1] + "\n")

    if strins!='':
        ins = strins.split(";")
        for r in ins:
            if len(r)>3:

                outf.write("INS\t" + r + "\n")
    outf.close()


def move_files_finaldestination(target, wdirOriginal, chain, OriginalTargetname):
    os.system('cp ' + wdirOriginal + target + 'Original.pdb ' + wdirOriginal + '_' + chain + '/' + OriginalTargetname + '_' + chain + '.pdb ')
    os.system('mv ' + wdirOriginal + target + 'Original.pdb ' + wdirOriginal + OriginalTargetname + '_' + chain + '.pdb ')
    os.system('mv ' + wdirOriginal + target + '.mapping ' + wdirOriginal + '_' + chain + '/' + OriginalTargetname + '_' + chain + '.mapping ')
    if os.path.isfile(wdirOriginal + target + '.pdb'):
        os.remove(wdirOriginal + target + '.pdb')
        os.remove(wdirOriginal + OriginalTargetname + '_' + chain + '.pdb')
    logging.debug("moving files final " + 'cp ' + wdirOriginal + target + 'Original.pdb ' + wdirOriginal + '_' + chain + '/' + OriginalTargetname + '_' + chain + '.pdb ')


def verify_region(target, sc, PDBdir, PDBregionDir):

    if os.path.isfile(sc + target + ".fastp"):
        with open(sc + target + ".fastp", 'r') as f:
            unitList = f.read().split("\n")
        unitList, rr ,insertions,ninsertions= sort_fastp_result(list(unitList))

        minUnit = unitList[0].split("\t")[0]
        maxUnit = unitList[len(unitList) - 1].split("\t")[1]
    created = create_aux_pdb(int(minUnit), int(maxUnit), PDBdir, PDBregionDir, target, target)
    if created:
        return True
    return False


def create_aux_pdb( min_, max_, origen, destination, target,   tempName):
        fpdbUnitpre = open(destination+tempName + "pre.pdb", 'w')
        fpdbUnitpost = open(destination+tempName + "post.pdb", 'w')
        if os.path.isfile(origen + target[:5] + '/' + target[:5] + ".pdb"):
            pdbFile = open(origen + target[:5] + '/' + target[:5] + ".pdb", 'r')
        elif os.path.isfile(origen + '/' + target[:5] +".pdb"):
            pdbFile = open(origen + '/' + target[:5] + ".pdb", 'r')
        else:
            pdbFile = open(origen + '/' + target, 'r')
        parser = PDBParser()
        k = 1
        countAuxpost,countAuxpre,countAux = 0.0,0.0,0.0
        structureAux = parser.get_structure(target[:4], pdbFile)
        for res in structureAux[0][target[4:5]].get_residues():
            if res.get_id()[0] == ' ':
                for a in res:
                    if len(a.get_name()) == 4:
                        spc = ''
                    elif len(a.get_name()) == 1:
                        spc = '   '
                    elif len(a.get_name()) == 2:
                        spc = '  '
                    else:
                        spc = ' '
                    spp = '     '
                    if int(res.get_id()[1]) not in range(int(min_), int(max_)+1):
                        if int(res.get_id()[1]) < int(min_):
                            if k >= 10000:
                                    fpdbUnitpre.write('ATOM  ' + '%5s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + '  ' + target[4:5] + '%4s'%str(res.get_id()[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                            else:
                                    fpdbUnitpre.write('ATOM   ' + '%4s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id()[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                            k += 1
                            if a.get_name() == 'C':
                                countAuxpre += 1
                        elif int(res.get_id()[1]) > int(max_) + 1:
                            if k >= 10000:
                                    fpdbUnitpost.write('ATOM  ' + '%5s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + '  ' + target[4:5] + '%4s'%str(res.get_id()[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                            else:
                                    fpdbUnitpost.write('ATOM   ' + '%4s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id()[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0])+'%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                            k += 1
                            if a.get_name() == 'C':
                                countAuxpost += 1
                    if a.get_name() == 'C':
                        countAux += 1

        fpdbUnitpost.close()
        fpdbUnitpre.close()
        pdbFile.close()
        if float(countAuxpre/countAux) > float(countAuxpost/countAux) and float(countAuxpre/countAux) > 0.3 :
            os.system('rm ' + destination + tempName + "post.pdb")
            os.system('mv ' + destination + tempName + "pre.pdb " + destination +tempName+".pdb")
            return True
        elif float(countAuxpre/countAux) < float(countAuxpost/countAux) and float(countAuxpost/countAux) > 0.3:
            os.system('rm ' + destination + tempName + "pre.pdb")
            os.system('mv ' + destination + tempName + "post.pdb " + destination + tempName+".pdb")
            return True
        else:
            os.system('rm ' + destination+ tempName+"pre.pdb")
            os.system('rm ' + destination+ tempName+"post.pdb")
        return False


def create_templates_folders(directory, TEMPLATESround):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    TEMPLATESround.append(directory)


def identify_class(templatesTMdataIII, templatesTMdataIV, templatesTMdataV, ord_, eliminateone):
    rmsdList , rmsdListsort,newList_rmsd = [], [], []
    rmsdList.append(float(templatesTMdataIII)+0.0001)
    rmsdList.append(float(templatesTMdataIV)+0.00011)
    rmsdList.append(float(templatesTMdataV)+0.000111)

    count = 0

    if eliminateone== True and  ord_ ==False:
        newList_rmsd = sorted(rmsdList, reverse=ord_)

        rmsdListsort = [rmsdList.index(newList_rmsd[0]), rmsdList.index(newList_rmsd[1]), rmsdList.index(newList_rmsd[2])]

        if templatesTMdataIII == 999:
            count = count + 1
        if templatesTMdataIV == 999:
            count = count + 1
        if templatesTMdataV == 999:
            count = count + 1

        if count == 1    :
            return list([int(rmsdListsort[0]), int(rmsdListsort[1])])
        elif count == 2:
            return list([int(rmsdListsort[0])])
        else:
            return list([int(rmsdListsort[0]), int(rmsdListsort[1]), int(rmsdListsort[2])])

    elif eliminateone == True and ord_ == True:
        rmsdListsort = [rmsdList.index(sorted(rmsdList, reverse=ord_)[0]),
                        rmsdList.index(sorted(rmsdList, reverse=ord_)[1]),
                        rmsdList.index(sorted(rmsdList, reverse=ord_)[2])]
        if templatesTMdataIII == 0 or templatesTMdataIII == 999:
            count = count + 1
        if templatesTMdataIV == 0 or templatesTMdataIV == 999:
            count = count + 1
        if templatesTMdataV == 0 or templatesTMdataV == 999:
            count = count + 1
        if count == 1:
            return list([int(rmsdListsort[0]), int(rmsdListsort[1])])
        elif count == 2:
            return list([int(rmsdListsort[2])])
        else:
            return list([int(rmsdListsort[0]), int(rmsdListsort[1]), int(rmsdListsort[2])])
    elif (eliminateone == False and ord_ == True ) or ( eliminateone == False and ord_ == False):
        rmsdListsort = [rmsdList.index(sorted(rmsdList, reverse=ord_)[0]),
                        rmsdList.index(sorted(rmsdList, reverse=ord_)[1]),
                        rmsdList.index(sorted(rmsdList, reverse=ord_)[2])]
        return list([int(rmsdListsort[0]), int(rmsdListsort[1]), int(rmsdListsort[2])])


def load_dictionary(Raphaelperiod, directory):
    f = open(directory+'Data/raphael.output')
    for lines in f:
        items = lines.split('\t')
        Raphaelperiod[items[0]] = eval(items[1])
    return Raphaelperiod


def avoid_holes(fastpListAA2, filename,dsspexe):
    secstr = ''
    p = PDBParser()

    structure = p.get_structure('unit', filename + '.pdb')
    model = structure[0]
    dssp = DSSP(model, filename + '.pdb',dssp=dsspexe)
    j =0

    fastpListAA = []
    for line in fastpListAA2:
        rr = line.split('\t')
        if j==0:
            prev_ini = int(rr[0])
            prev_end = int(rr[1])
            j=j+1
        elif 14 > int(rr[0]) - int(prev_end) > 1:

            secstr = ''
            done=False
            newstart = int(rr[0])
            newprevend = int(prev_end)
            for dsspvals in list(dssp.keys()):
                if int(dssp[dsspvals][0]) in range(int(prev_end),int(rr[0])):

                    if not done and dssp[dsspvals][2]=='-':
                        done = True

                        newprevend, newstart = dssp[dsspvals][0]-1, dssp[dsspvals][0]

            fastpListAA.append(str(prev_ini)+'\t'+str(newprevend))
            prev_ini = newstart
            prev_end = rr[1]
        else:
            fastpListAA.append(str(prev_ini)+'\t'+str(prev_end))
            prev_ini = int(rr[0])
            prev_end = int(rr[1])
    fastpListAA.append(str(prev_ini)+'\t'+str(prev_end))

    return fastpListAA


def sort_fastp_result(fastpListAA):
    fastpL = []
    insertion=''
    ninsertions=0
    for rr in list(fastpListAA):
        if rr != '':
            row = rr.split('\t')
            max_ = int(row[1])
            min_ = int(row[0])
            fastpL.append([min_, max_])

    fastpL = sorted(fastpL)
    fastpListAA = []
    region = str(fastpL[0][0])+' '
    j=0
    for rr in fastpL:
        if rr:
            fastpListAA.append(str(rr[0])+"\t"+str(rr[1]))
            if j==len(fastpL)-1:
                region = region + str(fastpL[j][1])

            else:
                if j<len(fastpL)-1 and fastpL[j+1][0]-fastpL[j][1]>1:
                    #region = region + str(fastpL[j][1])+';'+str(fastpL[j+1][0])+" "
                    insertion =insertion+str(fastpL[j][1]+1)+' '+str(fastpL[j+1][0]-1) +';'
                    ninsertions = int(fastpL[j+1][0] )-int( fastpL[j][1]) + ninsertions
        j += 1

    return fastpListAA, region, insertion,ninsertions

