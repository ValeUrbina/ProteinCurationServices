from predictorMethods import *
import glob
import ConfigParser
import sys, getopt
from time import time
import logging,timeit

from otherMethods import download_pdb_chains

argv=sys.argv[1:]

try:
    opts, args = getopt.getopt(argv,"hi:q:d:e:c:s:v:",["ifile=","qua=","dir=","exist=","chain=","software=","validate="])
except getopt.GetoptError:
    print ('python reupred.py -d <working path> -i <targetname> -c <chain> -q <all/first/chain> -e <True/False> -s <Softwarepath> -v <True/False>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print ('python reupred.py -d <working path> -i <targetname> -c <chain> -q <all/first/chain> -e <True/False> -s <softwarepath> -v <True/False>')
        sys.exit()
    elif opt in ("-i", "--ifile"):
        targetname = arg
    elif opt in ("-d", "--dir"):
        wdirOriginal = arg
        wdir = arg+"Temporary/"
    elif opt in ("-c", "--chain"):
        chain = arg
        if chain =='-':
            chain =''
        else:
            chain = chain.upper()
    elif opt in ("-q", "--qua"):
        quantity = arg
        quantity=quantity.lower()
    elif opt in ("-e", "--exist"):
        pdbExists = (arg)
        pdbExists = pdbExists.upper()
        #print pdbExists

    elif opt in ("-s", "--software"):
        directory=arg
    elif opt in ("-v", "--validate"):
        validdb = (arg)
        validdb = validdb.upper()
try:
    initial_time = time()
    #  if os.path.isfile(wdirOriginal + 'reupred.log'):
    #     os.system('rm '+wdirOriginal + 'reupred.log')
    LOG = wdirOriginal + 'reupred.log'
    logging.basicConfig(filename=LOG, level=logging.INFO, format='%(asctime)s | %(name)s | %(threadName)s | %(levelname)s | %(message)s', datefmt='%d/%m/%Y %H:%M:%S')
    logging.info("Starting .....")
    if pdbExists!='TRUE':
        targetname=targetname.lower()
    #print targetname[:8]
    if targetname[:8]=='reupred_':
        #print targetname[8:]

        os.system('cp ' + wdirOriginal + targetname + ' ' + wdirOriginal +  targetname[8:])
        targetname = targetname[8:]

    logging.info("Running predictor:" + ' ' + targetname + ' ' + wdirOriginal + ' ' + wdir + ' ' + chain + ' ' + quantity + ' ' + str(pdbExists) + ' ' + directory +' '+ str(validdb))
except Exception,e :
    logging.error(" an error occurred while creating the log file | "+ LOG)
    logginfile = open(LOG, 'a')
    logginfile.write('END')
    logginfile.close()
    raise ValueError(" an error occurred while creating the log file | "+ LOG)
if quantity == '' or  directory == '' or wdirOriginal=='' or targetname=='':
    logging.warning(
        "A required option for the predictor is missing. \nRemember usage: python reupred.py -d <working path> -i <targetname> -c <chain> -q <all/first/chain> -e <True/False> -s <Softwarepath> -v <True/False>")
    logging.error("A required option for the predictor is missing. ")
    logginfile = open(LOG, 'a')
    logginfile.write('END')
    logginfile.close()
    raise ValueError("A required option for the predictor is missing. ")
try:
    if not os.path.isfile(directory + 'CONFIGFILE'):
        logging.error("The configuration file is missing" | directory + 'CONFIGFILE')
        logginfile = open(LOG, 'a')
        logginfile.write('END')
        logginfile.close()
        raise ValueError("The configuration file is missing    " | directory + 'CONFIGFILE')
    configFile = directory + 'CONFIGFILE'
except Exception, e:
    logging.error("An error occurred while trying to obtain the config file  | "+ directory + 'CONFIGFILE')
    logginfile = open(LOG, 'a')
    logginfile.write('END')
    logginfile.close()
    raise ValueError("An error occurred while trying to obtain the config file | "+ directory + 'CONFIGFILE')
try:
    if os.path.isfile(configFile):
        config = ConfigParser.ConfigParser()
        config.readfp(open(configFile))
        DBfolder = config.get('ReUPred variables', 'DBfolder')
        folderUnits = config.get('ReUPred variables', 'finalUnits')
        SRULdir = directory+config.get('ReUPred variables', 'SRULdir')
        AlignFolder = config.get('ReUPred variables', 'AlignFolder')
        templateList = config.get('ReUPred variables', 'templateList')
        reportFile = config.get('ReUPred variables', 'reportFile')
        tempFolder = config.get('ReUPred variables', 'tempFolder')
        musdir = config.get('ReUPred variables', 'musdir')
        pdbDesc = config.get('ReUPred variables', 'pdbDesc')
        pdbFolderDB = config.get('ReUPred variables', 'pdbFolderDB')
        pdbFolderDBlocal = config.get('ReUPred variables', 'pdbFolderDBlocal')
        TMalignexe   = config.get('ReUPred variables', 'TMalignexe')
        Mustangexe   = config.get('ReUPred variables', 'Mustangexe')
        Dsspexe  = config.get('ReUPred variables', 'Dsspexe')



except Exception,e:
    logging.error("there is a problem with the CONFIGFILE error , a variable is missing |" +  directory + 'CONFIGFILE')
    logginfile = open(LOG, 'a')
    logginfile.write('END')
    logginfile.close()
    raise ValueError("there is a problem with the CONFIGFILE error , a variable is missing |" +  directory + 'CONFIGFILE')
try:
    PDBdir = wdirOriginal
    PDBcompfolder = wdirOriginal
    AllOriginalTargets = []
    PDBregionDir = wdir + 'Region/'
    ALIdir = wdir + AlignFolder
    ALIdirReg = wdir + AlignFolder
    OUTdir = wdir + DBfolder
    OUTdirTemp = OUTdir
    if os.path.isdir(wdirOriginal + 'Temporary'):
        shutil.rmtree(wdirOriginal + 'Temporary', ignore_errors=True)
    create_all_directories(PDBregionDir, ALIdir, OUTdir)
    AllOriginalTargetschains = []
    OriginalTargetname = targetname
    successpdb = False

except Exception, e:
    logging.error("there was a mistake while trying to create folders ")
    logginfile = open(LOG, 'a')
    logginfile.write('END')
    logginfile.close()
    raise ValueError("there was a mistake while trying to create folders ")

if pdbExists == True or pdbExists == 'TRUE':
    try:

        print wdirOriginal + targetname
        if os.path.isfile(wdirOriginal + targetname):
            os.system('cp ' + wdirOriginal + targetname + ' ' + wdirOriginal + 'temp')
            os.system('cp ' + wdirOriginal + targetname + ' ' + wdirOriginal + 'reupred_'+targetname )


            files = [wdirOriginal + "header", wdirOriginal + 'reupred_' + OriginalTargetname]
            if os.path.isfile(wdirOriginal + "header") and os.path.isfile(wdirOriginal + 'reupred_' + OriginalTargetname):
                concat = ''.join([open(f).read() for f in files])
                if concat:

                    auxfile = open(wdirOriginal + 'reupred_' + OriginalTargetname, 'w')

                    auxfile.write(concat)
                    auxfile.close()
            targetname = 'temp'

            chainobtained, successpdb,pdb_header = download_pdb_chains(targetname, PDBdir, True, quantity, chain,pdbFolderDB,wdirOriginal)
        else:
            logging.error("xThere is no input file, include it or change predictor option -e to False ")
            shutil.rmtree(wdirOriginal + 'Temporary', ignore_errors=True)
            end = timeit.timeit()
            sys.exit(1)
    except Exception, e:
            logging.error(
                "An error occurred while trying to obtain pdb information from user file |" + OriginalTargetname)
            logging.warning("There is an error on your input file, probably not a pdb " + OriginalTargetname)
            logginfile = open(LOG, 'a')
            logginfile.write('END')
            logginfile.close()
            raise ValueError(
                "An error occurred while trying to obtain pdb information from user file |" + OriginalTargetname)
else:
       # print "calling ",targetname, PDBdir, False, quantity, chain
        chainobtained, successpdb, pdb_header = download_pdb_chains(targetname, PDBdir, False, quantity, chain,pdbFolderDB,wdirOriginal)


if successpdb:
    obtain_pdb_input(AllOriginalTargets, AllOriginalTargetschains, quantity, targetname, chainobtained, PDBdir, wdirOriginal)
    logging.debug("All pdb files obtained correctly")
    try:
        logging.debug("defining  variables and creating folders")
        itisarepeat = False
        templatesWork = tempFolder
        templatesWorkIII = tempFolder + "III"
        templatesWorkIV = tempFolder + "IV"
        templatesWorkV = tempFolder + "V"
        Raphaelperiod = {}
        AllOriginalTargetsx = []
        raphOption = True
        templateListRed = 'reducedList'
        data_directory = directory + 'Data/'
        Raphaelperiod = load_dictionary(Raphaelperiod,directory)
        with open(data_directory + templateListRed + 'V', 'r') as f:
                templatesListRedV = f.read().split("\n")
        with open(data_directory+templateListRed + 'IV', 'r') as f:
                templatesListRedIV = f.read().split("\n")
        with open(data_directory + templateListRed + 'III', 'r') as f:
                templatesListRedIII = f.read().split("\n")
        with open(data_directory + 'predictorValues.txt', 'r') as f:
            reupredvalues = f.read().split("\n")
        with open(data_directory + 'predictorValuesParticular.txt', 'r') as f:
            reupredvaluesparticular = f.read().split("\n")
        MSA, MSAtable= True, True
        DSSPdirectory = wdir + 'dssp/'
        DSSPdirectoryRegion = wdir + 'dsspRegion/'
        inSecondRound, targetCounter, i = 0, 0, 0
        OUTdirDB = wdir + "repeatsDB_RN/"
        tempName = tempFolder + 'AUX'
        tempNameIII = tempFolder + 'IIIAUX'
        tempNameIV = tempFolder + 'IVAUX'
        tempNameV = tempFolder + 'VAUX'
        TEMPLATESround, TEMPLATESroundIII, TEMPLATESroundIV, TEMPLATESroundV, mustangList = [], [], [], [], []
        newTargetsList, templatesTMdata, templatesTMdataLowerRound, templatesTMdataSecondRound = [], [], [], []
        templatesTMdataThirdRound, withoutprediction_alltargets = [], []
        origenUnitPath = wdir + templatesWork + 'units/'
        origenUnitPathIII = wdir + templatesWorkIII + 'units/'
        origenUnitPathIV = wdir + templatesWorkIV + 'units/'
        origenUnitPathV = wdir + templatesWorkV + 'units/'
        create_folders_tree(musdir, wdir, TEMPLATESround, TEMPLATESroundIII, TEMPLATESroundIV, TEMPLATESroundV, DSSPdirectory, OUTdirDB, tempName, PDBregionDir, ALIdir, OUTdir)
        Predictorpass, min_size_pdb = 0, 30000
        TargetsToEvaluate = len(AllOriginalTargets)
        id_classes = ["III", "IV", "V"]
        firstpass, region, moretoevaluate= True, False, False
        logging.debug("defined variables and created folders correctly")
    except  Exception, e:
        logging.error("Unable to define variables or create folders" + str(e))
        logginfile = open(LOG, 'a')
        logginfile.write('END')
        logginfile.close()
        raise ValueError("Unable to define variables or create folders" + str(e))
else:
    logging.error("Unable to obtain pdb files | " +OriginalTargetname)
    shutil.rmtree(wdirOriginal + 'Temporary', ignore_errors=True)
    logginfile = open(LOG, 'a')
    logginfile.write('END')
    logginfile.close()
    raise ValueError("Unable to obtain pdb files | " +OriginalTargetname)
try:
    if validdb.upper() == "TRUE" :
        validatedbfile=True
    else:
        validatedbfile = False
    while not validatedbfile and successpdb and targetCounter < TargetsToEvaluate:
        prevTargetCounter = targetCounter
        target = AllOriginalTargets[targetCounter]
        chain = AllOriginalTargetschains[targetCounter]
        print target
       # print Predictorpass
        logging.info("Evaluating "+target)
        moretoevaluate = True
        predicted = False
        strins = ''
        if os.stat(PDBdir+target+".pdb").st_size > min_size_pdb and moretoevaluate:
            if Predictorpass >= 1:
                logging.info("Seems that it might have more regions")

                os.system('mv ' + PDBregionDir + target + ".pdb " + PDBdir + target + ".pdb ")
                shutil.rmtree(wdir)
                create_folders_tree(musdir, wdir,  TEMPLATESround, TEMPLATESroundIII, TEMPLATESroundIV, TEMPLATESroundV, DSSPdirectory, OUTdirDB, tempName, PDBregionDir, ALIdir, OUTdir)
            else:
                create_folders_tree(musdir, wdir, TEMPLATESround, TEMPLATESroundIII, TEMPLATESroundIV, TEMPLATESroundV,
                                    DSSPdirectory, OUTdirDB, tempName, PDBregionDir, ALIdir, OUTdir)
                os.system('cp ' + PDBdir + target + ".pdb " + PDBdir + target + "Original.pdb ")
                logging.debug("Trying to predict from original target")
            moretoevaluate = False
            region = True
            templatesTMdataIII, templatesTMdataIV, templatesTMdataV, order = [], [], [], []
            logging.debug("Evaluating target against SRUL")
            existaligns = False
            if Predictorpass==0:
           #     print "target",target
                listorder, existaligns, templatesTMdataNameIII ,templatesTMdataNameIV, templatesTMdataNameV = identify_order_predictor(target, PDBdir, ALIdir, SRULdir, DSSPdirectory, templatesListRedIII, templatesListRedIV, templatesListRedV,Dsspexe,TMalignexe)
                os.system('rm '+wdirOriginal+'Temporary/Align/*')

                final_time_reduced = time()
                exec_time_reduced =   final_time_reduced -initial_time
                logging.info( 'Reduced list'+ str(listorder)+str(exec_time_reduced))
            if listorder!=[]:
                order, templatesTMdataIII, templatesTMdataIV, templatesTMdataV, existaligns = found_order_predictor(target, PDBdir, ALIdir, SRULdir, DSSPdirectory, listorder, templatesTMdataIII, templatesTMdataIV, templatesTMdataV, templatesTMdataNameIII ,templatesTMdataNameIV, templatesTMdataNameV, directory, templateList,Dsspexe,TMalignexe)
                templatesTMdata_list = [templatesTMdataIII, templatesTMdataIV, templatesTMdataV]
                final_time = time()
                exec_time =   final_time-initial_time
                logging.info('Reduced list'+ str(listorder)+ str(exec_time))
           # print order
            if not existaligns and not itisarepeat:
                create_output_non_predicted(wdirOriginal + '_' + chain + '/', OriginalTargetname+'_'+chain, Raphaelperiod)
                os.system('cp '+wdirOriginal + '_' + chain + '/', OriginalTargetname+'_'+chain+' '+wdirOriginal + '_' + chain + '/', OriginalTargetname+'_'+chain+'.reupred')
                pdb_header['chain ' + chain] = "Predicted as Non Repeated"
                shutil.rmtree(wdirOriginal + 'Temporary', ignore_errors=True)
                if os.path.isfile(wdirOriginal + 'temp*'):
                    os.system('rm ' + wdirOriginal + 'temp*')
                logging.info(OriginalTargetname + " chain: " + chain + " is not a repeat")

            else:
                logging.info("Starting the prediction process.......")
                for sc in listorder:
                    logging.debug(" trying Class "+str(sc))
                    if not predicted and templatesTMdata_list[sc] != []:
                        predicted = Switcher().call_predictor(wdir, sc, target, templatesTMdata_list[sc], PDBdir, OUTdir, OUTdirDB, ALIdir, SRULdir, directory, templatesWorkIII, templatesWorkIV, templatesWorkV, TEMPLATESroundIII, TEMPLATESroundIV, TEMPLATESroundV, DSSPdirectory, OUTdir, Raphaelperiod, reupredvalues, reupredvaluesparticular, origenUnitPathIII, origenUnitPathIV, origenUnitPathV, region,TMalignexe,Mustangexe,Dsspexe)
                        #logging.debug("Result:",predicted, moretoevaluate)
                        if predicted:
                            logging.debug("in more to evaluate")
                            moretoevaluate = verify_region(target, OUTdir + id_classes[sc] + "/", PDBdir, PDBregionDir)
                            os.system('cp ' + OUTdirDB + target + '.db ' + wdirOriginal + '_' + chain + '/region' + str(Predictorpass))
                            logging.debug('cp ' + OUTdirDB + target + '.db ' + wdirOriginal + '_' + chain + '/region' + str(
                                Predictorpass))
                            itisarepeat = True

                if moretoevaluate:
                    targetCounter -= 1
                    Predictorpass += 1
                    predicted = False
                    logging.debug("analysing the rest of the protein")
                else:

                    if os.path.isfile(wdirOriginal + '_' + chain + '/region0'):
                        move_files_finaldestination(target, wdirOriginal,  chain, OriginalTargetname)
                        logging.info(
                            "Finished the prediction process..... chain: " + chain + ' ' + str(Predictorpass+1) + ' regions')
                       # print "Finished the prediction process"
                        itisarepeat= reevaluate_regions(wdirOriginal + '_' + chain + '/', target, OriginalTargetname + '_' + chain, Predictorpass, Raphaelperiod, chain,Mustangexe,TMalignexe,Dsspexe,wdirOriginal)
                        Predictorpass == 0
                        #print itisarepeat
                        shutil.rmtree(wdirOriginal + 'Temporary', ignore_errors=True)
                        if not itisarepeat:
                            logging.info("ReUPred does not identifies it as a repeat")

                            create_output_non_predicted(wdirOriginal + '_' + chain + '/',
                                                        OriginalTargetname + '_' + chain, Raphaelperiod)
                            logging.info("ReUPred does not identifies it as a repeat")
                            pdb_header['chain ' + chain] = "Predicted as Non Repeated"
                        else:
                            pdb_header['chain '+chain] = "Predicted as Repeated"
                            logging.info("Finished the prediction process..... chain: " + chain + ' ' + str(Predictorpass+1) + ' regions')
                        Predictorpass = 0
                        shutil.rmtree(wdirOriginal + 'Temporary', ignore_errors=True)
                    else:#if not itisarepeat:
                        logging.info("ReUPred does not identifies it as a repeat")
                        pdb_header['chain ' + chain] = "Predicted as Non Repeated"

                        create_output_non_predicted(wdirOriginal + '_' + chain + '/', OriginalTargetname+'_'+chain,Raphaelperiod)
                        Predictorpass == 0
                        shutil.rmtree(wdirOriginal + 'Temporary', ignore_errors=True)

        elif os.stat(PDBdir + target + ".pdb").st_size <= min_size_pdb or not itisarepeat:
                logging.info("ReUPred does not identifies it as a repeat")
                create_output_non_predicted(wdirOriginal + '_' + chain + '/', OriginalTargetname+'_'+chain, Raphaelperiod)
                pdb_header['chain ' + chain] = "Predicted as Non Repeated"
                Predictorpass == 0

                shutil.rmtree(wdirOriginal + 'Temporary', ignore_errors=True)
        targetCounter += 1

except Exception, e:
    logging.error("There was an error during the predicting process")
    logginfile = open(LOG, 'a')
    logginfile.write('END')
    logginfile.close()
    raise ValueError("There was an error during the predicting process")
try:
    if validatedbfile:
        while validatedbfile and successpdb and targetCounter < TargetsToEvaluate:
            try:
                prevTargetCounter = targetCounter
                target = AllOriginalTargets[targetCounter]
                chain = AllOriginalTargetschains[targetCounter]
                if os.path.isfile(wdirOriginal + '_' + chain + '/'+ OriginalTargetname + '_' + chain + '.db'):
                    if not os.path.isfile(wdirOriginal + '_' + chain + '/'+OriginalTargetname + '_' + chain+ '.pdb '):
                        if os.path.isfile( wdirOriginal + OriginalTargetname + chain + '.pdb'):
                            os.system('cp ' + wdirOriginal + OriginalTargetname + chain+ '.pdb ' + wdirOriginal + '_' + chain + '/'+OriginalTargetname + '_' + chain+ '.pdb ')
                        elif os.path.isfile(wdirOriginal + 'temp' + chain + '.pdb'):
                                os.system(
                                    'cp ' + wdirOriginal + 'temp' + chain + '.pdb '+ wdirOriginal + '_' + chain + '/' + OriginalTargetname + '_' + chain + '.pdb ')
                    logging.info("Evaluating " + target)

                    for foldername in os.listdir(wdirOriginal + '_' + chain):
                        if os.path.isdir(wdirOriginal + '_' + chain + '/' + foldername):
                            shutil.rmtree(wdirOriginal + '_' + chain + '/' + foldername)
                    reevaluate_regions_validation(wdirOriginal + '_' + chain + '/', OriginalTargetname + '_' + chain, chain,Mustangexe,TMalignexe,Dsspexe,wdirOriginal)
                else:
                    logging.warning("There is no db file for: " + OriginalTargetname + chain + '|' +  + OriginalTargetname + chain + '.db')
                    logginfile = open(LOG, 'a')
                    logginfile.write('END')
                    logginfile.close()
                    raise ValueError("A problem occurred during the validation process of " + OriginalTargetname + chain )
            except Exception, e:
                logging.error("A problem occurred during the validation process of " + OriginalTargetname + chain + str(e))
                logginfile = open(LOG, 'a')
                logginfile.write('END')
                logginfile.close()
                raise ValueError("A problem occurred during the validation process of " + OriginalTargetname + chain + str(e))
            targetCounter += 1
except Exception, e:
    logging.error("there was an error during the validation process | " + wdirOriginal + '_' + chain + '/' + OriginalTargetname + '_' + chain + '.db')
    logging.warning("there was an error during the validation process the file " + wdirOriginal + '_' + chain + '/' + OriginalTargetname + '_' + chain  + ".db does not exist")
    logginfile = open(LOG, 'a')
    logginfile.write('END')
    logginfile.close()
    raise ValueError("there was an error during the validation process | " + wdirOriginal + '_' + chain + '/' + OriginalTargetname + '_' + chain + '.db')


try:
    create_pdb_description_file(wdirOriginal ,pdb_header,OriginalTargetname)
    if os.path.isdir(wdirOriginal+'Temporary'):
        shutil.rmtree(wdirOriginal+'Temporary', ignore_errors=True)
    os.system('rm ' + wdirOriginal + '*.pdb')

    if os.path.isfile(wdirOriginal + 'temp'):
       os.system('rm '+wdirOriginal + 'temp*')

    if not bool(pdbExists) or pdbExists=='FALSE':
        os.system('cp ' + wdirOriginal + 'pdb'+OriginalTargetname + '.ent ' + wdirOriginal + 'reupred_' + OriginalTargetname)
        files = [wdirOriginal + "header", wdirOriginal + 'reupred_' + OriginalTargetname]
        concat = ''.join([open(f).read() for f in files])
        auxfile = open(wdirOriginal + 'reupred_' + OriginalTargetname, 'w')
        auxfile.write(concat)
        auxfile.close()
except Exception, e :
    logging.error("it was not possible to erase the temporary folder |" + wdirOriginal+"Temporary")
    logginfile = open(LOG, 'a')
    logginfile.write('END')
    logginfile.close()
    raise ValueError("it was not possible to erase the temporary folder |" + wdirOriginal + "Temporary")

logginfile=open(LOG,'a')
logginfile.write('END')
logginfile.close()

