from otherMethods import *
import os
import subprocess
import re
from Bio.PDB import *
from operator import itemgetter
import time
import logging
import argparse
from Bio.PDB import *
from alignmentManager import ProteinAlignmentSet
from operator import itemgetter


class Switcher(object):
    def call_predictor(self, wdir, argument, target, templatesTMdataOrig, PDBdir, OUTdir, OUTdirDB, ALIdir, SRULdir, directory, templatesWorkIII, templatesWorkIV, templatesWorkV, TEMPLATESroundIII, TEMPLATESroundIV, TEMPLATESroundV, DSSPdirectory, OUTdirTemp, Raphaelperiod, reupredvalues, reupredvaluesparticular, origenUnitPathIII, origenUnitPathIV, origenUnitPathV, regionF,TMalignexe,Mustangexe,Dsspexe):
        """Dispatch method"""
        # prefix the method_name with 'number_' because method names
        # cannot begin with an integer.
        method_name = 'call_class' + str(argument+3)+'_predictor'
        # Get the method from 'self'. Default to a lambda.
        method = getattr(self, method_name, lambda: "nothing")
        # Call the method as we return it
        logging.debug("Calling the predictor for class"+str(argument+3))
        return method(wdir, target, templatesTMdataOrig, PDBdir, OUTdir, OUTdirDB, ALIdir, SRULdir, directory, templatesWorkIII, templatesWorkIV, templatesWorkV, TEMPLATESroundIII, TEMPLATESroundIV,TEMPLATESroundV, DSSPdirectory, OUTdirTemp, Raphaelperiod, reupredvalues, reupredvaluesparticular, origenUnitPathIII, origenUnitPathIV, origenUnitPathV, regionF,TMalignexe,Mustangexe,Dsspexe)

    def call_class5_predictor(self, wdir, target, templatesTMdataOrig, PDBdir, OUTdir, OUTdirDB, ALIdir, SRULdir, directory, templatesWorkIII, templatesWorkIV, templatesWorkV, TEMPLATESroundIII, TEMPLATESroundIV, TEMPLATESroundV,  DSSPdirectory, OUTdirTemp, Raphaelperiod, reupredvalues, reupredvaluesparticular, origenUnitPathIII, origenUnitPathIV, origenUnitPathV, regionF,TMalignexe,Mustangexe,Dsspexe):
        if regionF:
            fold = "V/"
        else:
            fold = "/"
        inSecondRound = 0
        logging.debug("Initialized predictor object class 5 ")
        ReUPred = Predictor(wdir, reupredvalues[2], reupredvaluesparticular[2], SRULdir, directory, templatesWorkV, target, PDBdir, OUTdir + fold, OUTdirDB, OUTdirTemp, ALIdir, TEMPLATESroundV, DSSPdirectory, Raphaelperiod, templatesTMdataOrig, origenUnitPathV ,TMalignexe,Mustangexe,Dsspexe)
        logging.debug("Running repeat predictor V")
        predicted, infoValidate = ReUPred.repeat_predictor("V", inSecondRound, regionF)
        return predicted

    def call_class4_predictor(self, wdir, target, templatesTMdataOrig, PDBdir, OUTdir, OUTdirDB, ALIdir, SRULdir, directory, templatesWorkIII, templatesWorkIV, templatesWorkV, TEMPLATESroundIII, TEMPLATESroundIV, TEMPLATESroundV,  DSSPdirectory, OUTdirTemp, Raphaelperiod, reupredvalues, reupredvaluesparticular, origenUnitPathIII, origenUnitPathIV, origenUnitPathV, regionF,TMalignexe,Mustangexe,Dsspexe):
        if regionF:
            fold = "IV/"
        else:
            fold = "/"
        inSecondRound = 0
        logging.debug("Initialized predictor object class 4")
        ReUPred = Predictor(wdir, reupredvalues[1], reupredvaluesparticular[1], SRULdir, directory, templatesWorkIV, target, PDBdir, OUTdir + fold, OUTdirDB, OUTdirTemp, ALIdir, TEMPLATESroundIV, DSSPdirectory, Raphaelperiod, templatesTMdataOrig, origenUnitPathIV,TMalignexe,Mustangexe,Dsspexe)
        logging.debug("Running repeat predictor IV")
        logging.debug("path to align "+str(TMalignexe))
        predicted, infoValidate = ReUPred.repeat_predictor("IV", inSecondRound, regionF)
        return predicted

    def call_class3_predictor(self, wdir, target, templatesTMdataOrig, PDBdir, OUTdir, OUTdirDB, ALIdir, SRULdir, directory, templatesWorkIII, templatesWorkIV, templatesWorkV, TEMPLATESroundIII, TEMPLATESroundIV,TEMPLATESroundV,  DSSPdirectory, OUTdirTemp, Raphaelperiod, reupredvalues, reupredvaluesparticular, origenUnitPathIII, origenUnitPathIV, origenUnitPathV, regionF,TMalignexe,Mustangexe,Dsspexe):
        if regionF:
            fold = "III/"
        else:
            fold = "/"
        inSecondRound = 0
        logging.debug("Initialized predictor object class 3")
        ReUPred = Predictor(wdir, reupredvalues[0], reupredvaluesparticular[0], SRULdir, directory, templatesWorkIII, target, PDBdir, OUTdir + fold, OUTdirDB, OUTdirTemp, ALIdir, TEMPLATESroundIII,DSSPdirectory, Raphaelperiod, templatesTMdataOrig, origenUnitPathIII,TMalignexe,Mustangexe,Dsspexe )
        logging.debug("Running  repeat predictor III")
        predicted, infoValidate = ReUPred.repeat_predictor("III", inSecondRound, regionF)
        return predicted


class Predictor:




    def __init__(self, wdir, row, rowparticular, SRULdir, directory, templatesWork ,target, PDBdir, OUTdir, OUTdirDB, OUTdirTemp, ALIdir, TEMPLATESround, DSSPdirectory, Raphaelperiod, templatesTMdataOrig, origenUnitPath ,TMalignexe,Mustangexe,Dsspexe):
        self.TMalignexe, self.Mustangexe, self.Dsspexe= TMalignexe, Mustangexe, Dsspexe
        particulardata = rowparticular.split(",")
        self.PorclenUnitTarget1, self.PorclenUnitTarget2, self.LimlenUnitTarget1, self.LimlenUnitTarget2, \
        self.LimTMscore1, self.LimTMscore2, self.LimgapUnit1,self.Limgap1, self.LimRMSD1, self.LimgapUnit2, \
        self.Limgap2, self.LimRMSD2, self.porcMinLenUnit,  self.LimTMscore3, self.LimTMscore4, self.PorclenUnitTarget3,\
        self.PorclenTarget4, self.PorclenUnit4, self.LimlenUnitTarget3, self.LimlenUnitTarget4, self.LimgapUnit3, \
        self.Limgap3, self.LimRMSD3, self.LimgapUnit4, self.Limgap4, self.LimRMSD4, self.MinTMscoreInsertion, \
        self.procLenInsertion, self.minLenTemplate,self.coverageminfirstunit1,self.coverageminfirstunit2, \
        self.coverageminfirstunit3 = map(float, row.split(","))
        self.minTMscore, self.minalilen, self.porclenminPart, self.minCoverage,\
            self.porcLenPart1, self.tmscorePart1, self.porcLenUnitPart1, self.porcLenTargatPart1, self.coveragePart1, self.maxlen1,\
            self.porcLenPart2, self.tmscorePart2, self.porcLenUnitPart2, self.porcLenTargatPart2, self.coveragePart2, self.maxlen2,\
            self.porcLenPart3, self.tmscorePart3, self.coveragePart3, self.maxlen3,\
            self.porcLenPart4, self.tmscorePart4, self.coveragePart4, self.maxlen4 \
            =particulardata
        self.wdir, self.SRULdir, self.directory, self.templatesWork, self.target, self.PDBdir, self.OUTdir, self.OUTdirDB, \
        self.OUTdirTemp, self.alignName, self.TEMPLATESround , self.DSSPdirectory, \
         self.Raphaelperiod, self.templatesTMdataOrig,  \
         =  wdir, SRULdir, directory, templatesWork, target, PDBdir, OUTdir, OUTdirDB, OUTdirTemp, ALIdir, \
            TEMPLATESround, DSSPdirectory, Raphaelperiod,\
                              templatesTMdataOrig
        self.origenUnitPath = origenUnitPath
        logging.debug("All predictor variables initialized")

    def repeat_predictor(self,sclass ,inSecondRound,regionsSet):
        wdir, PDBdir, OUTdir, ALIdir, SRULdir, directory, templatesWork, TEMPLATESround, DSSPdirectory, OUTdirDB, OUTdirTemp, \
            target, templatesTMdataOrig, Raphaelperiod, origenUnitPath = self.wdir, self.PDBdir, \
         self.OUTdir, self.alignName, self.SRULdir, self.directory, self.templatesWork, self.TEMPLATESround, \
        self.DSSPdirectory, self.OUTdirDB, self.OUTdirTemp, self.target, self.templatesTMdataOrig,    \
        self.Raphaelperiod, self.origenUnitPath
        TMalignexe, Mustangexe, Dsspexe =self.TMalignexe,self.Mustangexe,self.Dsspexe
        enoughtemplates, difreg = True, True
        fastpListAA, fastpList, mustangList = [], [], []
        templateListOpt, targetsListOpt, rangeMINMAX, UsedTargetsList = [], [], [], []
        targetOriginal = target
        create_folders_target(wdir, templatesWork, target, PDBdir)
        while inSecondRound < 4 and enoughtemplates:
            if inSecondRound == 0 and os.path.isfile(PDBdir+target[:5]+'.pdb'):

                TEMPLATESroundInfo = []
                temp = inSecondRound
                infoValidate, strins = "", ""
                mustangFile = open(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.mus', 'w')
                mustangFile.write(">" + self.OUTdir + "\n")
                mustangFile.close()
                templatesTMdata, templatesTMdataSecondRound, templatesTMdataThirdRound, templatesTMdataLowerRound = [], [], [], []
                topTemplate, secontTemp, ThirdTemp = "", "", ""
                tmscore = 0.0
                templatesTMdata = self.add_data_to_templates_tmdata(templatesTMdata,sclass)
                logging.debug(templatesTMdata)
                i, it = 0, 1
                if len(templatesTMdata) > 0:
                    self.create_templates_list(templatesTMdata)

                templateListOpt, targetsListOpt, rangeMINMAX = [], [], []
                # calculates the first unit
                if len(templatesTMdataSecondRound) > 0:
                    templatesTMdataSecondRound.sort(cmp=None)
                    self.create_templates_list(templatesTMdataSecondRound)

                templatesTMdataThirdRound, templatesTMdataLowerRound = self.add_data_to_templates_tmdata_bigger_gaps(templatesTMdataThirdRound, templatesTMdataLowerRound)
                templatesTMdataThirdRound.sort(cmp=None)
                if len(templatesTMdataThirdRound) > 0:
                    self.create_templates_list(templatesTMdataThirdRound)

                if len(templatesTMdataLowerRound) > 0:
                    templatesTMdataLowerRound.sort(cmp=None)
                    self.create_templates_list(templatesTMdataLowerRound)

                templatesTMdata = self.put_all_templates_tmdata(templatesTMdata, templatesTMdataSecondRound, templatesTMdataThirdRound, templatesTMdataLowerRound, False)

                tmpcount = 0
                withoutprediction = []
                while tmpcount > 1 and tmpcount < 10 and len(templatesTMdata) > tmpcount:
                    withoutprediction.append(templatesTMdata[tmpcount][0])
                    tmpcount += 1


            if sclass == 'IV' or sclass == 'III':
                maxLen = 100
            else:
                maxLen = 100
            logging.debug(PDBdir + target[:5] + '.pdb')

            if os.path.isfile(PDBdir + target[:5] + '.pdb'):
                rangeMINMAX = []
                if len(templatesTMdata) > 0:
                    while rangeMINMAX == [] or int(rangeMINMAX[1]) - int(rangeMINMAX[0]) > maxLen and len(templatesTMdata)>0:
                        rangeMINMAX = []

                        self.calc_start_end_struct_alignment(templatesTMdata, inSecondRound, rangeMINMAX)

                    fastpList, InfoRanges = [], []
                    if rangeMINMAX != []:
                        pdbFile = open(wdir + templatesWork + target[:5] + '/' + target + '.pdb', 'r')
                        parser = PDBParser()
                        structure = parser.get_structure(target[:4], pdbFile)
                        pdbFile.close()

                        self.create_pdb_fragments(rangeMINMAX, i+1, structure, 1, targetsListOpt, templateListOpt, fastpList, target, targetOriginal, InfoRanges)
                        topTemplateActual = []
                        if len(templatesTMdata)>0 and i<len(templatesTMdata):
                            topTemplateActual.append(templatesTMdata[i][0])

                UsedTargetsList = []

                if len(templatesTMdata) > 0:

                    logging.debug("evaluating other units ")
                    lenFirstUnit = int(templatesTMdata[inSecondRound][7])
                    found = -1
                    oldtargetsListOpt, listEval, templatesTMdata2 = [], [], []
                    antlenTargetValue, antlenTemplatesValue, countTemp, lenTotrepeat, it, ConsitionSetPart = 0, 0, 0, 1, 1, 1
                    while( antlenTemplatesValue != len(templateListOpt) or (oldtargetsListOpt != targetsListOpt)) and targetsListOpt != []:
                        antlenTemplatesValue = len(templateListOpt)
                        oldtargetsListOpt = targetsListOpt
                        targetsListOpt = sorted(targetsListOpt)
                        antlenTargetValue = len(targetsListOpt)

                        for targetInOption in targetsListOpt:
                            if len(targetInOption) > 0:
                                for unit in templateListOpt:
                                    if len(unit) > 0 and os.path.isfile(wdir + templatesWork + targetInOption[:5] + '/' + targetInOption + '.pdb'):
                                        try: #aca

                                            tmOutput = self.execute_tmalign(targetInOption, unit, '.pdb')

                                            works = []
                                            ConsitionSet = 1
                                            templatesTMdata2 = self.add_data_to_templates_tmdata_particular(tmOutput, templatesTMdata2, fastpList, unit, temp, works, targetInOption, i, ConsitionSet, lenFirstUnit,sclass,Dsspexe)
                                            if len(templatesTMdata2) <= 2:
                                                ConsitionSet = 2
                                                templatesTMdata2 = self.add_data_to_templates_tmdata_particular(tmOutput, templatesTMdata2, fastpList, unit, temp, works, targetInOption, i, ConsitionSet, lenFirstUnit,sclass,Dsspexe)
                                            if len(templatesTMdata2) <= 2:
                                                ConsitionSet = 3
                                                templatesTMdata2 = self.add_data_to_templates_tmdata_particular(tmOutput, templatesTMdata2, fastpList, unit, temp, works, targetInOption, i, ConsitionSet,  lenFirstUnit,sclass,Dsspexe)
                                            if len(templatesTMdata2) <= 2:
                                                ConsitionSet = 4
                                                templatesTMdata2 = self.add_data_to_templates_tmdata_particular(tmOutput, templatesTMdata2, fastpList, unit, temp, works, targetInOption, i, ConsitionSet,  lenFirstUnit,sclass,Dsspexe)
                                            temp += 1
                                        except subprocess.CalledProcessError as grepexc:
                                            logging.error( "error code"+ grepexc.returncode+ grepexc.output)

                        if len(templatesTMdata2) > 0:
                            templatesTMdata2.sort(cmp=None)
                            templatesTMdata2.sort(key=itemgetter(1), reverse=True)# tmscore
                            templatesTMdata2.sort(key=itemgetter(7), reverse=True)# len
                            templatesTMdata2.sort(key=itemgetter(10))# rmsd
                            templatesTMdata2 = self.put_all_templates_tmdata_particular(templatesTMdata2, templatesTMdata2, templatesTMdata2, templatesTMdata2, True)

                            it += 1
                            rangeMINMAX = []
                            worked = False
                            targetInOption = templatesTMdata2[0][11]
                            #print "rest units",templatesTMdata2
                            while templatesTMdata2 != [] and not(worked):
                                    self.calc_start_end_struct_alignment_particular(templatesTMdata2, i, targetInOption, rangeMINMAX,TMalignexe)
                                    logging.debug(wdir + templatesWork + targetInOption[:5] + '/' + targetInOption + '.pdb '+str(rangeMINMAX))
                                    pdbFileOpt2 = open(wdir + templatesWork + targetInOption[:5] + '/' + targetInOption + '.pdb', 'r')
                                    parserOpt = PDBParser()
                                    if rangeMINMAX != []:
                                        structureOpt = parserOpt.get_structure(targetInOption[:4], pdbFileOpt2)
                                        pdbFileOpt2.close()
                                        worked, templateListOpt = self.create_pdb_fragments_particular(rangeMINMAX, i + 1, structureOpt, it, targetsListOpt, templateListOpt, fastpList, targetInOption, UsedTargetsList, InfoRanges,  targetOriginal, worked)

                                    if not worked:
                                        templatesTMdata2.remove(templatesTMdata2[0])
                                        worked = False

                                    lenTemplatesValue = len(templateListOpt)
                                    logging.debug("End of while")
                            countTemp = 0
                        templatesTMdata2 = []
                        targetsListOpt = list(set(targetsListOpt) - set(UsedTargetsList))
                        targetsListOpt = sorted(targetsListOpt)
                    logging.debug("after all units")
                    contFastp, w, totLostRes, count, ConsitionSetPart = 0, 0, 0, 0, 0
                    fastpListAA,fastpList = [], []
                    if os.path.isfile(wdir+templatesWork + target[:5] + '/' + target[:5] + '.fastp'):
                        with open(wdir+templatesWork + target[:5] + '/' + target[:5] + '.fastp', 'r') as flis:
                            fastpListAA = flis.read().split("\n")

                    fastpFinalLen = len(fastpListAA) - 1
                    numreg = 0
                    contFastp = 0
                    validate = 0
                    if len(fastpListAA) > 1:
                        origenUnitPath = wdir + templatesWork + 'units/'
                        origenPath = wdir + templatesWork + target[:5] + '/'
                        fastpListAA,strins = self.verify_holes(target, fastpListAA, sclass, strins)
                        logging.debug("after verify holes")
                        fffas = open(wdir + templatesWork + target[:5]+'/' + target[:5] + '.fastp', 'w')
                        for rp in fastpListAA:
                            fffas.write(str(rp)+"\n")
                        fffas.close()
                        actholes = self.calc_holes(fastpListAA)
                        logging.debug("after calc holes")
                        validate = self.contiguous_set(fastpListAA)
                        logging.debug("after contiguos")
                        contFastp = validate * (len(fastpListAA) - 2)
                        destinationPath = self.OUTdir

                        fastpListAA, region ,insertions,ninsertions= sort_fastp_result(fastpListAA)
                        logging.debug("after sort")
                        numreg = region.count(";")

                        mustangOK, changed, fastpListAA = evaluatemustangResult(fastpListAA, wdir + templatesWork + target[:5]+'/', targetOriginal, target[4],self.Mustangexe)
                        logging.debug("after mustang")
                        #print mustangOK, changed, fastpListAA
                        if mustangOK and changed:
                            actholes = self.calc_holes(fastpListAA)
                            validate = self.contiguous_set(fastpListAA)
                            contFastp = validate * (len(fastpListAA) - 2)

                            fastpListAA, region ,insertions,ninsertions= sort_fastp_result(fastpListAA)

                            fffas = open(wdir + templatesWork + target[:5] + '/' + target[:5] + '.fastp', 'w')
                            for rp in fastpListAA:
                                fffas.write(str(rp) + "\n")
                            fffas.close()
                            actholes=actholes-ninsertions
                    NumUnits=3
                    if mustangOK and contFastp >= 1 and (len(fastpListAA) >= NumUnits and (validate >= float(0.6) and actholes <= 18) or (len(fastpListAA) >= NumUnits and validate >= float(0.57) and actholes <= 45)or (validate >= float(0.75) and len(fastpListAA) >= NumUnits)):
                        #if not changed:
                        logging.debug("after conditions to result")
                        fastpListAA2 = avoid_holes(fastpListAA, wdir + templatesWork + target[:5] + '/' + target[:5],self.Dsspexe)
                        logging.debug("after avoid holes")
                        fastpListAA, region,insertions,ninsertions = sort_fastp_result(fastpListAA2)
                        logging.debug("after sort")
                        print "Start\tEnd\n"
                        strUnits = ''
                        for row in fastpListAA:
                            if len(row) > 0:
                                print row.split("\t")[0], "\t", row.split("\t")[1]
                                strUnits = strUnits + ";" + row.split("\t")[0] + " " + row.split("\t")[1]
                        print "Selected template", templatesTMdata[inSecondRound][0][:5]
                        selectedTemplate = templatesTMdata[inSecondRound][0]
                        mustangList.append(OUTdir + target[:5] + '.mus')
                        logging.debug("befor identify class ")

                        if sclass == "V":
                            txt, sc = identify_classification_class_v(target, directory, templatesTMdata, inSecondRound, strUnits, region, False, strins)
                        elif sclass == "IV":
                            txt,  sc = identify_classification_class_iv(target, directory, templatesTMdata, inSecondRound, strUnits, region, False, strins)
                        else:
                            txt,  sc = identify_classification_class_iii( target, directory,   templatesTMdata, inSecondRound, strUnits, region, False, strins)
                        logging.debug("after identify"+ txt+" "+sc)
                        create_final_output(target, OUTdirDB, region, fastpListAA, Raphaelperiod, insertions, sc, regionsSet,selectedTemplate)
                        logging.debug("after create final output")
                        aliList, templatesTMdata, templatesTMdataSecondRound, templatesTMdataThirdRound, templatesTMdataLowerRound = [], [], [], [], []
                        inSecondRound, rangeMinTopTemplate, rangeMinSecTemplate, rangeMinThirdTemplate, rangeMaxTopTemplate, rangeMaxSecTemplate, rangeMaxThirdTemplate = 0, 0, 0, 0, 0, 0, 0
                        move_files(origenUnitPath, target, destinationPath, origenPath)
                        logging.debug("after move files")
                        if os.path.isfile(destinationPath + target[:5] + '*_*'):
                            os.remove(destinationPath + target[:5] + '*_*')
                        logging.debug("after Result "+infoValidate)
                        return True, infoValidate
                    else:
                        #print "en el else"
                        mustangFile = open(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.mus', 'w')
                        mustangFile.write(">" + self.OUTdir + "\n")
                        mustangFile.close()
                        if len(fastpListAA) > 3:
                            fffas = open(wdir + templatesWork + target[:5]+'/' + target[:5] + '.fastp', 'w')
                            for rp in fastpListAA:
                                fffas.write(str(rp) + "\n")
                            fffas.close()
                            actholes = self.calc_holes(fastpListAA)
                            validate = self.contiguous_set(fastpListAA)
                            TEMPLATESroundInfo.append([inSecondRound, actholes, validate, templatesTMdata[inSecondRound][0][:5], fastpListAA, strins])
                        if inSecondRound == 3:
                            if os.path.isfile(origenUnitPath + target[:5] + '*_* '):
                                os.remove(origenUnitPath + target[:5] + '*_* ')
                            if len(fastpListAA) > 2:
                                destinationPath = TEMPLATESround[2] + '/'
                                move_files(origenUnitPath, target, destinationPath, origenPath)
                            if difreg:
                                infoValidate = self.contiguous_set_all(target[:5], TEMPLATESround, OUTdirTemp, TEMPLATESroundInfo, templatesTMdata, sclass, strins)
                            validate = False
                            #print "validate False",infoValidate
                            inSecondRound = 0
                            it = 1
                            rangeMINMAX, targetsListOpt, templateListOpt, templatesTMdata, templatesTMdataSecondRound, templatesTMdataThirdRound, templatesTMdataLowerRound = [], [], [], [], [], [], []
                            return validate, infoValidate
                        else:
                            mustangFile = open(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.mus', 'w')
                            mustangFile.write(">" + self.OUTdir + "\n")
                            mustangFile.close()
                            if len(fastpListAA) > 2:
                                if inSecondRound == 0:
                                    destinationPath = TEMPLATESround[3] + '/'
                                if inSecondRound == 1:
                                    destinationPath = TEMPLATESround[0] + '/'
                                if inSecondRound == 2:
                                    destinationPath = TEMPLATESround[1] + '/'
                                if os.path.isfile(origenUnitPath + target[:5] + '*_* '):
                                    os.remove(origenUnitPath + target[:5] + '*_* ')
                                move_files(origenUnitPath, target, destinationPath, origenPath)
                            it = 1
                            strins = ""
                            templateListOpt, targetsListOpt, rangeMINMAX = [], [], []
                            inSecondRound += 1
                            temp = inSecondRound * 1000
                            if (len(templatesTMdata)-1) < inSecondRound:
                                inSecondRound = 0
                                aliList = []
                                it = 1
                                rangeMINMAX, targetsListOpt, templateListOpt, templatesTMdata, templatesTMdataSecondRound, templatesTMdataThirdRound, templatesTMdataLowerRound = [], [], [], [], [], [], []
                                return False, infoValidate
                else:

                    inSecondRound = 0
                    aliList = []
                    return False,infoValidate

    def add_data_to_templates_tmdata(self, templatesTMdata,sc):
        wdir, target = self.wdir,  self.target
        if sc =='III':
            MinLenSubclass = 22
        elif sc == 'IV':
            MinLenSubclass = 26
        else:
            MinLenSubclass = 31

        for line in self.templatesTMdataOrig:
            if line:
                row = list(line.get_info_align())

                if row and row[0] != ''  :
                    contig = 0
                    unit, lenTarget, start, end, nicePosition, lenUnit, alignmentLength, RmsdValue, sequenceIdentity, tmScoreO, tmscoreAct, tmScoreAVG, align, sequenceRelation, align3, gapsPerc, gapsUnitPerc, coverage = (row[1]), int(row[2]), int(row[3]), int(row[4]), int(row[5]), int(row[6]), float(row[7]), float(row[8]), float(row[9]), float(row[10]), float(row[11]), float(row[12]), row[13], row[14], row[15], float(row[16]), float(row[17]), float(row[18])
                    if float(RmsdValue) < float(self.LimRMSD2):
                        if float(tmscoreAct) >= float(self.LimTMscore2):
                            if float(self.LimlenUnitTarget2) > float(lenUnit) >= float(self.minLenTemplate)*float(self.PorclenUnitTarget2) and float(self.LimlenUnitTarget2) > lenTarget >= float(self.minLenTemplate)*float(self.PorclenUnitTarget2) and float(coverage) > float(self.coverageminfirstunit1)  and end - start >= MinLenSubclass:

                                templatesTMdata.append([unit, round(tmscoreAct, 2), start, end, gapsPerc, gapsUnitPerc, contig, lenUnit, nicePosition, RmsdValue, len(align[:start])-align[:start].count('-'), coverage])

        return templatesTMdata

    def add_data_to_templates_tmdata_bigger_gaps(self, templatesTMdataThirdRound, templatesTMdataLowerRound):
        wdir, SRULdir, templatesWork, target, PDBdir, targetInOption = self.wdir, self.SRULdir, self.templatesWork, self.target, self.PDBdir, self.target
        temp = 1
        for line in self.templatesTMdataOrig:
            if line:
                row = list(line.get_info_align())
                if row and row[0] != '':
                    contig = 0
                    target, unit, lenTarget, start, end, nicePosition, lenUnit, alignmentLength, RmsdValue, sequenceIdentity, tmScoreO, tmscoreAct, tmScoreAVG, align, sequenceRelation, align3, gapsPerc, gapsUnitPerc, coverage = row[0], (row[1]), int(row[2]), int(row[3]), int(row[4]), int(row[5]), int(row[6]), float(row[7]), float(row[8]), float(row[9]), float(row[10]), float(row[11]), float(row[12]), row[13], row[14], row[15], float(row[16]), float(row[17]), float(row[18])
                    if float(RmsdValue)<self.LimRMSD4 and unit[:5] != target[:5] and (tmscoreAct > float(0.2) and lenTarget >= self.minLenTemplate*float(1.2)) or (lenTarget >= self.minLenTemplate*float(1.15)):
                        auxAlign = align3[int(start):int(end)]
                        unitTemp = unit + '_' + str(temp)
                        hasBigGap, lostResiduesE, lostResiduesI = self.evaluate_align(auxAlign)
                        if hasBigGap and (lostResiduesE > 0 or lostResiduesI > 0):
                            works = []
                            works, resNum = self.create_new_temp_unit_gaps(unit, unitTemp, lostResiduesE, lostResiduesI)
                            if os.path.isfile(wdir + templatesWork + "units/" + unitTemp) and os.stat(wdir + templatesWork + "units/" + unitTemp).st_size > 5000:
                                unit = unitTemp
                                temp += 1
                                tmOutput = self.execute_tmalign(targetInOption, unitTemp, '.pdb')
                                values, tmscoreAct, rmsd, RmsdValue, start, end, gaps, gapsUnit, gapsUnitPerc, gapsPerc, lenUnit, auxAlign, lenTarget, auxAlignX, valuesLine, TMscoreLine2, align, align3, aliLenValue, coverage = self.assign_values_tmalign(tmOutput)
                                coverage = abs(float((aliLenValue - gapsUnit) + (aliLenValue - gaps))/float(int(end - start + 1) + int(lenUnit)))
                                if (lenUnit >= self.minLenTemplate*1.1 and lenTarget >= self.minLenTemplate*1.1) and (re.search("^-*[A-Z]{1,4}-{6,}[A-Z]{1,}.+",auxAlignX) or re.search(".+[A-Z]{1,}-{6,}[A-Z]{1,4}-*$",auxAlignX)):
                                    unitTemp2 = unit + '_' + str(temp)
                                    unit = unitTemp
                                    lostResiduesI, lostResiduesE = 0, 0
                                    hasBigGap, lostResiduesE, lostResiduesI = self.evaluate_align(auxAlignX)
                                    if hasBigGap and (lostResiduesE > 0 or lostResiduesI > 0):
                                        works,resNum = self.create_new_temp_unit_gaps(unit, unitTemp2, lostResiduesE, lostResiduesI)
                                        if os.path.isfile(wdir + templatesWork + "units/" + unitTemp2) and os.stat(wdir + templatesWork + "units/" + unitTemp2).st_size > 5000:
                                            unit = unitTemp2
                                            temp += 1
                                            tmOutput = self.execute_tmalign(targetInOption, unitTemp2, '.pdb ')
                                            values, tmscoreAct, rmsd, RmsdValue, start, end, gaps, gapsUnit, gapsUnitPerc, gapsPerc, lenUnit, auxAlign, lenTarget, auxAlignX, valuesLine, TMscoreLine2, align, align3, aliLenValue, coverage = self.assign_values_tmalign(tmOutput)
                                            coverage = abs(float((aliLenValue - gapsUnit) + (aliLenValue - gaps))/float(int(end - start + 1) + int(lenUnit)))
                        if tmscoreAct > self.LimTMscore3 and self.LimlenUnitTarget3 > lenTarget >= self.minLenTemplate*self.PorclenUnitTarget3 and lenUnit < self.LimlenUnitTarget3 and coverage >= self.coverageminfirstunit2:
                            templatesTMdataThirdRound.append([unit, round(tmscoreAct, 2), start, end, gapsPerc, gapsUnitPerc, contig, lenUnit, nicePosition, RmsdValue, len(align[:start]) - align[:start].count('-'), coverage])

                        if tmscoreAct > self.LimTMscore4 and self.LimlenUnitTarget4 > lenTarget >= self.minLenTemplate * self.PorclenTarget4 and self.LimlenUnitTarget4 > lenUnit >= self.minLenTemplate*self.PorclenUnit4 and coverage >= self.coverageminfirstunit3:
                            templatesTMdataLowerRound.append([unit, round(tmscoreAct, 2), start, end, gapsPerc, gapsUnitPerc, contig, lenUnit, nicePosition, RmsdValue, len(align[:start])-align[:start].count('-'), coverage])
        return templatesTMdataThirdRound, templatesTMdataLowerRound

    def create_templates_list(self, templatesTMdata):
        h = 0
        countTop = len(templatesTMdata) - 1
        templatesTMdata.sort(cmp=None)
        templatesTMdata.sort(key=itemgetter(7), reverse=True)
        countTopList = len(templatesTMdata) - 2
        for h in range(countTopList):
            if 0 >= h <= len(templatesTMdata)-1:
                valueTemp = templatesTMdata[h][2]
                m = countTop
                while m >= h + 1:
                    valueTempOther = templatesTMdata[m][2]
                    if float(valueTemp - 3) <= float(valueTempOther) <= float(valueTemp + 3):
                            templatesTMdata.remove(templatesTMdata[m])
                            countTop = len(templatesTMdata) - 1
                    m -= 1
            countTopList = len(templatesTMdata) - 2
        templatesTMdata.sort(cmp=None)
        if len(templatesTMdata) > 0:
            templatesTMdata.sort(key=itemgetter(1), reverse=True)
        return templatesTMdata

    def put_all_templates_tmdata(self, templatesTMdata, templatesTMdataSecondRound, templatesTMdataThirdRound, templatesTMdataLowerRound, particular):
        if not particular:
            templatesTMdata = templatesTMdata + templatesTMdataSecondRound + templatesTMdataThirdRound + templatesTMdataLowerRound
        lst = []
        new = []
        i = 0
        while i < len(templatesTMdata):
                if particular:
                    start = templatesTMdata[i][2]
                else:
                    start = templatesTMdata[i][10]
                if (start in lst) or (start + 1 in lst) or (start + 2 in lst) or (start + 3 in lst) or (start - 1 in lst) or (start - 2 in lst) or (start - 3 in lst)  :
                    templatesTMdata.remove(templatesTMdata[i])
                else:
                    lst.append(start)
                    new.append(templatesTMdata[i])
                    i += 1
        return new

    def calc_start_end_struct_alignment_particular(self, templatesTMdata2, i, target, rangeMINMAX,TMalignexe):
        wdir, templatesWork, directory = self.wdir, self.templatesWork, self.directory
        logging.debug("In calc start end particular")
        logging.debug(self.TMalignexe+' ' + wdir + templatesWork + target[:5] + '/' + target + '.pdb ' + wdir + templatesWork + "units/" + templatesTMdata2[0][0] + '  -o ' + wdir + templatesWork + target[:5] + '/TMscore' + target + templatesTMdata2[0][0] + ' -a T > ' + wdir + 'temp')

        if os.path.isfile(wdir + templatesWork + target[:5] + '/' + target + '.pdb') and os.path.isfile(wdir + templatesWork + "units/" + templatesTMdata2[0][0]):
            os.system(self.TMalignexe+' ' + wdir + templatesWork + target[:5] + '/' + target + '.pdb ' + wdir + templatesWork + "units/" + templatesTMdata2[0][0] + '  -o ' + wdir + templatesWork + target[:5] + '/TMscore' + target + templatesTMdata2[0][0] + ' -a T > ' + wdir + 'temp2')
        logging.debug("mensaje 1")

        if os.path.isfile(wdir+templatesWork + target[:5] + '/TMscore' + target + templatesTMdata2[0][0] + '_all_atm'):
            co = 0
            pdbFile = open(wdir+templatesWork + target[:5] + '/TMscore' + target + templatesTMdata2[0][0] +'_all_atm', 'r')
            parser = PDBParser()
            structureAux = parser.get_structure(target[:4], pdbFile)
            c = 0
            for ch in structureAux.get_chains():
                    if c == 0:
                            chx = ch.get_id()
                            c += 1
            for res in structureAux[0][chx].get_residues():
                if res.get_id()[0] == ' ':
                    for a in res:
                        if co == 0:
                            min_Possible = (res.get_id()[1])
                            co += 1
                        else:
                            max_Possible = (res.get_id()[1])

        if os.path.isfile(wdir+templatesWork + target[:5] + '/TMscore' + target + templatesTMdata2[0][0] + '_atm'):
            co = 0
            pdbFile = open(wdir+templatesWork + target[:5] + '/TMscore' + target + templatesTMdata2[0][0] +'_atm', 'r')
            parser = PDBParser()
            structureAux = parser.get_structure(target[:4], pdbFile)
            c = 0
            for ch in structureAux.get_chains():
                    if c == 0:
                            chx = ch.get_id()
                            c += 1
            for res in structureAux[0][chx].get_residues():
                if res.get_id()[0] == ' ':
                    for a in res:
                        if co == 0:
                            min_ = (res.get_id()[1])
                            co += 1
                        else:
                            max_ = (res.get_id()[1])
            if max_Possible-max_< 5:
                max_=max_Possible
            if min_-min_Possible< 5:
                min_ = min_Possible
            rangeMINMAX.append(min_)
            rangeMINMAX.append(max_)
            pdbFile.close()

    def create_pdb_fragments_particular(self, rangeMINMAX, j, structure, it, targetsListOpt, templateListOpt, fastpList, target, UsedTargetsList , InfoRanges, targetOriginal, worked):
        wdir, templatesWork, directory = self.wdir, self.templatesWork, self.directory

        spp = ' '
        spc = ''
        k = 1
        co = 0
        countLen = 0
        #print rangeMINMAX[0],rangeMINMAX[1]
        min_ = int(rangeMINMAX[0])
        max_ = int(rangeMINMAX[1])
        fastpFile = open(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.fastp', 'a')
        lenFile = open(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.len', 'a')
        if os.path.isfile(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.mus'):
            mustangFile = open(wdir+templatesWork +targetOriginal[:5] + '/' + targetOriginal[:5] + '.mus', 'a')
        else:
            mustangFile = open(wdir+templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.mus', 'w')
            mustangFile.write(">" + self.OUTdir + "\n")
        fpdbUnit = open(wdir + templatesWork + 'units/' + target[:5] + 'repeat' + str(it), 'w')
        initR, initL, last = -100, -100, -100
        right, left , coL, coR, countLeft, countRight = 0, 0, 0, 0, 0, 0
        fpdbUnitPre = open(wdir + templatesWork + target[:5] + "/" + target + 'pre.pdb', 'w')
        fpdbUnitPost = open(wdir + templatesWork + target[:5] + "/" + target + 'post.pdb', 'w')
        for res in structure[0][target[4:5]].get_residues():
            if co == 0:
                init = int(res.get_id()[1])
                co += 1
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
                    if int(res.get_id()[1]) in range(init, int(min_)):
                        left = 1
                        if coL == 0:
                            initL = int(res.get_id()[1])
                            coL += 1
                        if k >= 10000:
                                fpdbUnitPre.write('ATOM  ' + '%5s'%str(k) + '  ' + a.get_name() + spc+'%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id( )[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        else:
                                fpdbUnitPre.write('ATOM   ' + '%4s'%str(k) + '  '+ a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id( )[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        k += 1
                        if a.get_name() == 'C':
                            countLeft += 1
                    elif int(res.get_id()[1]) in range(int(min_), int(max_)+1):
                        if a.get_name() == 'C':
                            countLen += 1
                        if k >= 10000:
                                fpdbUnit.write('ATOM  ' + '%5s'%str(k) + '  '+a.get_name()+spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id()[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        else:
                                fpdbUnit.write('ATOM   ' + '%4s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id()[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        k += 1
                    else:
                        last = int(res.get_id()[1])
                        right = 1
                        if coR == 0:
                            initR = int(res.get_id()[1])
                            coR += 1
                        if k >= 10000:
                                fpdbUnitPost.write('ATOM  ' + '%5s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id( )[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        else:
                                fpdbUnitPost.write('ATOM   ' + '%4s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id( )[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        k += 1
                        if a.get_name() == 'C':
                            countRight += 1
        fpdbUnitPre.close()
        fpdbUnitPost.close()
        fpdbUnit.close()
        #print "after write pdb",countLeft, countRight

        if countLeft > 6:
            targetsListOpt.append(target + 'pre')
            if left == 1 and initL != -100:
                InfoRanges.append([target + 'pre.pdb', initL, int(min_) - 1])
        else:
            InfoRanges.append([target + 'pre.pdbXX', init, int(min_) - 1])

        if countRight > 6:
            targetsListOpt.append(target + 'post')
            InfoRanges.append([target + 'post.pdb', initR, last])
        else:
            if right == 1:
                InfoRanges.append([target + 'post.pdbXX', initR, last])
        fpdbUnit.close()
        #print "befor condition", initR,initL,right,left
        if (countLen > 6) and (((right == 1 or left == 1) and (initR != -100 or initL != -100)) or (right == 0 and left == 0)):
            concat = ''
            if left == 1 and init <= int(min_)-1 and (int(min_)-init) <= 6 and (self.is_after_of_repeat(InfoRanges, init - 1)):
                files = [wdir + templatesWork + target[:5] + '/' + target + 'pre.pdb', wdir + templatesWork+"units/" + target[:5] + 'repeat' + str(it)]
                concat = ''.join([open(ff).read() for ff in files])
                file = open(wdir + templatesWork + "units/" + target[:5] + 'repeat' + str(it), "w")
                file.write(concat)
                file.close()
                min_ = int(min_) - (int(min_) - init)
                os.remove(wdir + templatesWork + target[:5] + '/' + target + 'pre.pdb')
                if target+'pre' in targetsListOpt:
                    targetsListOpt.remove(target+'pre')
            if right == 1 and initR <= last and (last - initR) <= 6 and self.is_first_of_repeat(InfoRanges, last + 1):
                files = [wdir + templatesWork + "units/" + target[:5] + 'repeat' + str(it), wdir + templatesWork + target[:5] + '/' + target + 'post.pdb']
                concat = ''.join([open(ff).read() for ff in files])
                file = open(wdir + templatesWork + "units/" + target[:5] + 'repeat' + str(it), "w")
                file.write(concat)
                file.close()
                max_ += (last - initR + 1)
                os.remove(wdir + templatesWork + target[:5] + '/' + target + 'post.pdb')
                if target + 'post' in targetsListOpt:
                    targetsListOpt.remove(target + 'post')
            templateListOpt.append(target[:5] + 'repeat' + str(it))
            worked = True or worked
            mustangFile.write("+" + target[:5] + 'repeat' + str(it) + "\n")
            InfoRanges.append([target + 'repeat', int(min_), int(max_)])
            fastpFile.write(str(min_) + "\t" + str(max_) + "\n")
            lenFile.write(str(countLen) + "\n")
            UsedTargetsList.append(target)
            fastpList.append([min_, max_])
        else:
            worked = False or worked
            if os.path.isfile(wdir + templatesWork + "units/" + target[:5] + 'repeat' + str(it)):
                os.remove(wdir + templatesWork +"units/" + target[:5] + 'repeat' + str(it))
            if os.path.isfile(wdir + templatesWork + target[:5] + '/' + target[:5] + 'repeat' + str(it)):
                os.remove(wdir + templatesWork + target[:5] + '/' + target[:5] + 'repeat' + str(it))
        fastpFile.close()
        lenFile.close()
        mustangFile.close()
        return worked,templateListOpt

    def verify_holes(self, target, fastpListAA, sc,strins):
        wdir, directory, templatesWork, PDBdir = self.wdir, self.directory, self.templatesWork, self.PDBdir
        fastpL, holesL = [], []
        fastOrig = fastpListAA
        sum, valmin = 0, 999
        for rr in list(fastpListAA):
            if rr != '':
                row = rr.split('\t')
                max_ = int(row[1])
                min_ = int(row[0])
                sum += (max_ - min_)
                if (max_-min_) < valmin:
                    valmin = (max_ - min_)
                fastpL.append([min_, max_])
            else:
                fastpListAA.remove('')
        avgLen = sum/len(fastpL)
        fastpL = sorted(fastpL)
        nholes = 1
        j = 1
        for row in range(len(fastpL) - 1):

            if int(fastpL[j][0])-int(fastpL[j-1][1]) > 1:
                nholes += 1
            j += 1

        if sc != "V" and nholes > len(fastpL) * 0.7:
            return fastOrig, strins
        j, i, itt = 1, 1, 1

        for row in range(len(fastpL) - 1):
            if nholes <= len(fastpL) * 0.7:
                worked = False
                if int(fastpL[j][0])-int(fastpL[j-1][1]) > 1:
                    if (int(fastpL[j][0]) - int(fastpL[j-1][1])-1) < valmin * 0.85:
                        tmstrr = str(int(fastpL[j-1][1])+1) + " " + str(int(fastpL[j][0])-1)
                        fastpL, fastpListAA, worked = self.find_insertion_position(fastpL, j, fastpListAA, int(fastpL[j][0]) - int(fastpL[j-1][1]), itt, target, wdir, templatesWork, valmin)
                        if worked:
                            strins = strins + tmstrr + ";"
                        itt += 1
                    else:
                        holesL.append([int(fastpL[j-1][1]) + 1, int(fastpL[j][0]) - 1])
                    i += 1
            j += 1
        new, it = 1, 1
        tempName = 'AUXrepeat'
        for row in holesL:
            itsUnit = False
            max_ = int(row[1])
            min_ = int(row[0])
            countAux = self.create_aux_unit(min_, max_, wdir + templatesWork, wdir + templatesWork + "units/", target, it, tempName)
            for count in range(len(fastpL)):
                if not itsUnit and countAux > 6 and os.path.isfile(wdir + templatesWork + 'units/' + target[:5] + 'AUXrepeat' + str(it)) and os.path.isfile(wdir + templatesWork + "units/" + target[:5] + 'repeat' + str(count+1)):
                    tmOutput = self.execute_tmalign( target[:5] + 'repeat' + str(count+1), target[:5] + 'AUXrepeat' + str(it), '')
                    dataList = tmOutput.split('\n')
                    TMscoreLine2 = float(dataList[18].split(' ')[1])
                    TMscoreLine1 = float(dataList[17].split(' ')[1])
                    alilen= float((dataList[16].split(',')[0]).split("=")[1])
                    ali = dataList[23]
                    ali2 = dataList[25]
                    countend = 0
                    countstart = 0
                    if TMscoreLine1 >= float(0.5) and float(ali.count('-')) >= (0.5 * valmin) and (re.search("-+$", ali) or re.search("^-+", ali)):
                        if re.search("-+$", ali):
                            countend = len(ali)-re.search("-+$", ali).start()
                        else:
                            countstart = re.search("^-+", ali).end()
                        if max_ - min_ - alilen >= valmin * 0.85 and (countend > valmin * 0.85 or countstart > valmin * 0.85):
                            itsUnit = True
                            if countstart > valmin * 0.85:
                                holesL.append([min_, min_ + re.search("^-+", ali).end()])
                                holesL.append([int(min_ + re.search("^-+", ali).end() + 1), max_])
                            elif countend > valmin * 0.85:
                                holesL.append([min_, min_+re.search("-+$", ali).start() - 1])
                                holesL.append([re.search("-+$", ali).start()+min_, max_])
                    else:
                        if TMscoreLine2 >= float(0.35) and (alilen >= valmin*0.8 or alilen >= countAux * 0.8):
                            itsUnit = True
                            os.system("mv " + wdir + templatesWork + "units/" + target[:5] + 'AUXrepeat' + str(it) +' ' + wdir + templatesWork + "units/" + target[:5] + 'repeat' + str(len(fastpL) + new))
                            if os.path.isfile(wdir + templatesWork + target[:5] + '/' + target[:5] + '.mus'):
                                mustangFile = open(wdir + templatesWork + target[:5] + '/' + target[:5] + '.mus', 'a')
                            else:
                                mustangFile = open(directory + templatesWork + target[:5] + '/' + target[:5] + '.mus', 'w')
                                mustangFile.write(">" + self.OUTdir + "\n")
                            mustangFile.write("+" + target[:5] + 'repeat' + str(len(fastpL) + new) + "\n")
                            new += 1
                            fastpListAA.append(str(min_) + '\t' + str(max_))
                            fastpL.append([min_, max_])
                            mustangFile.close()
                            if re.search("(-){3}", ali) and re.search("(-){3}", ali).start() > 0.2 * countAux:
                                strins = self.obtain_str_insertions(ali, ali2, [min_, max_], strins)
            if os.path.isfile(wdir + templatesWork + "units/" + target[:5] + 'AUXrepeat' + str(it)):
                os.system("rm " + wdir + templatesWork + "units/" + target[:5] + 'AUXrepeat' + str(it))
            it += 1
        if len(strins) == 1:
            return fastpListAA, strins
        else:
            return fastpListAA, strins

    def find_insertion_position(self, fastpL, i, fastpListAA, value, it, target, directory, templatesWork, valmin):
        tempName = "AUXINS"
        itsUnit = False
        MaxTMscoreLine2first, MaxTMscoreLine2second = 0, 0
        itsUnit = False
        if not itsUnit:
            resCount2 = self.create_aux_unit(fastpL[i-1][0], (fastpL[i-1][1]) + int(value) - 1, directory + templatesWork, directory + templatesWork + "units/", target, it, tempName)
            for count in range(len(fastpL)):
                if count != fastpListAA.index(str(fastpL[i-1][0])+"\t"+str(fastpL[i-1][1])) and not(itsUnit) and os.path.isfile(directory + templatesWork + 'units/' + target[:5] + tempName + str(it)) and os.path.isfile(directory + templatesWork + "units/" + target[:5] + 'repeat' + str(count + 1)):
                        tmOutput = self.execute_tmalign(target[:5] + 'repeat' + str(count+1), target[:5] + tempName + str(it), '')
                        dataList = tmOutput.split('\n')
                        TMscoreLine2 = float(dataList[18].split(' ')[1])
                        aliLen = float(((dataList[16].split(',')[0]).split('='))[1])
                        if TMscoreLine2 > MaxTMscoreLine2second and TMscoreLine2 >= self.MinTMscoreInsertion and (resCount2 * self.procLenInsertion <= aliLen or valmin <= aliLen):
                            MaxTMscoreLine2second = TMscoreLine2
                            Maxpossec = i
        itsUnit = False
        if not itsUnit:
            resCount = self.create_aux_unit((fastpL[i][0])-int(value) + 1, fastpL[i][1], directory + templatesWork, directory + templatesWork + "units/", target, it, tempName)
            for count in range(len(fastpL)):
                if count != fastpListAA.index(str(fastpL[i][0]) + "\t" + str(fastpL[i][1])) and not itsUnit and os.path.isfile(directory + templatesWork + 'units/' + target[:5] + tempName + str(it)) and os.path.isfile(directory + templatesWork + "units/" + target[:5] + 'repeat' + str(count + 1)):
                        tmOutput = self.execute_tmalign(target[:5] + 'repeat' + str(count + 1), target[:5] + tempName + str(it), '')
                        dataList = tmOutput.split('\n')
                        aliLen =float(((dataList[16].split(',')[0]).split('='))[1])
                        TMscoreLine2 = float(dataList[18].split(' ')[1])
                        if TMscoreLine2 > MaxTMscoreLine2first and TMscoreLine2 >= self.MinTMscoreInsertion and (resCount * self.procLenInsertion <= aliLen or valmin <= aliLen):
                            MaxTMscoreLine2first = TMscoreLine2
                            Maxposfirst = i
        if MaxTMscoreLine2first > MaxTMscoreLine2second:
            i = Maxposfirst
            fastpListAA[fastpListAA.index(str(fastpL[i][0]) + "\t" + str(fastpL[i][1]))] = str((fastpL[i][0]) - int(value) + 1) + '\t' + str(fastpL[i][1])
            fastpL[i] = ([fastpL[i][0] - int(value) + 1, fastpL[i][1]])
            self.create_aux_unit(fastpL[i][0], fastpL[i][1], directory+templatesWork, directory + templatesWork + "units/", target, 1 + fastpListAA.index(str(fastpL[i][0]) + "\t" + str(fastpL[i][1])), 'repeat')
        elif MaxTMscoreLine2second != 0:
            i = Maxpossec
            fastpListAA[fastpListAA.index(str(fastpL[i-1][0]) + "\t" + str(fastpL[i - 1][1]))] = str(fastpL[i - 1][0]) + '\t' + str((fastpL[i - 1][1]) + int(value) - 1)
            fastpL[i - 1] = ([fastpL[i - 1][0], fastpL[i - 1][1] + int(value) - 1])
            self.create_aux_unit(fastpL[i - 1][0],fastpL[i - 1][1], directory + templatesWork, directory + templatesWork + "units/", target, 1 + fastpListAA.index(str(fastpL[i - 1][0]) + "\t" + str(fastpL[i - 1][1])), 'repeat')
        return fastpL, fastpListAA, itsUnit

    def calc_holes(self, fastpListAA):
        counterHoles, end = 0, 0
        fastpListAux = []
        for pair in fastpListAA:
            if len(pair) > 2:
                dat = pair.split('\t')
                fastpListAux.append([int(dat[0]), int(dat[1])])
        fastpListAux.sort(key=itemgetter(0))
        count = 0
        for pair in fastpListAux:
            if count == 0:
                end = int(pair[1])
                count += 1
            else:
                counterHoles += int(pair[0])-end
                end = int(pair[1])
                count += 1
        return counterHoles

    def contiguous_set(self, fastpListAA):
        fastpListAux = []
        contFastp, i = 0, 0
        for pair in fastpListAA:
            if len(pair) > 0:
                dat = pair.split('\t')
                fastpListAux.append([int(dat[0]), int(dat[1])])
        fastpListAux.sort(key=itemgetter(0))
        for pair in fastpListAux:
            if i == 0:
               start = int(pair[0])
               end = int(pair[1])
            elif len(fastpListAux) <= 1:
                    return 0
            elif (start - int(pair[1]) <= 1 and (start - int(pair[1]) > 0)) or ((int(pair[0]) - end <= 1) and (int(pair[0]) - end > 0)):
                    contFastp += 1
            start = int(pair[0])
            end = int(pair[1])
            i += 1
        if i >= 2:
            return float(float(contFastp)/float(i - 1))
        else:
            return 0

    def contiguous_set_all(self, target, TEMPLATESround, destinationPath, info, templatesTMdata, clas, strins):
        wdir, directory = self.wdir, self.directory
        min_ = 9999
        done = False
        max_, maxLen = 0, 0
        maxlensel=[]
        fastpListAASel = []
        maxstrins = ""
        maxindx = 0
        ind=-1
        for i in info:
            validate = i[2]
            actholes = i[1]
            strins = i[5]
            inSecondRound = i[3]
            fastpListAA = list(i[4])
            if validate > max_ and actholes < min_ and (len(fastpListAA) >= 4 and ((validate >= float(0.6) and actholes <= 18) or (len(fastpListAA) >= 4 and validate >= float(0.6) and actholes <= 40) or (validate >= float(0.75) and len(fastpListAA) >= 4))):
                max_ = validate
                min_ = actholes
                ind = int(i[0])
                fastpListAASel = list(i[4])
            if len(fastpListAA) > maxLen:
                maxlensel = i[4]
                maxLen = len(fastpListAA)
                maxindx = int(i[0])
                maxstrins = i[5]
        strUnits = ' '
        fastpL = []
        j = 0
        if ind >= 0:
            if fastpListAASel == [] :
                done = True
                strins = maxstrins
                ind = int(maxindx)
                fastpListAASel2 = maxlensel
                for rr in (fastpListAASel2):
                    if rr != '':
                        row = rr.split('\t')
                        max_ = int(row[1])
                        min_ = int(row[0])
                        fastpL.append([min_, max_])
                if fastpL[0][0]:
                    region = str(fastpL[0][0]) + '-'
                fastpListAA = []
                for r in fastpListAASel2:
                    rr = r.split("\t")
                    fastpListAA.append(str(rr[0]) + "\t" + str(rr[1]))
                    if j == len(fastpL) - 1:
                        region = region + str(fastpL[j][1])
                    else:
                        if j < len(fastpL) - 1 and fastpL[j + 1][0] - fastpL[j][1] > 1:
                            region = region + str(fastpL[j][1]) + ';' + str(fastpL[j + 1][0]) + " "
                    j += 1
                    fastpListAASel = list(fastpListAA)
            if fastpListAASel != [] and done == False:

                fastpListAASel, region ,insertions,ninsertions = sort_fastp_result(fastpListAASel)
            if fastpListAASel != [] :
                for row in fastpListAASel:
                    if len(row) > 0:
                        strUnits = strUnits + ";" + row.split("\t")[0]+" " + row.split("\t")[1]
                strUnits = strUnits[1:]
                result, sc = "", ""
                if clas == 'III':
                    xxx = True
                    result, sc = identify_classification_class_iii(target, directory,   templatesTMdata, ind, strUnits, region, xxx, strins)
                if clas == 'IV':
                    result,sc = identify_classification_class_iv(target, directory,   templatesTMdata, ind, strUnits, region, True, strins)
                if clas == 'V':
                    result, sc = identify_classification_class_v(target, directory,  templatesTMdata, ind, strUnits, region, True, strins)

                origenUnitPath = TEMPLATESround[ind] + "/"
                origenPath = TEMPLATESround[ind] + "/"
                move_files(origenUnitPath, target, destinationPath, origenPath)
                if os.path.isfile(destinationPath + target[:5] + '*_*'):
                    os.remove(destinationPath + target[:5] + '*_*')
                return result
        return ""

    def evaluate_align(self, auxAlign):
        lostResiduesI, lostResiduesE, posFirstRes, lostResiduesE, posFirstRes = 0, 0, 0, 0, 0
        IsGap = False
        if re.search("^-*[A-Z]{1,6}-{6,}[A-Z]{1,}.+", auxAlign) or re.search(".+[A-Z]{1,}-{6,}[A-Z]{1,6}-*$", auxAlign) or re.search("-{20,}[A-Z]{1,7}$", auxAlign):
            IsGap = True
        if IsGap:
            if re.search("^-*[A-Z]{1,6}-{6,}[A-Z]{1,}.+", auxAlign):
                posFirstRes = re.search("^-*[A-Z]{1,6}-{6,}[A-Z]{1,}.+", auxAlign).start()
                posFirstGap = re.search("-", auxAlign).start()
                lostResiduesI = posFirstGap-posFirstRes
            if re.search(".+[A-Z]{1,}-{6,}[A-Z]{1,6}-*$", auxAlign):
                posFirstLastRes = re.search("-{6,}[A-Z]{1,6}-*$", auxAlign).start()
                lostResiduesE = (len(auxAlign[posFirstLastRes:]))-(auxAlign[posFirstLastRes:].count('-'))
            if re.search("-{20,}[A-Z]{1,7}$", auxAlign):
                posFirstLastRes = re.search("-{6,}[A-Z]{1,7}$", auxAlign).start()
                lostResiduesE = (len(auxAlign[posFirstLastRes:])) - (auxAlign[posFirstLastRes:].count('-'))
        return IsGap, lostResiduesE, lostResiduesI

    def execute_tmalign(self, targetInOption,  unitTemp, ext ):
        wdir, templatesWork, directory, PDBdir = self.wdir, self.templatesWork, self.directory, self.PDBdir
        if ext == '.pdb':
            if os.path.isfile(wdir + templatesWork + targetInOption[:5] + '/' + targetInOption + '.pdb'):
                tmOutput = subprocess.check_output(self.TMalignexe + ' '  + wdir + templatesWork + targetInOption[:5] + '/' + targetInOption + ext + ' '+wdir + templatesWork + "units/" + unitTemp+' -a T  ', shell=True)
            else:
                tmOutput = subprocess.check_output(self.TMalignexe + ' '  + PDBdir + targetInOption + ext + ' ' +wdir + templatesWork + "units/" + unitTemp + ' -a  T ', shell=True)
        else:
            tmOutput = subprocess.check_output(self.TMalignexe + ' '  + wdir + templatesWork + "units/" + targetInOption+' ' + wdir + templatesWork + 'units/' + unitTemp + ' -a  T ', shell=True)
        return tmOutput

    def create_aux_unit(self, min_, max_, origen, destination, target, it, tempName):

        fpdbUnit = open(destination + tempName + str(it), 'w')
        if os.path.isfile(origen + target[:5] + target[:5] + ".pdb"):
            pdbFile = open(origen + target[:5] + target[:5] + ".pdb", 'r')
        elif os.path.isfile(origen + target[:5] +".pdb"):
            pdbFile = open(origen + target[:5] + ".pdb", 'r')
        elif os.path.isfile(origen + target):
            pdbFile = open(origen + target, 'r')
        else:
            pdbFile = open(origen + target + '/' + target[:5] + ".pdb", 'r')
        parser = PDBParser()
        k = 1
        countAux = 0
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
                    if int(res.get_id()[1]) in range(int(min_), int(max_)+1):
                        if k >= 10000:
                                fpdbUnit.write('ATOM  ' + '%5s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + '  ' + target[4:5] + '%4s'%str(res.get_id()[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        else:
                                fpdbUnit.write('ATOM   ' + '%4s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id()[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        k += 1
                        if a.get_name() == 'C':
                            countAux += 1
        fpdbUnit.close()
        return countAux

    def assign_values_tmalign(self, tmOutput):
        dataList = tmOutput.split('\n')
        lenline, lenline2, valuesLine, TMscoreLine1, TMscoreLine2, TMscoreLine3, align, align2, align3 = dataList[13], dataList[14], dataList[16], dataList[17], dataList[18], dataList[19], dataList[23], dataList[24], dataList[25]
        values = valuesLine.split(",")
        tmscoreAct = float(TMscoreLine2.split(' ')[1])
        rmsd = values[1].split("=")
        RmsdValue = float(rmsd[1])
        aliLen = values[0].split("=")
        aliLenValue = float(aliLen[1])
        start = re.search('[A-Z]', align3).start()
        if start < 15 :
            start = 0
        end = re.search("[A-Z]-*$", align3).start()+1
        gaps = align[int(start):int(end)].count('-')
        gapsUnit = align3[int(start):int(end)].count('-')
        gapsUnitPerc = float((float(gapsUnit)/float(int(end)-int(start))))
        gapsPerc = float((float(gaps)/float(int(end)-int(start))))
        lenUnit = int(end)-int(start)-int(gapsUnit)-int(gaps)
        auxAlign = align3[int(start):int(end)]
        lenTarget = len(align[int(start):int(end)])-gaps
        auxAlignX = align3[int(start):int(end)]
        lenlin2 = lenline2.split(":")[1]
        lenline2 = lenlin2.split("r")[0]
        lenlin3 = lenline.split(":")[1]
        lenline3 = lenlin3.split("r")[0]
        coverage = abs(float(2*int(aliLenValue))/float((int(lenline2)+int(lenTarget))))

        return values, tmscoreAct, rmsd, RmsdValue, start, end, gaps, gapsUnit, gapsUnitPerc, gapsPerc, lenUnit, auxAlign, lenTarget, auxAlignX, valuesLine, TMscoreLine2, align, align3, aliLenValue, coverage

    def verify_second_structure(self, directory, templatesWork, targetInOption, start, end, gaps,Dsspexe):
        os.system(Dsspexe+' -i ' + directory + templatesWork + targetInOption[:5] + "/" + targetInOption + '.pdb -o ' + directory + templatesWork + targetInOption + '.dssp')
        if os.path.isfile(directory + templatesWork + targetInOption + '.dssp'):
            with open(directory + templatesWork + targetInOption + '.dssp', 'r') as f:
                dsspinfo = f.read().split("\n")
            ignored = False
            dssp = ""
            for row in dsspinfo:
                if ignored:
                    dssp += row[16:17]
                elif row[:3] == "  #":
                    ignored = True
        else:
            dssp = ""
        inispace, endspace = 0, 0

        if ( len(dssp)> 0) and ((start > 0 and ((dssp[start - 1] == "H" and dssp[start] == "H" or
                                        dssp[start - 1] == "E" and dssp[start] == "E" or
                                        dssp[start - 1] == "T" and dssp[start] == "H" or
                                        dssp[start - 1] == "H" and dssp[start] == "T"))) or (
                    ((len(dssp) > end - gaps - 1) and (len(dssp) > end - gaps + 1)) or
            (end < len(dssp) and end +1 - gaps < len(dssp) and (dssp[end + 1 - gaps] == "H" and dssp[end - gaps] == "H" or
                                        dssp[end + 1 - gaps] == "E" and dssp[end - gaps] == "E" or
                                    dssp[end + 1 - gaps] == "T" and dssp[end - gaps] == "H" or
                                dssp[end + 1 - gaps] == "H" and dssp[end - gaps] == "T")))):
            if start - 1 > 0:
                inispace = len(dssp[:start - 1]) - re.search(" ", dssp[:start - 1]).end()
            if end + 1 < (len(dssp)):
                endspace = re.search(" ", dssp[end + 1:]).start()
            if inispace > 2:
                inispace = 0
            if endspace > 2:
                endspace = 0
            return 1, inispace, endspace
        else:
            return 0, inispace, endspace

    def add_data_to_templates_tmdata_particular(self, tmOutput, templatesTMdata2, fastpList, unit, temp, works, targetInOption, i, ConsitionSetPart, lenFirstUnit,sc,Dsspexe):
        wdir, PDBdir, directory, templatesWork = self.wdir, self.PDBdir, self.directory, self.templatesWork
        flagInsert = 0
        values, tmscoreAct, rmsd, RmsdValue, start, end, gaps, gapsUnit, gapsUnitPerc, gapsPerc, lenUnit, auxAlign, lenTarget, auxAlignX, valuesLine, TMscoreLine2, align, align3, aliLenValue, coverage = self.assign_values_tmalign(tmOutput)
        if sc =='III':
            maxRMSDbyClass = 2.0
            minTMscore = 0.25
        else:
            maxRMSDbyClass = 4.0
            minTMscore = 0.28
        if float(TMscoreLine2.split()[1]) > float(self.minTMscore) and aliLenValue >= (
            float(self.minLenTemplate) * float(self.minalilen)) and RmsdValue < maxRMSDbyClass and float(TMscoreLine2.split()[1]) > minTMscore:

            auxAlign = align3[int(start):int(end)]
            contig = self.is_contiguos(start, end, fastpList)
            lenUnit = int(end) - int(start) - int(gapsUnit) - int(gaps)
            diffEnd = len(align) - end
            if start < diffEnd and start < 4:
                nicePosition = 1
            elif start > diffEnd and diffEnd < 4:
                nicePosition = 1
            else:
                nicePosition = 0
            auxAlign = align3[int(start):int(end)]

            auxAlignX = auxAlign
            auxAligntar = align[int(start):int(end)]
            unitTemp = unit + '_' + str(temp)
            lostResiduesI ,lostResiduesE = 0,0

            if ConsitionSetPart != 4 and (re.search("[A-Z]-*[A-Z]{1,5}$", auxAlign)or re.search("^-*[A-Z]{1,5}-*[A-Z]{1,5}", auxAlign)):
                if re.search("[A-Z]-*[A-Z]{1,5}$", auxAlign):
                    inicio=re.search("[A-Z]-*[A-Z]+$", auxAlign).start()
                    fin = re.search("[A-Z]+$", auxAlign).start()
                    lastgaps = auxAlign[inicio:].count('-')
                    if lastgaps > 3:

                        lostResiduesE=len(auxAlign[fin:])

                if start!=0 and re.search("^-*[A-Z]{1,5}-*[A-Z]{1,5}", auxAlign):
                    inicio = re.search("[A-Z]{1,5}-*[A-Z]+", auxAlign).start()
                    fin =  re.search("[A-Z]{1,5}-*[A-Z]", auxAlign).end()
                    lastgapsS = auxAlign[inicio:fin].count('-')
                    if lastgapsS > 3 :

                        lostResiduesI = len(auxAlign[inicio:fin])-lastgapsS-1

                if     lostResiduesI !=0 or lostResiduesE!=0:
                        lostResiduesI=0
                        works, resNum = self.create_new_temp_unit_gaps_part(unit, unitTemp, lostResiduesE, lostResiduesI,wdir, templatesWork)
                        unit = unitTemp
                        temp += 1
                        if works[0] == 1 and os.path.isfile(wdir + templatesWork + 'units/' + unitTemp) and resNum > 5:
                            works = []
                            tmOutput = self.execute_tmalign(targetInOption, unitTemp, '.pdb')

                            values, tmscoreAct, rmsd, RmsdValue, start, end, gaps, gapsUnit, gapsUnitPerc, gapsPerc, lenUnit, auxAlign, lenTarget, auxAlignX, valuesLine, TMscoreLine2, align, align3, aliLenValue, coverage = self.assign_values_tmalign(tmOutput)
                            auxAlignX = align3[int(start):int(end)]
                            auxAligntar = align[int(start):int(end)]


            if re.search("-{10,}[A-Z]{1,7}$", auxAlign) and ConsitionSetPart == 4:
                lostResiduesI, lostResiduesE = 0, 0
                hasBigGap, lostResiduesE, lostResiduesI = self.evaluate_align(auxAlign)
                if hasBigGap and (lostResiduesE > 0 or lostResiduesI > 0):
                    works, resNum = self.create_new_temp_unit_gaps_part(unit, unitTemp, lostResiduesE, lostResiduesI, wdir, templatesWork)
                    unit = unitTemp
                    temp += 1
                    if works[0] == 1 and os.path.isfile(wdir + templatesWork + 'units/' + unitTemp) and resNum > 5:
                        works = []
                        tmOutput = self.execute_tmalign(targetInOption, unitTemp, '.pdb')

                        values, tmscoreAct, rmsd, RmsdValue, start, end, gaps, gapsUnit, gapsUnitPerc, gapsPerc, lenUnit, auxAlign, lenTarget, auxAlignX, valuesLine, TMscoreLine2, align, align3, aliLenValue, coverage = self.assign_values_tmalign(tmOutput)
                        auxAlignX = align3[int(start):int(end)]
                        auxAligntar = align[int(start):int(end)]
                        if re.search("^-*[A-Z]{1,4}-{6,}[A-Z]{1,}.+", auxAlignX) or re.search(
                                ".+[A-Z]{1,}-{6,}[A-Z]{1,4}-*$", auxAlignX):
                            unitTemp2 = unit + '_' + str(temp)
                            unit = unitTemp
                            lostResiduesI, lostResiduesE = 0, 0
                            hasBigGap, lostResiduesE, lostResiduesI = self.evaluate_align(auxAlignX)
                            if hasBigGap and (lostResiduesE > 0 or lostResiduesI > 0):
                                works, resNum = self.create_new_temp_unit_gaps_part(unit, unitTemp2, lostResiduesE, lostResiduesI, directory, templatesWork)
                                unit = unitTemp2
                                temp += 1
                                if works[0] == 1 and os.path.isfile(wdir + templatesWork + 'units/' + unitTemp2) and os.stat(wdir + templatesWork + 'units/' + unitTemp2).st_size > 5000:
                                    works = []
                                    tmOutput = self.execute_tmalign(targetInOption, unitTemp2, '.pdb')
                                    values, tmscoreAct, rmsd, RmsdValue, start, end, gaps, gapsUnit, gapsUnitPerc, gapsPerc, lenUnit, auxAlign, lenTarget, auxAlignX, valuesLine, TMscoreLine2, align, align3, aliLenValue, coverage = self.assign_values_tmalign(
                                        tmOutput)
                                    auxAlignX = align3[int(start):int(end)]
                                    auxAligntar = align[int(start):int(end)]
            unitTemp3 = unit + '_' + str(temp)
            isokSecStruc, inival, endval = self.verify_second_structure(wdir, templatesWork, targetInOption, start, end, gaps,Dsspexe)
            if isokSecStruc == 1:
                works, resNum = self.create_new_temp_unit_gaps_part(unit, unitTemp3, endval, -inival, wdir, templatesWork)
                unit = unitTemp3
                temp += 1
                if works[0] == 1 and os.path.isfile(wdir + templatesWork + 'units/' + unitTemp3) and os.stat(wdir + templatesWork + 'units/' + unitTemp3).st_size > 5000:
                    works = []
                    tmOutput = self.execute_tmalign(targetInOption, unit, '.pdb')

                    values, tmscoreAct, rmsd, RmsdValue, start, end, gaps, gapsUnit, gapsUnitPerc, gapsPerc, lenUnit, auxAlign, lenTarget, auxAlignX, valuesLine, TMscoreLine2, align, align3, aliLenValue, coverage = self.assign_values_tmalign(tmOutput)
                    auxAlignX = align3[int(start):int(end)]
                    auxAligntar = align[int(start):int(end)]
            if (works == [] or works[0] == 1) and lenTarget >= float(self.minLenTemplate) * float(
                    self.porclenminPart) and coverage >= float(self.minCoverage) and gapsUnitPerc < float(0.5) and lenUnit>15:
                isokSecStruc = 0
                if isokSecStruc == 0 and ConsitionSetPart == 1:
                    if lenUnit >= float(self.minLenTemplate) * float(self.porcLenPart1) and float(
                            TMscoreLine2.split()[1]) > float(self.tmscorePart1) and lenUnit >= float(
                            lenFirstUnit) * float(self.porcLenUnitPart1) and lenTarget >= float(
                            lenFirstUnit) * float(self.porcLenTargatPart1) and coverage > float(
                            self.coveragePart1) and int(end - start) <= int(self.maxlen1) and float(RmsdValue) < float(3.0) :

                        templatesTMdata2.append(
                            [unit, float(TMscoreLine2.split()[1]), start, end, gapsPerc, gapsUnitPerc, contig,
                             lenUnit, nicePosition, flagInsert, RmsdValue, targetInOption,
                             len(align[:start]) - align[:start].count('-'), coverage, auxAlignX, auxAligntar])
                elif isokSecStruc == 0 and ConsitionSetPart == 2:
                    if lenUnit >= float(self.minLenTemplate) * float(self.porcLenPart2) and float(
                            TMscoreLine2.split()[1]) > float(self.tmscorePart2) and lenUnit >= lenFirstUnit * float(
                            self.porcLenUnitPart2) and lenTarget >= float(lenFirstUnit) * float(
                            self.porcLenTargatPart2) and coverage > float(self.coveragePart2) and int(
                                    end - start) <= int(self.maxlen2) and float(RmsdValue) < float(3.0) :

                        templatesTMdata2.append(
                            [unit, float(TMscoreLine2.split()[1]), start, end, gapsPerc, gapsUnitPerc, contig,
                             lenUnit, nicePosition, flagInsert, RmsdValue, targetInOption,
                             len(align[:start]) - align[:start].count('-'), coverage, auxAlignX, auxAligntar])
                elif ConsitionSetPart == 3:
                    if lenUnit >= float(self.minLenTemplate) * float(self.porcLenPart3) and float(
                            TMscoreLine2.split()[1]) > float(self.tmscorePart3) and float(coverage) > float(
                            self.coveragePart3) and int(end - start) <= int(self.maxlen3) and float(RmsdValue) < float(3.0) :

                        templatesTMdata2.append(
                            [unit, float(TMscoreLine2.split()[1]), start, end, gapsPerc, gapsUnitPerc, contig,
                             lenUnit, nicePosition, flagInsert, RmsdValue, targetInOption,
                             len(align[:start]) - align[:start].count('-'), coverage, auxAlignX, auxAligntar])
                elif ConsitionSetPart == 4:
                    if float(lenTarget) >= float(self.minLenTemplate) * float(self.porcLenPart4) and float(
                            TMscoreLine2.split()[1]) > float(self.tmscorePart4) and float(coverage) > float(
                            self.coveragePart4) and int(end - start) <= int(self.maxlen4) and float(RmsdValue) < float(3.0) :

                        templatesTMdata2.append(
                            [unit, float(TMscoreLine2.split()[1]), start, end, gapsPerc, gapsUnitPerc, contig,
                             lenUnit, nicePosition, flagInsert, RmsdValue, targetInOption,
                             len(align[:start]) - align[:start].count('-'), coverage, auxAlignX, auxAligntar])
            #print "valores", unit, str(TMscoreLine2.split()[1]), start, end, gapsPerc, gapsUnitPerc, contig,lenUnit, nicePosition, flagInsert, RmsdValue, targetInOption, len(align[:start]) - align[:start].count('-'), coverage, auxAlignX, auxAligntar
        return templatesTMdata2

    def calc_start_end_struct_alignment(self, templatesTMdata, i, rangeMINMAX):
        wdir, SRULdir, directory, templatesWork, target = self.wdir, self.SRULdir, self.directory, self.templatesWork, self.target
        # print 'TMalign '+directory+templatesWork+target[:5]+'/'+target+'.pdb '+SRULdir+templatesTMdata[i][0]
        logging.debug("In calc start end particular")
        if os.path.isfile(SRULdir + templatesTMdata[i][0]) and os.path.isfile(wdir + templatesWork + target[:5] + '/' + target + '.pdb'):
            os.system(self.TMalignexe + ' ' + wdir + templatesWork + target[:5] + '/' + target + '.pdb ' + SRULdir +templatesTMdata[i][0] + '  -o ' + wdir + templatesWork + target[:5] + '/TMscore' + target + templatesTMdata[i][0] + ' -a T > ' + wdir + '/temp')
            logging.debug(self.TMalignexe + ' ' + wdir + templatesWork + target[:5] + '/' + target + '.pdb ' + SRULdir +
                      templatesTMdata[i][0] + '  -o ' + wdir + templatesWork + target[:5] + '/TMscore' + target +
                      templatesTMdata[i][0] + ' -a T > ' + wdir + '/temp')

        else :
            if os.path.isfile(wdir + templatesWork + target[:5] + '/' + target + '.pdb') and os.path.isfile(wdir + templatesWork + 'units/' + templatesTMdata[i][0]):
                os.system(self.TMalignexe + ' '  + wdir + templatesWork + target[
                                                              :5] + '/' + target + '.pdb ' + wdir + templatesWork + 'units/' + templatesTMdata[i][0] + '  -o ' + wdir + templatesWork + target[:5] + '/TMscore' + target + templatesTMdata[i][0] + ' -a T > ' + wdir + '/temp')
                logging.debug(self.TMalignexe + ' '  + wdir + templatesWork + target[
                                                              :5] + '/' + target + '.pdb ' + wdir + templatesWork + 'units/' + templatesTMdata[i][0] + '  -o ' + wdir + templatesWork + target[:5] + '/TMscore' + target + templatesTMdata[i][0] + ' -a T > ' + wdir + '/temp')

        if os.path.isfile(wdir + templatesWork + target[:5] + '/TMscore' + target + templatesTMdata[i][0] + '_atm'):
            co, c = 0, 0

            pdbFile = open(wdir + templatesWork + target[:5] + '/TMscore' + target + templatesTMdata[i][0] + '_atm', 'r')
            parser = PDBParser()
            structureAux = parser.get_structure(target[:4], pdbFile)
            for ch in structureAux.get_chains():
                if c == 0:
                    chx = ch.get_id()
                    c += 1
            for res in structureAux[0][chx].get_residues():
                if res.get_id()[0] == ' ':
                    for a in res:
                        if co == 0:
                            min_ = (res.get_id()[1])
                            co += 1
                        else:
                            max_ = (res.get_id()[1])
            pdbFile.close()
            rangeMINMAX.append(min_)
            rangeMINMAX.append(max_)
            if max_ - min_ > 100:
                templatesTMdata.remove(templatesTMdata[i])
            print templatesTMdata
    def create_new_temp_unit_gaps(self, unit, unitTemp, lostResiduesE, lostResiduesI):
        wdir,SRULdir, directory, templatesWork=self.wdir, self.SRULdir, self.directory, self.templatesWork
        dire = ''
        countAux = 0
        works = []
        if os.path.isfile(SRULdir+unit):
            dire = SRULdir
        elif os.path.isfile(wdir+templatesWork + 'units/' + unit):
            dire = wdir + templatesWork + 'units/'
        if dire != '':
            pdbFile = open(dire+unit, 'r')
            parser = PDBParser()
            co, count = 0, 0
            structureAux = parser.get_structure(unit[:4], pdbFile)
            for res in structureAux[0][unit[4:5]].get_residues():
                if co == 0:
                    groupNumI = int(res.get_id()[1])
                if res.get_id()[0] == ' ':
                    for a in res:
                        groupNumE = int(res.get_id()[1])
                if a.get_name() == 'C':
                    count += 1
            pdbFile.close()
            groupTocutE = groupNumE - int(lostResiduesE)
            groupTocutI = groupNumI + lostResiduesI - 1
            it = ''
            countAux = self.create_aux_unit(groupTocutI, groupTocutE+1, dire, wdir + templatesWork + "units/", unit, it, unitTemp)
            if countAux > self.minLenTemplate*self.porcMinLenUnit:
                works.append(1)
                return works, countAux
            else:
                if os.path.isfile(wdir + templatesWork + "units/" + unitTemp):
                    os.remove(wdir + templatesWork + "units/" + unitTemp)
                works.append(0)
                return works, countAux
        else:
            works.append(0)
            return works, countAux

    def create_new_temp_unit_gaps_part(self, unit, unitTemp, lostResiduesE, lostResiduesI, directo, templatesWork) :
        dire = ''
        countAux = 0
        works = []
        if os.path.isfile(directo+templatesWork + 'units/' + unit):
            dire = directo + templatesWork + 'units/'
        if dire != '':
            pdbFile = open(dire+unit, 'r')
            parser = PDBParser()
            co = 0
            structureAux = parser.get_structure(unit[:4], pdbFile)
            for res in structureAux[0][unit[4:5]].get_residues():
                if co == 0:
                    groupNumI = int(res.get_id()[1])
                    co += 1
                if res.get_id()[0] == ' ':
                    for a in res:
                        groupNumE = int(res.get_id()[1])
            groupTocutE = groupNumE - int(lostResiduesE)
            groupTocutI = groupNumI + int(lostResiduesI)
            pdbFile.close()
            it = ''
            countAux = self.create_aux_unit(groupTocutI, groupTocutE, dire, directo + templatesWork + "units/", unit, it, unitTemp)
            if countAux > self.minLenTemplate*self.porcMinLenUnit:
                works.append(1)
                return works, countAux
            else:
                if os.path.isfile(directo + templatesWork + "units/" + unitTemp):
                    os.remove(directo + templatesWork + "units/" + unitTemp)
                works.append(0)
                return works, countAux
        else:
            works.append(0)
            return works, countAux

    def create_pdb_fragments(self, rangeMINMAX, j, structure, it, targetsListOpt, templateListOpt, fastpList, target, targetOriginal, InfoRanges):
        wdir, directory, templatesWork = self.wdir, self.directory, self.templatesWork
        spp, spc = ' ', ''
        k = 1
        min_ = int(rangeMINMAX[0])
        max_ = int(rangeMINMAX[1])
        fastpFile = open(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.fastp', 'w')
        if os.path.isfile(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.mus'):
            mustangFile = open(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.mus', 'a')
        else:
            mustangFile = open(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.mus', 'w')
            mustangFile.write(">" + self.OUTdir + "\n")
        lenFile = open(wdir + templatesWork + targetOriginal[:5] + '/' + targetOriginal[:5] + '.len', 'w')
        fpdbUnit = open(wdir + templatesWork + 'units/' + target[:5] + 'repeat' + str(it), 'w')
        countLen, countLeft, countRight = 0, 0, 0
        fpdbUnitPre = open(wdir + templatesWork + target[:5] + "/" + target + 'pre.pdb', 'w')
        fpdbUnitPost = open(wdir + templatesWork + target[:5] + "/" + target + 'post.pdb', 'w')
        for res in structure[0][target[4:5]].get_residues():
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
                    if int(res.get_id()[1]) in range(int(1), int(min_)):
                        if a.get_name() == 'C':
                                countLeft += 1
                        if k >= 10000:
                                fpdbUnitPre.write('ATOM  ' + '%5s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id( )[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) +'  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        else:
                                fpdbUnitPre.write('ATOM   ' + '%4s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id( )[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')

                    elif int(res.get_id()[1]) in range(int(min_), int(max_)+1):
                        if a.get_name() == 'C':
                            countLen += 1
                        if k >= 10000:
                                fpdbUnit.write('ATOM  ' + '%5s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id( )[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        else:
                                fpdbUnit.write('ATOM   ' + '%4s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id( )[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                    else:
                        if a.get_name() == 'C':
                            countRight += 1
                        if k >= 10000:
                                fpdbUnitPost.write('ATOM  ' + '%5s'%str(k) + '  ' + a.get_name() + spc+'%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id( )[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                        else:
                                fpdbUnitPost.write('ATOM   ' + '%4s'%str(k) + '  ' + a.get_name() + spc + '%3s'%res.get_resname() + ' ' + target[4:5] + '%4s'%str(res.get_id( )[1]) + '    ' + '%8s'%str('%.3f'%a.get_coord()[0]) + '%8s'%str('%.3f'%a.get_coord()[1]) + '%8s'%str('%.3f'%a.get_coord()[2]) + '  ' + str('%.2f'%a.get_occupancy()) + '%6s'%str('%.2f'%a.get_bfactor()) + spp + a.get_id()[0] + '\n')
                    k += 1
        fpdbUnitPre.close()
        fpdbUnitPost.close()
        fpdbUnit.close()
        if countLeft > 6:
            targetsListOpt.append(target + 'pre')
        else:
            os.remove(wdir+templatesWork + target[:5] + '/' + target + 'pre.pdb')
        if countRight > 6:
            targetsListOpt.append(target + 'post')
        else:
            os.remove(wdir+templatesWork + target[:5] + '/' + target + 'post.pdb')
        if countLen > 6:
            templateListOpt.append(target[:5] + 'repeat' + str(it))
            fastpFile.write(str(min_) + "\t" + str(max_) + "\n")
            mustangFile.write("+" + target[:5] + 'repeat' + str(it) + "\n")
            lenFile.write(str(countLen) + "\n")
            InfoRanges.append([target + 'repeat', int(min_), int(max_)])
            fastpList.append([min_, max_])
        else:
            if os.path.isfile(wdir + templatesWork + target[:5] + '/' + target[:5] + 'repeat' + str(it)):
                os.remove(wdir + templatesWork + target[:5] + '/' + target[:5] + 'repeat' + str(it))
            if os.path.isfile(wdir + templatesWork + "units/" + target[:5] + 'repeat' + str(it)):
                os.remove(wdir + templatesWork + "units/" + target[:5] + 'repeat' + str(it))
        fastpFile.close()
        lenFile.close()
        mustangFile.close()

    def is_contiguos(self, start, end, fastpList):
        for pair in fastpList:
            if (start - int(pair[1]) <= 3 and (start - int(pair[1]) > 0)) or ((int(pair[0]) - end <= 3) and (int(pair[0]) - end > 0)):
                return 0
        return 1

    def is_first_of_repeat(self, InfoRanges, first):
        for row in InfoRanges:
            line = row
            if line[0][len(line[0]) - 3:] == "eat" and line[1] == first:
                return True
        return False

    def is_after_of_repeat(self, InfoRanges, last):
        for row in InfoRanges:
            line = row
            if line[0][len(line[0]) - 3:] == "eat" and line[2] == last:
                return True
        return False

    def obtain_str_insertions(self,alitar,alitmp,rangeMINMAX,strins):
        prevcount = 0
        if alitmp.count('---') >= 1:
            auxaliins = alitmp
            auxaliinstar = alitar
            while auxaliins.count('---') >= 1:
                iniGaps = int(re.search("---", auxaliins).start())
                cantGapsTar = auxaliinstar[:iniGaps].count(" ")
                ini_insertion = rangeMINMAX[0] + iniGaps - cantGapsTar + prevcount
                if re.search("[A-Z]", auxaliins[iniGaps:]):
                    endgaps = (re.search("[A-Z]", auxaliins[iniGaps:]).start()) + iniGaps
                else:
                    endgaps = len(auxaliins[iniGaps:]) + iniGaps
                cantGapsTar2 = auxaliinstar[:endgaps].count(" ")
                end_insertion = endgaps + rangeMINMAX[0] - cantGapsTar2 - 1 + prevcount
                strins = strins + str(ini_insertion) + " " + str(end_insertion) + ";"
                auxaliins = auxaliins[endgaps:]
                auxaliinstar = auxaliinstar[endgaps:]
                prevcount = endgaps
        return strins

    def put_all_templates_tmdata_particular(self, templatesTMdata, templatesTMdataSecondRound, templatesTMdataThirdRound, templatesTMdataLowerRound, particular):
        if not particular:
            templatesTMdata = (templatesTMdata + templatesTMdataSecondRound + templatesTMdataThirdRound + templatesTMdataLowerRound)
        lst = []
        i = 0

        while i < len(templatesTMdata):
                if particular:
                    start = templatesTMdata[i][2]
                else:
                    start = templatesTMdata[i][10]
                if (start in lst) or (start + 1 in lst) or (start + 2 in lst) or (start + 3 in lst) or (start - 1 in lst) or (start - 2 in lst) or (start - 3 in lst):
                    templatesTMdata.remove(templatesTMdata[i])
                else:
                    lst.append(start)
                    i += 1
        return templatesTMdata

