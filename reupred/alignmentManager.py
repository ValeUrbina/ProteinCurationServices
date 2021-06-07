import time
import Bio
from threading import Thread
import argparse
import os
import shutil
import subprocess
import otherMethods
import re
from Bio.PDB import *
from operator import itemgetter
import logging

class Alignment:
    #'Common base class for alignments'
    targetCount = 0

    def __init__(self, row):
        target, template, start, end, gapsTar, gapsTem, nicePos, len1, len2, aliLen, RMSD, seqID, TMscoreUnit1, TMscoreUnit2, TMscoreUnitAvg, SeqUnit1, SeqRel, SeqUnit2, coverage = '', '', 0, 0, 0.0, 0.0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, "", "", "", 0.0
        if row and (len(row)) > 20:
            rowTemp = row.split('\t')
            if len(rowTemp) > 17:
                target, template, start, end, gapsTar, gapsTem, nicePos, len1, len2, aliLen, RMSD, seqID, TMscoreUnit1, TMscoreUnit2, TMscoreUnitAvg, SeqUnit1, SeqRel, SeqUnit2, coverage = rowTemp
        self.target, self.template, self.unitOlength, self.start, self.end, self.nicePos, self.unitDlength, self.alignmentLength, self.alignmentRMSD, self.sequenceIdentity, self.tmScoreO, self.tmScoreD, self.tmScoreAVG, self.sequenceO, self.sequenceRelation, self.sequenceD, self.gapsTarget, self.gapsTemplate, self.coverage = target, template, len1, start, end, nicePos, len2, aliLen, RMSD, seqID, TMscoreUnit1, TMscoreUnit2, TMscoreUnitAvg, SeqUnit1, SeqRel, SeqUnit2, gapsTar, gapsTem, coverage

    def get_average_tmscore(self):

        return float(self.tmScoreD)

    def get_coverage(self):
        return float(self.coverage)

    def get_info_align(self):
        return self.target, self.template, self.unitOlength, self.start, self.end, self.nicePos, self.unitDlength, self.alignmentLength, self.alignmentRMSD, self.sequenceIdentity, self.tmScoreO, self.tmScoreD, self.tmScoreAVG, self.sequenceO, self.sequenceRelation, self.sequenceD, self.gapsTarget, self.gapsTemplate, self.coverage

    def get_rmsd(self):
        return float(self.alignmentRMSD)

    def display_alignments(self):
        print self.target, self.template, self.unitOlength, self.start, self.end, self.nicePos, self.unitDlength, self.alignmentLength, self.alignmentRMSD, self.sequenceIdentity, self.tmScoreO, self.tmScoreD, self.tmScoreAVG, self.sequenceO, self.sequenceRelation, self.sequenceD, self.gapsTarget, self.gapsTemplate, self.coverage

    def get_template(self):
        return self.template


class ProteinAlignmentSet:

    def __init__(self, prot, DSSPdirectory, pdbFolder,Dsspexe,TMalignexe):
        self.protein = prot
        self.Alignments = []
        self.Dsspexe = Dsspexe
        self.TMalignexe = TMalignexe
        #self.dssp = self.calculate_dssp(DSSPdirectory, pdbFolder)
       # print "xxx",pdbFolder+prot

        self.dssp = self.aligned_dssp_protein(pdbFolder+prot,False)
        #for dsspvals in list(dssp.keys()):
         #   secstr = secstr + dssp[dsspvals][2]
        self.maxTM = 0.0
        self.minRMSD = 10.0

    def aligned_dssp_protein(self,filename, dic):
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
            dssp = DSSP(model, filename + '.pdb',dssp=self.Dsspexe)
            for dsspvals in list(dssp.keys()):
                secstr = secstr + dssp[dsspvals][2]
            if dic == True:
                return dssp
            else:
                return secstr
        except Exception, e:
            logging.error("Error while trying to obtain the dssp of the aligned units | " + filename + ".pdb" + str(e))
            return ''

    def add_data_to_templates_tmdata_align(self, tmOutput):
            templatesTMdata = []
            dataList = tmOutput.split('\n')
            lines = tmOutput.split('\n')
            aliLen = lines[16].split(',')[0].split("=")[1]
            RMSD = lines[16].split(',')[1].split("=")[1]
            SeqId = float(lines[16].split(',')[2].split("=")[2])
            lenline2 = lines[14].split(':')[1].split("r")[0]
            tmscore2 = float(lines[18].split(' ')[1])
            tmscore1 = float(lines[17].split(' ')[1])
            tmscore3 = float(lines[19].split(' ')[1])
            align = lines[23]
            align2 = lines[24]
            align3 = lines[25]
            if re.search('[A-Z]', align3):
                start = re.search('[A-Z]', align3).start()

            else:
                start = 0
            if re.search("-+$", align3):
                end = re.search("-+$", align3).start()

            else:
                end = len(align3)
            gaps = align[start:end].count('-')
            gapsUnit = align3[start:end].count('-')
            gapsTem = float((float(gapsUnit) / float(int(end) - int(start))))
            gapsTar = float((float(gaps) / float(int(end) - int(start))))
            coverage,contig = 0,0
            aliLen_ = end - start + 1
            diffEnd = len(align) - end
            if start < diffEnd and start < 4:
                    nicePosition = 1
            else:
                if start > diffEnd and diffEnd < 4:
                    nicePosition = 1
                else:
                    nicePosition = 0
            lenTarget = len(align[int(start):int(end)]) - gaps
            lenUnit = int(end) - int(start) - int(gapsUnit) - int(gaps)
            if aliLen_ >= gaps and (aliLen_ >= gapsUnit):
                coverage = float(abs(float(2 * int(aliLen)) / float((int(lenline2) + int(lenTarget)))))
            return start, end, gapsTar, gapsTem, nicePosition, lenTarget, lenUnit, aliLen, float(RMSD), SeqId, tmscore1, tmscore2, tmscore1, align, align2, align3, coverage

    def calculate_dssp(self, DSSPdirectory, pdbFolder):
        entry = self.protein
        if not os.path.isfile(DSSPdirectory + entry + '.dssp'):
            os.system(self.Dsspexe+' -i ' + pdbFolder+entry[:5] + '.pdb -o ' + DSSPdirectory + entry + '.dssp')
        if os.path.isfile(DSSPdirectory + entry + '.dssp'):
            with open(DSSPdirectory + entry + '.dssp', 'r') as f:
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
        return dssp

    def display(self):
        for row in self.Alignments:
            print row

    def get_alignments(self):
        return  self.Alignments

    def get_alignments_name(self):
        listnames = []
        for ali in self.Alignments:
            listnames.append(ali.get_template())
        return listnames

    def get_max_tmscore(self):
        return self.maxTM

    def calc_max_tmscore(self):
        auxList = []
        sum = 0
        for ali in self.Alignments:

            auxList.append(ali.get_average_tmscore())
            sum = ali.get_average_tmscore()+sum
        if float(len(self.Alignments)) == 0:
            return (auxList), 0.0
        else:
            return (auxList), float(float(sum) / float(len(self.Alignments)))

    def calc_max_coverage(self):
        auxList = []
        sum = 0
        for ali in self.Alignments:
            auxList.append(ali.get_coverage())
            sum = ali.get_coverage() + sum
        return (auxList), float(float(sum) / float(len(self.Alignments)))

    def get_min_rmsd(self):
        return self.minRMSD

    def calc_min_rmsd(self):
        auxList = []
        for ali in self.Alignments:
            auxList.append(ali.get_rmsd())
        if auxList:
            return min(auxList)
        else:
            return 999

    def get_all_rmsd(self, auxList):
        for ali in self.Alignments:
            auxList.append(ali.get_rmsd())
        return auxList

    def get_all_tmscore(self, auxList):
        for ali in self.Alignments:
            auxList.append(ali.get_average_tmscore())
        return auxList

    def load_alignments(self, templatesList, pdbFolder, alignFolder, SRULdirectory, sc,exeTMalign):
        AlignmentsData = []
        entry = self.protein

        if entry and ((not os.path.isfile(alignFolder + sc + "AliTMalign" + entry[:5])) or os.stat(alignFolder + sc + "AliTMalign" + entry[:5]).st_size < 1):
            t = Thread(target=self.calc_alignments, args=(templatesList, pdbFolder, entry, SRULdirectory, AlignmentsData, alignFolder, sc,exeTMalign))
            t.start()
            while t.is_alive():
                pass
        if os.path.isfile(alignFolder + sc + "AliTMalign" + entry[:5]):
            with open(alignFolder + sc + "AliTMalign" + entry[:5], 'r') as f:
                AlignmentsData = f.read().split("\n")
            f.close()
        self.Alignments = [Alignment(row) for row in AlignmentsData if row]
        self.maxTM = self.calc_max_tmscore()
        self.minRMSD = self.calc_min_rmsd()

    def calc_alignments(self, templatesList, pdbFolder, entry, SRULdirectory, AlignmentsData, alignFolder, sc,TMalignexe):
        logging.debug("Calculating alignments")
        logfile = open(alignFolder + sc + "xxxAliTMalign" + entry[:5], 'w')
        AliTMalignFile = open(alignFolder + sc + "AliTMalign" + entry[:5], 'w')
        if sc == 'III':
            min_value_coverage = 0.8
            min_value_rmsd=2.3
            minTMscore = 0.35
        elif sc == 'IV':
            min_value_coverage = 0.6
            min_value_rmsd = 2.8
            minTMscore = 0.4
        elif sc == 'V':
            min_value_coverage = 0.6
            min_value_rmsd = 3.5
            minTMscore = 0.4
        processes = []
        for file in templatesList:
            f = os.tmpfile()

            bashCommand = self.TMalignexe +' ' + pdbFolder + entry[:5] + '.pdb ' + SRULdirectory + file + ' -a T'
            p = subprocess.Popen(bashCommand, shell=True, stdout=subprocess.PIPE)
            tmOutput = p.communicate()[0]
            processes.append((p, tmOutput, file))
        for p, tmOutput,unit in processes:
            p.wait()
            if len(tmOutput) > 100:
                lines = tmOutput.split('\n')
              #  print sc, float(lines[18].split(' ')[1]), minTMscore
                if float(lines[18].split(' ')[1]) >= minTMscore:
                    startA, endA, gapsTar, gapsTem, nicePos, len1, len2, aliLen, RMSD, seqID, TMscoreUnit1, TMscoreUnit2, TMscoreUnitAvg, SeqUnit1, SeqRel, SeqUnit2, coverage = self.add_data_to_templates_tmdata_align(tmOutput)
                    #if sc == 'IV':

                    start = startA - SeqUnit1[:startA].count("-")
                    end = endA - SeqUnit1[:endA].count("-")

                    if (sc =='V' and coverage > 0 and   float(coverage) >= float(min_value_coverage) and float(RMSD) <= float(min_value_rmsd) )or((sc=='III' or sc=='IV') and coverage > 0  and float(coverage) >= float(min_value_coverage) and float(RMSD) <= float(min_value_rmsd) ):#and \
                           # ((start >= 2 and (self.dssp[start - 2] != 'E' and self.dssp[start - 2] != 'H'))or (start >= 1 and (self.dssp[start - 1] != 'E' and self.dssp[start - 1] != 'H')) or (self.dssp[start] != 'E' and self.dssp[start] != 'H') or start == 0 or (self.dssp[start + 1] != 'E' and self.dssp[start + 1] != 'H') or (self.dssp[start + 2] != 'E' and self.dssp[start + 2] != 'H')) and \
                            #(((end + 1) < len(SeqUnit1) and (self.dssp[end + 1] != 'E' and self.dssp[end + 1] != 'H')) or (self.dssp[end] != 'E' and self.dssp[end] != 'H') or end == len(self.dssp) or (self.dssp[end - 1] != 'E' and self.dssp[end - 1] != 'H') or (self.dssp[end - 2] != 'E' and self.dssp[end - 2] != 'H') ) ):
                        if self.minRMSD > RMSD:
                            self.minRMSD = RMSD
                        if self.maxTM < TMscoreUnit2:
                            self.maxTM = TMscoreUnit2
                        AlignmentsData.append([entry, unit, startA, endA, gapsTar, gapsTem, nicePos, len1, len2, aliLen, RMSD, seqID, TMscoreUnit1, TMscoreUnit2, TMscoreUnitAvg, SeqUnit1, SeqRel, SeqUnit2, coverage])
        for line in AlignmentsData:
            entry, unit, startA, endA, gapsTar, gapsTem, nicePos, len1, len2, aliLen, RMSD, seqID, TMscoreUnit1, TMscoreUnit2, TMscoreUnitAvg, SeqUnit1, SeqRel, SeqUnit2, coverage = line
            AliTMalignFile.write(
            entry + "\t" + unit + "\t" + str(startA) + "\t" + str(endA) + "\t" + str(gapsTar) + "\t" + str(
                gapsTem) + "\t" + str(nicePos) + "\t" + str(len1) + "\t" + str(len2) + "\t" + str(aliLen) + "\t" + str(
                RMSD) + "\t" + str(seqID) + "\t" + str(TMscoreUnit1) + "\t" + str(TMscoreUnit2) + "\t" + str(
                TMscoreUnitAvg) + "\t" + SeqUnit1 + "\t" + SeqRel + "\t" + SeqUnit2 + "\t" + str(coverage) + "\n")
        AliTMalignFile.close()

