import os
import to_fasta
import sendemail
import newproject
import pfco
import deletesequences
import deletefragments
import wholesequence
import efetch
import desc
import createalign
import extendterminal
import pfnew
import pfbuild
import pfmake
import nextduf
import pfci
import missing
import overlap
import pfamoutput
import repeatsdb
import hmmer
import deleteleft
import deleteright
from spcleaner.alignment_curator import main as spcleaner
from fastapi import FastAPI
from fastapi.responses import FileResponse
from pydantic import BaseModel

app = FastAPI()
base = '/home/ubuntu/tesis/pfam_curation/'
pfamseq_path = '/home/ubuntu/tesis/pfam_curation/seqlib/pfamseq'
#base = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/'
#pfamseq_path = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/seqlib/pfamseq'


class List(BaseModel):
    acc_numbers_ids: list


# Para crear un nuevo proyecto
@app.get("/newproject/")
def newproject_service(curator_id: str, project_name: str, email: str, pfam_code: str):
    resault = newproject.main(base, curator_id, project_name, email, pfam_code)
    if os.path.exists(resault["project_path"]):
        return resault["directory"]
    else:
        return {"error": "Error in creating the directory"}


# Para descargar una familia de proteínas con pfco
@app.get("/getalign/")
def getalign_service(project_name: str, pf_code: str):
    file_path = pfco.main(base, project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such ALIGN.fasta file"}


# Para ejecutar la funcion deleteleft
@app.get("/deleteleft/")
def deleteleft_service(project_name: str, pf_code: str, column: int, file_name: str, outputfile_name: str):
    file_path = deleteleft.main(
        base, project_name, pf_code, column, file_name, outputfile_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such delseqlefttALIGN.fasta file"}


# Para ejecutar la funcion deleteright
@app.get("/deleteright/")
def deleteright_service(project_name: str, pf_code: str, column: int, file_name: str, outputfile_name: str):
    file_path = deleteright.main(
        base, project_name, pf_code, column, file_name, outputfile_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such delseqrightALIGN.fasta file"}


# Para ejecutar la funcion deletesequences
@app.get("/deletesequences/")
def deletesequences_service(project_name: str, pf_code: str, accnumbids_list: List, file_name: str, outputfile_name: str):
    file_path = deletesequences.main(
        base, project_name, pf_code, accnumbids_list.acc_numbers_ids, file_name, outputfile_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such delseqALIGN.fasta file"}


# Para ejecutar la funcion deletefragments
@app.get("/deletefragments/")
def deletefragments_service(project_name: str, pf_code: str, file_name: str, outputfile_name: str):
    file_path = deletefragments.main(
        base, project_name, pf_code, file_name, outputfile_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such delseqALIGN.fasta file"}


# Para ejecutar la funcion wholeseq
@app.get("/wholeseq/")
def wholeseq_service(project_name: str, pf_code: str, file_name: str, outputfile_name: str):
    file_path = wholesequence.main(
        base, project_name, pf_code, file_name, outputfile_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such wholeALIGN.fasta file"}


# Para retornar el archivo efetch
@app.get("/efetch/")
def efetch_service(project_name: str, pf_code: str, accnumber: str):
    file_path = efetch.main(base, project_name, pf_code, accnumber)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such efetch.txt file"}


# Para retornar el archivo DESC
@app.get("/desc/")
def desc_service(project_name: str, pf_code: str):
    file_path = desc.main(base, project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such DESC file"}


# Para ejecutar la funcion createalign
@app.get("/createalign/")
def createalign_service(project_name: str, pf_code: str, method: str, file_name: str, outputfile_name: str):
    file_path = createalign.main(
        base, project_name, pf_code, method, file_name, outputfile_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such newALIGN.fasta file"}


# Para ejecutar la funcion extend
@app.get("/extend/")
def extend_service(project_name: str, pf_code: str, n: str, c: str, method: str, file_name: str, outputfile_name: str):
    file_path = extendterminal.main(
        base, project_name, pf_code, n, c, method, file_name, outputfile_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such extendALIGN.fasta file"}


# Para ejecutar la funcion pfnew: PARA PROTEÍNAS SIN FAMILIA
@app.get("/pfnew/")
def pfnew_service(project_name: str, accnumdir: str):
    file_path = pfnew.main(base, project_name, accnumdir)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"pfnew error"}


# Para ejecutar la funcion pfbuild
@app.get("/pfbuild/")
def pfbuild_service(project_name: str, pf_code: str):
    file_path = pfbuild.main(
        base, project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such HMM file"}


# Para ejecutar la funcion pfmake
@app.get("/pfmake/")
def pfmake_service(project_name: str, pf_code: str, tmayus: str, tminus: str, e: str, outputfile_name: str):
    file_path = pfmake.main(base, project_name, pf_code,
                            tmayus, tminus, e, outputfile_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such ALIGN.fasta file"}


# Para ejecutar la funcion nextduf
@app.get("/nextduf/")
def nextduf_service(project_name: str, pf_code: str):
    file_path = nextduf.main(base, project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such nextDUF.txt file"}


# Para ejecutar la funcion pfci: PARA PROTEÍNAS CON FAMILIA
@app.get("/pfci/")
def pfci_service(project_name: str, pf_code: str, options: str, description: str):
    answer = pfci.main(base, project_name, pf_code, options, description)
    return answer


# Para ejecutar la funcion missing
@app.get("/missing/")
def missing_service(project_name: str, pf_code: str):
    result = missing.main(base, project_name, pf_code)
    return result


# Para ejecutar la funcion overlap
@app.get("/overlap/")
def overlap_service(project_name: str, pf_code: str):
    file_path = overlap.main(base, project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such overlap file"}


# Para ejecutar la funcion pfamoutput
@app.get("/pfamoutput/")
def pfamoutput_service(project_name: str, pf_code: str, evalue: float):
    result = pfamoutput.main(base, project_name, pf_code, evalue)
    return {"cutoffvalues": result["cutoffvalues"], "sequences": result["sequences"], "domains": result["domains"]}


# Para ejecutar la funcion spcleaner
@app.get("/spcleaner/")
def spcleaner_service(pfam_code: str, project_name: str, file_input_name: str, file_output_name: str):
    file_path = os.path.join(base, 'pfam_data', project_name, pfam_code)
    output_path = spcleaner(pfam_code, file_path,
                            file_input_name, file_output_name)
    if os.path.exists(output_path):
        fasta_path = os.path.join(
            base, 'pfam_data', project_name, pfam_code, file_output_name + '.fasta')
        to_fasta.main(output_path, fasta_path)
    else:
        return {"error": "There's no such spcleaner file"}
    if os.path.exists(fasta_path):
        return FileResponse(fasta_path)
    else:
        return {"error": "There's no such spcleaner file"}


# Para reupred
@app.get("/reupred/")
def reupred_service():
    file_path = 'reupred/4gg4_A.reupred'
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no units found"}


# Para repeatsdb
@app.get("/repeatsdb/")
def repeatsdb_service(directory: str, pfam_code: str, pdb_id: str):
    resault = repeatsdb.main(base, directory, pfam_code, pdb_id)
    return resault


# Para hmmer
@app.get("/hmmer/")
def hmmer_service(project_name: str, pf_code: str, file_name: str):
    file_path = hmmer.main(base, pfamseq_path, project_name,
                           pf_code, file_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such domtblout.txt file"}


# Para el retornar cualquier archivo
@app.get("/returnfile/")
def returnfile(project_name: str, pfam_code: str, file_name: str):
    file_path = os.path.join(
        base, 'pfam_data', project_name, pfam_code, file_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no " + file_name + 'file'}


# Para enviar un email
@app.get("/sendemail/")
def sendemail_service(project_name: str, operation_id: str, receiver_mail: str, pfam_code: str):
    resault = sendemail.main(
        base, operation_id, receiver_mail, project_name, pfam_code)
    return resault
