from fastapi import FastAPI
from fastapi.responses import FileResponse
from pydantic import BaseModel
from spcleaner.alignment_curator import main as spcleaner
import os
import to_fasta
import newproject
import pfco
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

app = FastAPI()
base = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/'


class Protein(BaseModel):
    pf_code: str
    dir: str


# Para el undo y redo
@app.get("/returnfile/")
def returnfile(project_name: str, pfam_code: str, file_name: str):
    file_path = os.path.join(
        base, 'pfam_data', project_name, pfam_code, file_name)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no " + file_name + 'file'}


# Para descargar una familia de proteínas con pfco
@app.get("/newproject/")
def newproject_service(curator_id: str, project_name: str, email: str, pfam_code: str):
    resault = newproject.main(curator_id, project_name, email, pfam_code)
    if os.path.exists(resault["project_path"]):
        return resault["directory"]
    else:
        return {"error": "Error in creating the directory"}


# Para descargar una familia de proteínas con pfco
@app.get("/getalign/")
def getalign_service(project_name: str, pf_code: str):
    file_path = pfco.main(project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such ALIGN.fasta file"}


# Para ejecutar la funcion wholeseq
@app.get("/wholeseq/")
def wholeseq_service(project_name: str, pf_code: str, seed: str):
    file_path = wholesequence.main(project_name, pf_code, seed)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such wholeALIGN.fasta file"}


# Para retornar el archivo efetch
@app.get("/efetch/")
def efetch_service(project_name: str, pf_code: str, accnumber: str):
    file_path = efetch.main(project_name, pf_code, accnumber)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such efetch.txt file"}


# Para retornar el archivo DESC
@app.get("/desc/")
def desc_service(project_name: str, pf_code: str):
    file_path = desc.main(project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such DESC file"}


# Para ejecutar la funcion createalign
@app.get("/createalign/")
def createalign_service(project_name: str, pf_code: str, seed: str, method: str):
    file_path = createalign.main(project_name, pf_code, seed, method)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such newALIGN.fasta file"}


# Para ejecutar la funcion extend
@app.get("/extend/")
def extend_service(project_name: str, pf_code: str, seed: str, n: str, c: str, method: str):
    file_path = extendterminal.main(project_name, pf_code, seed, n, c, method)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such extendALIGN.fasta file"}


# Para ejecutar la funcion pfnew: PARA PROTEÍNAS SIN FAMILIA
@app.get("/pfnew/")
def pfnew_service(project_name: str, accnumdir: str):
    file_path = pfnew.main(project_name, accnumdir)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"pfnew error"}


# Para ejecutar la funcion pfbuild
@app.get("/pfbuild/")
def pfbuild_service(project_name: str, pf_code: str):
    file_path = pfbuild.main(project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such HMM file"}


# Para ejecutar la funcion pfmake
@app.get("/pfmake/")
def pfmake_service(project_name: str, pf_code: str, tmayus: str, tminus: str, e: str):
    file_path = pfmake.main(project_name, pf_code, tmayus, tminus, e)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such ALIGN.fasta file"}


# Para ejecutar la funcion nextduf
@app.get("/nextduf/")
def nextduf_service(project_name: str, pf_code: str):
    file_path = nextduf.main(project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such nextDUF.txt file"}


# Para ejecutar la funcion pfci: PARA PROTEÍNAS CON FAMILIA
@app.get("/pfci/")
def pfci_service(project_name: str, pf_code: str, options: str, description: str):
    answer = pfci.main(project_name, pf_code, options, description)
    return answer


# Para ejecutar la funcion missing
@app.get("/missing/")
def missing_service(project_name: str, pf_code: str):
    result = missing.main(project_name, pf_code)
    return result


# Para ejecutar la funcion overlap
@app.get("/overlap/")
def overlap_service(project_name: str, pf_code: str):
    file_path = overlap.main(project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such overlap file"}


# Para ejecutar la funcion pfamoutput
@app.get("/pfamoutput/")
def pfamoutput_service(project_name: str, pf_code: str, evalue: float):
    result = pfamoutput.main(project_name, pf_code, evalue)
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
    resault = repeatsdb.main(directory, pfam_code, pdb_id)
    return resault


# Para hmmer
@app.get("/hmmer/")
def hmmer_service(project_name: str, pf_code: str):
    file_path = hmmer.main(project_name, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such domtblout.txt file"}
