from fastapi import FastAPI, File, UploadFile
from fastapi.responses import FileResponse
from typing import Optional
from pydantic import BaseModel
import os
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

app = FastAPI()
base = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/'


class Protein(BaseModel):
    pf_code: str
    dir: str


# Para descargar una familia de proteínas con pfco
@app.get("/getalign/")
def getalign(path: str, pf_code: str):
    file_path = pfco.main(path, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such ALIGN.fasta file"}


# Para ejecutar la funcion wholeseq
@app.get("/wholeseq/")
def wholeseq(path: str, pf_code: str, seed: str):
    file_path = wholesequence.main(path, pf_code, seed)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such wholeALIGN.fasta file"}


# Para retornar el archivo efetch
@app.get("/efetch/")
def efetch(path: str, pf_code: str, accnumber: str):
    file_path = efetch.main(path, pf_code, accnumber)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such efetch.txt file"}


# Para retornar el archivo DESC
@app.get("/desc/")
def desc(path: str, pf_code: str):
    file_path = desc.main(path, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such DESC file"}


# Para ejecutar la funcion createalign
@app.get("/createalign/")
def createalign(path: str, pf_code: str, seed: str, method: str):
    file_path = createalign.main(path, pf_code, seed, method)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such newALIGN.fasta file"}


# Para ejecutar la funcion extend
@app.get("/extend/")
def extend(path: str, pf_code: str, seed: str, n: str, c: str, method: str):
    file_path = extendterminal.main(path, pf_code, seed, n, c, method)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such extendALIGN.fasta file"}


# Para ejecutar la funcion pfnew: PARA PROTEÍNAS SIN FAMILIA
@app.get("/pfnew/")
def pfnew(path: str, accnumdir: str):
    file_path = pfnew.main(path, accnumdir)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"pfnew error"}


# Para ejecutar la funcion pfbuild
@app.get("/pfbuild/")
def pfbuild(path: str, pf_code: str):
    file_path = pfbuild.main(path, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such HMM file"}


# Para ejecutar la funcion pfmake
@app.get("/pfmake/")
def pfmake(path: str, pf_code: str, tmayus: str, tminus: str, e: str):
    file_path = pfmake.main(path, pf_code, tmayus, tminus, e)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such ALIGN.fasta file"}


# Para ejecutar la funcion nextduf
@app.get("/nextduf/")
def nextduf(path: str, pf_code: str):
    file_path = nextduf.main(path, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such nextDUF.txt file"}


# Para ejecutar la funcion pfci: PARA PROTEÍNAS CON FAMILIA
@app.get("/pfci/")
def pfci(path: str, pf_code: str, options: str, description: str):
    answer = pfci.main(path, pf_code, options, description)
    return answer


# Para ejecutar la funcion missing
@app.get("/missing/")
def missing(path: str, pf_code: str):
    result = missing.main(path, pf_code)
    if os.path.exists(result["missing_path"]):
        missing_file = FileResponse(result["missing_path"])
    else:
        missing_file = "There's no such missing file"
    if os.path.exists(result["found_path"]):
        found_file = FileResponse(result["found_path"])
    else:
        found_file = "There's no such found file"

    return {"result": result["result"], "missing_file": missing_file, "found_file": found_file}


# Para ejecutar la funcion overlap
@app.get("/overlap/")
def overlap(path: str, pf_code: str):
    file_path = overlap.main(path, pf_code)
    if os.path.exists(file_path):
        return FileResponse(file_path)
    else:
        return {"error": "There's no such overlap file"}


# Para ejecutar la funcion pfamoutput
@app.get("/pfamoutput/")
def pfamoutput(path: str, pf_code: str, evalue: float):
    result = pfamoutput.main(path, pf_code, evalue)
    return {"cutoffvalues": result["cutoffvalues"], "sequences": result["sequences"], "domains": result["domains"]}
