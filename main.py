from fastapi import FastAPI, File, UploadFile
from fastapi.responses import FileResponse
from typing import Optional
from pydantic import BaseModel
import os
import pfco
import wholesequence

app = FastAPI()
path = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/pfam_data/'


class Protein(BaseModel):
    pf_code: str
    path: str


# Para descargar una familia de proteínas con pfco
@app.get("/getalign/")
def getalign(path: str, pf_code: str):
    file_path = pfco.main(path, pf_code)
    return FileResponse(file_path)

# Para ejecutar la funcion wholeseq


@app.get("/wholeseq/")
def wholeseq(path: str, pf_code: str, seed: str):
    file_path = wholesequence.main(path, pf_code)
    return FileResponse(file_path)
