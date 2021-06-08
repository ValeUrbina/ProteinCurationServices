import sys
import os
from datetime import datetime

# create a HMM over the existing alignment
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT = "mkdir /home/valeria/Documentos/Tesis_2/Docker/pfam_curation/pfam_data/{}/"


def main(PATH, curator_id, project_name, email, pfam_code):

    # fisrt create the project id
    now = datetime.now()
    timestamp = datetime.timestamp(now)
    project_id = round(timestamp)

    # then delete the white spaces from the project name
    project_name = project_name.replace(" ", "")

    # now we create the directory name of the project
    directory = project_name + '-' + curator_id + \
        '-' + pfam_code + '-' + str(project_id)

    # se crea la carpeta si no existe
    try:
        os.mkdir(os.path.join(PATH, 'pfam_data', directory))
        project_path = os.path.join(PATH, 'pfam_data', directory)
    except:
        directory = "There is already a proyect with the same name"

    return {"project_path": project_path, "directory": directory}


if __name__ == "__main__":
    main()
