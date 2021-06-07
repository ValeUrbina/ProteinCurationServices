# SPCleaner

## Flow

#### pfam_download/

1. **pfam_download.py**
   - Input example: PF03377
2. WHOLE SEQUENCE FALTA

#### curation/

(applied to non repeats alignments)
PASO PREVIO: HACER WHOLESEQ PARA OBTENER LA SECUENCIA COMPLETA DEL ALINEAMIENTO DE LA FAMILIA DESCARGADA

1. **divide_alignments.py** -> divide segmentos de la secuencia (unidad de repetición). hace un mappeo de secuencia pfam a uniprot, realiza cortes y crea un nuevo alineamiento que corresponde a una unidad de repetición.
   - test input file:
     - output for previous step
     - ALIGN
2. **stockholm_to_fasta.py** -> cambio de formato Stockholm a FASTA
   - test input file:
     - output for previous step
     - ALIGN_repeat
3. **mafft_align.py** -> aplica el alineamieto MAFT para conservar el alineamiento
   - test input file:
     - output for previous step
     - prueba_fasta
4. **fasta_to_stockholm.py** -> Regresa al formato original
   - test input file:
     - output for previous step
     - alineado_fasta

(applied to repeat and non repeat alignments, I think...)

5. **alignment_curator.py** -> incrementa la conservación del alineamiento teniendo en cuenta las propiedades que tiene como unidad de repetición usando la matriz de blossum y los parámetros (para que se vea bien el alineamiento) elimina ciertas secuencias

   - test input file:
     - output for previous step
     - alineado_stockholm

   -> PARA TRABAJOS FUTUROS: CAMBIAR LOS PERCENTILES: PERCENTIL MÍNIMO Y MÁXIMO DE LAS LONGITUDES DE LAS UNIDADES DE REPETICIÓN PARA UNIFORMIZAR ESTAS UNIDADES, EN DIVIDE_ALIGNMENT: (LINEA 50) if len(short_seq) < 5: LONGITUD MÍNIMA PERMITIDA DE UNA UNIDAD DE REPETICIÓN

-> PARA TRABAJOS FUTUSO: QUE LA APARIENCIA SEA CONFIGURABLE

#### TRIM/

1. trim.py -> como a veces quedan huecos o vacíos en las columnas de las secuencias, el programa se encarga de eliminar la secuencia que
   provoca esos vacíos para mejorar el alineamiento

#### hmm_maker/

1. **hmm_profile_maker.py** -> crea un perfil de markov para el alineamiento
   - test input file:
     - output for previous step
     - test
2. **hmm_logo_maker.py** -> a partir del perfil crea el logo
   - test input file:
     - output for previous step
     - hmm_logo_maker_input
