## Interpretação de Variantes

```python
MeuDrive = "/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"
amostra = "cap-ngse-b-2019-chr11"

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

caminho_arquivo = f"{MeuDrive}/{amostra}/output/{amostra}.clinvar.ann.txt"
df = pd.read_csv(caminho_arquivo, sep="\t")
df
```

![image](https://github.com/user-attachments/assets/ebe87341-131e-47a0-a197-496a6497782e)

```python
patho = df['CLNSIG'].isin(["Pathogenic", "Likely_pathogenic", "Uncertain_significance"])
df[patho]
```

![image](https://github.com/user-attachments/assets/b1ec22e4-f5ca-450a-879e-1326e984a737)

Após a anotação de variantes, com auxílio de ferramentas como o Pandas, para melhor visualização, filtramos o resultado final por "CLNSIG" (entrada no Clinvar) "Pathogenic", "Likely_pathogenic", "Uncertain_significance", para que todas as variantes com entrada e sem entrada fossem chamadas para análise manual. 

Variantes encontradas: 7

Durante análise manual, observamos majoritariamente os campos "POS" (posição cromossômica), "CLNDN" (entrada no OMIM), "CLNSIG" (entrada no ClinVar), "GENEINFO" (nome do gene) e "ANN" para extrair informações como HGVS, transcrito, etc. e chegamos no seguinte resultado:


