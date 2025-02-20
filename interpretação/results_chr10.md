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

Após a anotação de variantes, com auxílio de ferramentas como o Pandas, para melhor visualização, filtramos o resultado final por "CLNSIG" (entrada no Clinvar) "Pathogenic", "Likely_pathogenic", "Uncertain_significance", para que todas essas variantes fossem chamadas para análise manual. 

Variantes encontradas: 7

Durante análise manual, observamos majoritariamente os campos "POS" (posição cromossômica), "CLNDN" (entrada no OMIM), "CLNSIG" (entrada no ClinVar), "GENEINFO" (nome do gene) e "ANN" para extrair informações como HGVS, transcrito, etc. e chegamos no seguinte resultado:

## TBATA ✖️

HGVS: TBATA(NM_001318241.2):c.666T>C
p.(Ala222=)

chr10:70777180-70777180

- Frequência: 29% gnomAD Exome
- VAF: HET

Essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara.

## CHST3 ✖️

HGVS: CHST3(NM_004273.5):c.*66G>T
p.?

chr10:72008537-72008537

- Frequência: 0.14% gnomAD Exome
- VAF: HET

ACMG: Benigna
- BP7	Está em uma região UTR. (Variante sinônima ou não codificadora que não está localizada em uma região de splicing e não se prevê que tenha consequências que alterem o splicing). -2

## TMEM ✖️

HGVS: TMEM254(NM_001270367.1):c.378C>A
p.(Phe126Leu)

chr10:80090851-80090851

- Frequência: 0.019% gnomAD Exome
- VAF: HET

Não encontramos nenhuma patologia relacionada a esse gene nos banco de dados analisados (Varsome, Decipher, Clinvar). https://www.ncbi.nlm.nih.gov/clinvar/variation/3327028/

ACMG: VUS
- PM2	Frequência extremamente baixa em bancos de dados populacionais do gnomAD +1
- BP4	Múltiplas linhas de evidências computacionais sugerem nenhum impacto no gene ou no produto genético (conservação, evolução, impacto de splicing, etc.) -1

## NRG3 ✖️

HGVS: NRG3(NM_001370081.1):c.1986C>T
p.(Ser662=)

chr10:82985500-82985500

- Frequência: 32% gnomAD Exome
- VAF: HET

Essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara.

## NFKB2 ✔️

HGVS: NFKB2(NM_001322934.2):c.2557C>T
p.(Arg853Ter)

chr10:102402138-102402138

- Frequência: sem dados no gnomAD Exome e Genome
- VAF: 50% - HET
- OMIM: Immunodeficiency, common variable, 10 (AD)

ACMG: Patogênica (14 pontos)
- PS2: Relatado de novo no ClinVar nos seguintes envios: SCV000773792: ​​''Em pelo menos um indivíduo, a variante foi observada como de novo.' +1
- PM2: Frequência populacional 0,0%. Variante não encontrada nos genomas gnomAD +1
- PVS1: A perda de função é um mecanismo conhecido da doença: 15 variantes nulas patogênicas foram relatadas em ClinVar para este gene ( chr10:104161900:CA>C :HG19, chr10:104161908:CCAGCA>C :HG19, chr10:104161893:CCCGAGACA>C :HG19, chr10:104161895:C>T :HG19...), em 7 exons diferentes, dos quais 4 variantes neste exon (22) pontuação observada/esperada do gnomAD. +8
- PS3: Estudos experimentais mostraram que esse sinal de parada translacional prematuro afeta a função do NFKB2 (PMID: 24140114, 28778864) +4

A alteração NFKB2 p.Arg853Ter está altamente associada à imunodeficiência primária, explicando a maioria dos fenótipos listados, especialmente infecções recorrentes, hipogamaglobulinemia e sintomas gastrointestinais. 

## NT5C2 ✖️

HGVS: NT5C2(NM_001351173.2):c.141G>C
p.(Lys47Asn)

chr10:103139440-103139440

- Frequência: 0.017% no gnomAD Exome
- VAF: 41% - HET
- OMIM: Hereditary Spastic paraplegia 45 (AR)

O fenótipo associado não é condizente com o caso clínico.

## TCERG1L ✖️

HGVS: TCERG1L(NM_174937.4):c.721G>A
p.(Ala241Thr)

chr10:131260394-131260394

- Frequência: 0.0006% no gnomAD Exome
- VAF: HET

Não encontramos nenhuma patologia relacionada a essa variante nos banco de dados analisados (Varsome, Decipher, Clinvar).
