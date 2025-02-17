**Resultados chr1** 

## FastQC

R1:

![image](https://github.com/user-attachments/assets/8dcbe6da-39b6-496b-88ca-bed255ee69ca)

Interpretação:

R2: 

![image](https://github.com/user-attachments/assets/304c4228-00c3-4107-884a-564a8bbb050d)

Interpretação: 


## Interpretação de Variantes

```python
MeuDrive = "/content/drive/MyDrive/AtividadeFinalGerminativo_chr1"
amostra = "cap-ngse-b-2019-chr1"

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

caminho_arquivo = f"{MeuDrive}/{amostra}/output/{amostra}.clinvar.ann.txt"
df = pd.read_csv(caminho_arquivo, sep="\t")
df
```

![image](https://github.com/user-attachments/assets/8ab4e6e2-41cc-4554-8725-80dbe5c7dcc8)

```python
patho = df['CLNSIG'].isin(["Pathogenic", "Likely_pathogenic", "Uncertain_significance"])
df[patho]
```

![image](https://github.com/user-attachments/assets/4aa4e116-d070-4d64-bde2-2717c17b33d3)

Após a anotação de variantes, com auxílio de ferramentas como o Pandas, para melhor visualização, filtramos o resultado final por "CLNSIG" (entrada no Clinvar) "Pathogenic", "Likely_pathogenic", "Uncertain_significance", para que todas as variantes com entrada e sem entrada fossem chamadas para análise manual. 

Variantes encontradas: 11

Durante análise manual, observamos majoritariamente os campos "POS" (posição cromossômica), "CLNDN" (entrada no OMIM), "CLNSIG" (entrada no ClinVar), "GENEINFO" (nome do gene) e "ANN" para extrair informações como HGVS, transcrito, etc. e chegamos no seguinte resultado:

## PERM1 ✖️

HGVS: PERM1(NM_001394713.1):c.2330T>C
p.(Val777Ala)

chr1:976215-976215

- Frequência: 79% gnomAD Exome

Essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara e não há evidências funcionais suficientes para classificá-la como patogênica.

ACMG: Benigna
- BA1 Stand Alone (frequência muito alta em controles populacionais saudáveis)

## TNFRSF1B ✖️

HGVS: TNFRSF1B(NM_001066.3):c.587T>G
p.(Met196Arg)

chr1:12192898-12192898

- Frequência: 23% gnomAD Exome

Essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara. 

ACMG: Benigna
- BA1 Stand Alone (frequência muito alta em controles populacionais saudáveis)

## AK2

Encontramos duas entradas nesse gene. 

1.
HGVS: AK2(NM_001625.4):c.614G>A
p.(Gly205Glu)

chr1:33013287-33013287

- Frequência: variante não passou nos critérios de qualidade do gnomAD Exome; 0.0563% no gnomad Genome.
- OMIM: sem dados catalogados para análise. 
- VAF: 88% = variante em homozigose. 

2.
HGVS: AK2(NM_001625.4):c.602A>T
p.(Tyr201Phe)

chr1:33013299-33013299

- Frequência: 0.0028% no gnomad Exome
- OMIM: sem dados catalogados para análise. 
- VAF: 88% = variante em homozigose. 

## CDCA8 ✖️

HGVS:CDCA8(NM_001256875.2):c.799-8_799-6del
p.?

chr1:37708312-37708312

- Frequência: 21% gnomAD Exome

Essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara. 

ACMG: Benigna
- BA1 Stand Alone (frequência muito alta em controles populacionais saudáveis)

## PCSK9 *atenção para achado secundário ✖️

HGVS: PCSK9(ENST00000302118.5):c.996+44A>G
p.? - variante presente em região não-codificante

chr1:55056233-55056233

- Frequência: 49% no gnomad Exome

Essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara. 

## ABCA4 ✔️

HGVS: ABCA4(NM_000350.3):c.6146del
p.(Lys2049ArgfsTer12)?

chr1: 94005442-94005442

- Frequência: <0.001% no gnomAD Exome
- VAF: 48% - variante em heterozigose. 
- LOF: Deleçao frameshift que gera um códon de parada prematuro.

```A variante c.6146del resulta em uma frameshift (mudança do quadro de leitura), levando à formação de uma proteína truncada e potencialmente não funcional. Esta variante está associada a doenças oculares autossômicas recessivas, como a Doença de Stargardt. O paciente está em heterozigose para esta variante, o que sugere que ele é portador assintomático. Não se espera que o paciente desenvolva a condição associada à variante. Contudo, é importante observar que o paciente pode transmitir a variante para futuros filhos, e se ambos os pais forem portadores, existe 25% de chance de um filho herdar as duas cópias da mutação e apresentar a condição. Recomenda-se aconselhamento genético para discussão das opções de testes para familiares.```

## STRIP1

HGVS: STRIP1(NM_033088.4):c.650+10G>A
p.? variante presente em região não-codificante

chr1: 110040713-110040713

- Frequência: <0.001% no gnomAD Exome

## CASQ2 *atenção para achado secundári

HGVS: CASQ2(NM_001232.4):c.926A>G
p.(Asp309Gly)

chr1: 115705205-115705205

- Frequência: <0.001% no gnomAD Exome

## HSD3B2

HGVS: HSD3B2(NM_000198.4):c.809T>C
p.(Ile270Thr)

chr1: 119422310-119422310

- Frequência: <0.001% no gnomAD Exome

## FCGR3A ✖️

HGVS: FCGR3A(NM_001127592.2):c.509T>A
p.(Leu170His)

chr1: 161548543-161548543

- Frequência: 6% no gnomAD Exome

Essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara. 
