**Resultados chr1** 

## FastQC

R1:

![image](https://github.com/user-attachments/assets/8dcbe6da-39b6-496b-88ca-bed255ee69ca)

R2: 

![image](https://github.com/user-attachments/assets/304c4228-00c3-4107-884a-564a8bbb050d)

Interpretação: A análise de qualidade por base (Per Base Sequence Quality) apresentou valores de Phred variando entre 28 e 34 ao longo dos reads, indicando uma qualidade geral alta. Como valores de Phred acima de 30 representam uma taxa de erro inferior a 0,1%, os dados são considerados confiáveis para análises downstream


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

Após a anotação de variantes, com auxílio de ferramentas como o Pandas, para melhor visualização, filtramos o resultado final por "CLNSIG" (entrada no Clinvar) "Pathogenic", "Likely_pathogenic", "Uncertain_significance", para que todas essas variantes fossem chamadas para análise manual. 

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

## AK2 ✔️ achado complementar

Encontramos duas entradas nesse gene. 

1.
HGVS: AK2(NM_001625.4):c.614G>A
p.(Gly205Glu)

chr1:33013287-33013287

- Frequência: variante não passou nos critérios de qualidade do gnomAD Exome; 0.0563% no gnomad Genome.
- OMIM: disgenesia reticular (AR)
- VAF: 88% = variante em homozigose.

ACMG: Prov. Benigna (-3 pontos)
- BS1: Frequência alélica é maior do que o esperado para a doença. -4 (valor baseado no threshold do arquivo https://humgenomics.biomedcentral.com/articles/10.1186/s40246-023-00549-6/tables/2,  (< 0.01%))
- PP3: Evidências computacionais mostram efeito deletério (REVEL SCORE - Pathogenic Moderate) +1

2.
HGVS: AK2(NM_001625.4):c.602A>T
p.(Tyr201Phe)

chr1:33013299-33013299

- Frequência: 0.0028% no gnomad Exome
- OMIM: disgenesia reticular (AR)
- VAF: 88% = variante em homozigose.

ACMG: VUS (3 pontos)
PM2: Frequência alélica menor do que o esperado para a doença (< 0.01%) +2 > https://humgenomics.biomedcentral.com/articles/10.1186/s40246-023-00549-6/tables/2
PP3: Evidências computacionais mostram efeito deletério (REVEL SCORE - Pathogenic Moderate) +1

O Gene AK2 está relacionado com disgenesia reticular (AR), que é uma doença genética rara que afeta o sistema imunológico. É uma forma grave de imunodeficiência primária em crianças. 

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

## ABCA4 ✔️ achado primário

HGVS: ABCA4(NM_000350.3):c.6146del
p.(Lys2049ArgfsTer12)?

chr1: 94005442-94005442

- Frequência: <0.001% no gnomAD Exome
- VAF: 48% - variante em heterozigose. 
- LOF: Deleçao frameshift que gera um códon de parada prematuro.
- OMIM: degeneração macular relacionada à idade, tipo 2 - (AD)

ACMG: Patogênica (14 pontos)
- PVS1: Perda de função é um mecanismo conhecido da doença (o gene tem 1.029 variantes LOF patogênicas relatadas). O exon afeta 1 domínio funcional: domínio ABCA4_HUMAN da proteína UniProt 'transportador ABC 2'. O exon contém 39 variantes patogênicas. A região truncada contém 214 variantes patogênicas. +8
- PS4: Relatado no ClinVar em casos afetados nos seguintes envios: SCV001337922. ClinVar classifica esta variante como Patogênica, 2 estrelas (revisado em set. de 24, 3 envios), citando 5 artigos (31144483, 29975949, 25312043, 10958761 e 31814693). +4
- PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD. +2

A variante c.6146del resulta em uma frameshift (mudança do quadro de leitura), levando à formação de uma proteína truncada e potencialmente não funcional.
A degeneração macular relacionada à idade (DMRI) causa lesão progressiva da mácula, a área central e mais vital da retina, tendo como consequência a perda de visão.

## STRIP1 ✖️

HGVS: STRIP1(NM_033088.4):c.650+10G>A
p.? variante presente em região não-codificante

chr1: 110040713-110040713

- Frequência: <0.001% no gnomAD Exome
- VAF: 61% - variante em heterozigose.

Variante em região intrônica. Não há informações sobre ela no OMIM e as informações presentes no Clinvar (https://www.ncbi.nlm.nih.gov/clinvar/variation/930957/) não condizem com a clínica do paciente. 

## CASQ2 *atenção para achado secundário  ✖️
HGVS: CASQ2(NM_001232.4):c.926A>G
p.(Asp309Gly)

chr1: 115705205-115705205

- Frequência: <0.001% no gnomAD Exome
- VAF: 58% - variante em heterozigose.
- OMIM: Ventricular tachycardia, catecholaminergic polymorphic, 2 (AR)

ACMG: VUS (2 pontos)
- PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD. +1
- PP3: Variante missense com REVEL Score maior que 0.7 +1

O fenótipo associado não é condizente com o caso clínico. Além disso, a variante está em heterozigose e a condição é autossomica recessiva. 

## HSD3B2 ✖️

HGVS: HSD3B2(NM_000198.4):c.809T>C
p.(Ile270Thr)

chr1: 119422310-119422310

- Frequência: 0.073% no gnomAD Exome
- VAF: 50% - variante em heterozigose.
- OMIM: Adrenal hyperplasia, congenital, due to 3-beta-hydroxysteroid dehydrogenase 2 deficiency (AR)

ACMG: VUS
- BP4 -1

O fenótipo associado não é condizente com o caso clínico. Além disso, a variante está em heterozigose e a condição é autossomica recessiva. 

## FCGR3A ✖️

HGVS: FCGR3A(NM_001127592.2):c.509T>A
p.(Leu170His)

chr1: 161548543-161548543

- Frequência: 6% no gnomAD Exome

Essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara. 
