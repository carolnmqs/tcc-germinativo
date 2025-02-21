## Interpretação de variantes

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

![image](https://github.com/user-attachments/assets/cc10d674-98ba-4dde-a0a4-3b2fb60a3a72)



```python
patho = df['CLNSIG'].isin(["Pathogenic", "Likely_pathogenic", "Uncertain_significance", "Pathogenic/Likely_pathogenic"])
filtered_df = df[patho]
df[patho]

# filtered_df.to_csv('filtered_chr11.csv',index=False)
```
![image](https://github.com/user-attachments/assets/a2352455-f3cf-4d78-8ed6-9f7a5500753b)

Após a anotação de variantes, com auxílio de ferramentas como o Pandas, para melhor visualização, filtramos o resultado final por "CLNSIG" (entrada no Clinvar) "Pathogenic", "Likely_pathogenic", "Uncertain_significance" e "Pathogenic/Likely_pathogenic", para que todas essas variantes fossem chamadas para análise manual. 

Variantes encontradas: 4

Durante análise manual, observamos majoritariamente os campos "POS" (posição cromossômica), "CLNDN" (entrada no OMIM), "CLNSIG" (entrada no ClinVar), "GENEINFO" (nome do gene) e "ANN" para extrair informações como HGVS, transcrito, etc. e chegamos no seguinte resultado:

## MUC5AC ✖️

HGVS: MUC5AC(NM_001304359.2):c.10301C>T
p.(Pro3434Leu)

chr11:1188446-1188446

- Frequência: 15% gnomAD Exome

Essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara e não há evidências funcionais suficientes para classificá-la como patogênica.

ACMG: Benigna
- BA1 Stand Alone (frequência muito alta em controles populacionais saudáveis)

## TPP1 ✔️ complementar

HGVS: TPP1(NM_000391.4):c.509-1G>A
p.?

chr11:6617154-6617154

- Frequência: 0.003% gnomAD Exome
- VAF: 52% - variante em heterozigose.
- OMIM: Ceroid lipofuscinosis, neuronal e Spinocerebellar ataxia (ambos AR)

ACMG: Provavelmente Patogênica (12 pontos)
- PVS1: Perda de função é um mecanismo conhecido da doença (o gene tem 153 variantes LOF patogênicas relatadas). +8
- PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD. +2
- PM3: Casos afetados relatados nas seguintes submissões: IDs A interrupção deste local de splicing foi observada em indivíduos com lipofuscinose ceroide neuronal (PMID: 10330339, 10356316, 22832778).'' ClinVar: 207574 Relatado no ClinVar em casos afetados nas seguintes submissões: SCV000696665: ''Foi relatado em vários indivíduos afetados com LINCL com atividade enzimática residual quase ausente.'+2

`Doença relacionada a perda de fução do gene TPP1: A lipofuscinose ceróide neuronal (LCN) é um grupo de doenças neurodegenerativas raras, hereditárias e autossômicas recessivas. Também conhecida como doença de Batten, a LCN afeta o desenvolvimento cognitivo e motor. Principais sintomas Perda de visão, Convulsões, Declínio das capacidades mentais e motoras, Atraso na linguagem.

Paciente é heterozigoto para essa variante. Obs. pode ser que exista outra variante sem entrada no clinvar para esse gene que não foi encontrada em nossa busca.

## TECTA ✖️

HGVS: TECTA(NM_005422.4):c.4337C>G
p.(Thr1446Arg)

chr11:121157872-121157872

- Frequência: 0.002% gnomAD Exome
- VAF: 49% - variante em heterozigose.
- OMIM: Deafness, autosomal dominant 8/12

ACMG: VUS
- PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD
- BP4: Variante missense com REVEL score abaixo de 0.4

O fenótipo associado não é condizente com o caso clínico.

## MEN1 *atenção para achado secundário ✔️

HGVS: MEN1(NM_001370259.2):c.1548dup
p.(Lys517GlufsTer14)

chr11:64804619-64804619

- Frequência: <0.001% gnomAD Exome
- VAF: 52% - variante em heterozigose.
- OMIM: Multiple endocrine neoplasia 1 (AD)

ACMG: Patogênica (14 pontos)
- PVS1: Mudança de quadro, perda de função da Proteina +8
- PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD. +2

Como a variante está na lista de achados secundários da ACMG e foi classificada como patogênica, ela também será incluída no laudo.
