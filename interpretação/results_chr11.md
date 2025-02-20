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

## TPP1

HGVS: TPP1(NM_000391.4):c.509-1G>A
p.?

chr11:6617154-6617154

- Frequência: 0.003% gnomAD Exome
- VAF: 52% - variante em heterozigose.
- OMIM: 

## TECTA ✖️

HGVS: TECTA(NM_005422.4):c.4337C>G
p.(Thr1446Arg)

chr11:121157872-121157872

- Frequência: 0.002% gnomAD Exome
- VAF: 49% - variante em heterozigose.
- OMIM: Deafness, autosomal dominant 8/12

O fenótipo associado não é condizente com o caso clínico.

## MEN1 *atenção para achado secundário

HGVS: MEN1(NM_001370259.2):c.1548dup
p.(Lys517GlufsTer14)

chr11:64804619-64804619

- Frequência: <0.001% gnomAD Exome
- VAF: 52% - variante em heterozigose.
- OMIM: Multiple endocrine neoplasia 1 (AD)

