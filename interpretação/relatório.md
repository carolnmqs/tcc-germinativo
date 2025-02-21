## 1 Caso Clínico de Imunodeficiência e Distúrbios Endócrinos
PT-BR: Um paciente de 4 anos de idade apresentou histórico de infecções recorrentes do trato respiratório, otite média crônica e diarréia crônica. A avaliação imunológica revelou hipogamaglobulinemia, e foi iniciada terapia de reposição de imunoglobulina intravenosa. Aos 6 anos, o paciente desenvolveu alopecia e apresentou perda de peso, falta de apetite e pressão arterial baixa. Uma avaliação endócrina revelou deficiência do hormônio adrenocorticotrófico. Aos 9 anos, o paciente relatou dificuldades crescentes com a visão. O exame oftalmoscópico mostrou reflexo foveal bilateral reduzido, mas sem outras alterações fundoscópicas aparentes.

ENG-US: A 4-year-old patient presented with a history of recurrent respiratory tract infections, chronic otitis media, and chronic diarrhea. Immunological evaluation revealed hypogammaglobulinemia, and intravenous immunoglobulin replacement therapy was initiated. By the age of 6, the patient developed alopecia and presented with weight loss, lack of appetite, and low blood pressure. An endocrine evaluation revealed deficiency of adrenocorticotropic hormone. By the age of 9, the patient reported increasing difficulties with vision. Ophthalmoscope examination showed bilateral reduced foveal reflex but no other apparent fundoscopic changes

## 2 Descrição detalhada da metodologia utilizada


A metodologia utilizada neste estudo compreende uma série de etapas estruturadas para a análise de variantes germinativas a partir de dados de sequenciamento de nova geração (NGS). O fluxo de trabalho foi implementado na plataforma Google Colab, utilizando ferramentas amplamente reconhecidas na bioinformática.

## 2.1 Configuração do Ambiente
Inicialmente, foi criado um diretório no Google Drive para armazenar os dados e arquivos intermediários. Dentro deste diretório, foram criadas subpastas para organização dos arquivos de entrada e saída, incluindo:
- fastq (arquivos brutos de sequenciamento);
- bam (arquivos de alinhamento);
- vcf (arquivos de variantes);
- logs (arquivos de log de execução);
- reference (genoma de referência).

## 2.2 Instalação de Ferramentas
As ferramentas necessárias foram instaladas e registradas em arquivos de log:
- BWA para alinhamento de sequências;
- FastQC para controle de qualidade das leituras;
- Samtools para manipulação de arquivos BAM/SAM;
- Bedtools para manipulação de arquivos BED;
- BGZip e Tabix para compactação e indexação de arquivos VCF;
- GATK para chamadas de variantes;
- Picard para remoção de duplicatas;
- snpEff para anotação funcional das variantes;
- MultiQC para geração de relatórios.

## 2.3 Preparação do Genoma de Referência
O genoma de referência utilizado foi o hg38, obtido do UCSC Genome Browser. Foram realizados os seguintes passos:
- Download do cromossomo 1, 10 e 11 no formato FASTA;
- Indexação do genoma com BWA para permitir alinhamentos rápidos;
- Geração de um índice de referência com Samtools;
- Criação de um dicionário de sequência utilizando Picard.

## 2.4 Processamento de Dados de Sequenciamento
O arquivo FASTQ correspondente à amostra foi transferido para a estrutura de diretórios e processado da seguinte forma:
- Controle de Qualidade: O FastQC foi utilizado para avaliar a qualidade das leituras;
- Alinhamento: O alinhamento foi realizado utilizando BWA-MEM, adicionando informações do grupo de leitura (Read Group) para compatibilidade com as ferramentas subsequentes.
- Conversão e Indexação: O arquivo SAM gerado pelo BWA foi convertido para o formato BAM ordenado utilizando Samtools;
- O BAM foi indexado para permitir acesso rápido às leituras alinhadas.

## 2.5 Filtragem e Manipulação de Arquivos
Para garantir a qualidade dos dados alinhados, foram aplicadas etapas de filtragem:
- Conversão para BED e posterior fusão e ordenação dos intervalos com Bedtools;
- Remoção de duplicatas utilizando Picard MarkDuplicates;
- Nova indexação do BAM final.

## 2.6 Chamadas de Variantes
As variantes germinativas foram chamadas utilizando o GATK HaplotypeCaller, resultando em um arquivo VCF contendo as variantes identificadas. Este arquivo foi posteriormente comprimido e indexado para facilitar a manipulação.

## 2.7 Anotação Funcional das Variantes
A anotação das variantes foi realizada com snpEff, utilizando o banco de dados GRCh38.p14 para obtenção de informações funcionais. Além disso:
- Foram incorporadas anotações do ClinVar e dbSNP para identificar variantes de interesse clínico;
- O SnpSift foi utilizado para mapear as variantes contra o ClinVar e enriquecer as informações clínicas.

## 2.8 Extração de Informações e Geração de Relatórios
Os dados do arquivo VCF anotado foram convertidos para uma tabela utilizando GATK VariantsToTable, extraindo informações relevantes, como:
- Cromossomo e posição;
- Qualidade da chamada;
- Tipo da variante;
- Identificação no dbSNP e ClinVar;
- Informações sobre a cobertura e genotipagem;
- Frequência alélica na população.
- A análise exploratória dos dados foi conduzida utilizando Pandas, incluindo:
- Estatísticas descritivas das variantes;
- Contagem de diferentes tipos de variantes;
- Filtragem de variantes com significado clínico relevante.

## 2.9 Conclusão
O pipeline desenvolvido no Google Colab permitiu a análise eficiente de variantes germinativas presentes nos cromossomos 1, 10 e 11. As ferramentas utilizadas proporcionaram um fluxo de trabalho reprodutível e bem documentado, permitindo futuras análises com novos conjuntos de dados.

## 3 Resultados e Interpretação clínica das variantes

## 3.1 Resultado da análise do chr1
Para auxiliar na visualização e no processamento dos dados gerados após o processo de Anotação das Variantes, utilizamos a ferramenta Pandas (Imagem 1). 

<img width="1096" alt="Captura de Tela 2025-02-21 às 17 32 17" src="https://github.com/user-attachments/assets/1d6940a2-d741-4564-930d-90d5e3501d48" />

Imagem 1. Tabela de variantes genéticas identificadas no chr1 a partir da análise do Pipeline Germinativo, exibida em um dataframe do Pandas. As colunas incluem informações como cromossomo (CHROM), posição cromossômica (POS), qualidade das variantes (QUAL), tipo de variante (TYPE) e anotações clínicas do ClinVar (CLNDN, CLNSIG, CLNHGVS). Fonte: imagem de autoria própria.

Em seguida, filtramos os resultados obtidos pelo campo "CLNSIG" (classificação no ClinVar), selecionando apenas as variantes classificadas como "Pathogenic", "Likely_pathogenic", "Uncertain_significance" ou “Pathogenic/Likely_pathogenic" (Imagem 2). Dessa forma, observamos 11 variantes.

<img width="1191" alt="Captura de Tela 2025-02-21 às 17 33 05" src="https://github.com/user-attachments/assets/f90b7626-fc37-4069-ba57-b908bc9855d0" />

Imagem 2. Tabela de variantes genéticas anotadas no chr1, exibida em um dataframe do Pandas. Esta tabela inclui colunas identificadores como ID e ALLELEID, bem como informações clínicas do ClinVar: condição associada (CLNDN), classificação clínica (CLNSIG), tipo de variante (CLNVC), gene associado (GENEINFO) e nomenclatura HGVS (CLNHGVS). Fonte: imagem de autoria própria.

Por fim, foram encaminhadas para análise manual, focada, principalmente, nos campos: "POS" – posição cromossômica, "CLNDN" – condição associada (OMIM), "CLNSIG" – classificação no ClinVar, "GENEINFO" – nome do gene e "ANN" –  anotações funcionais e classificação das variantes com base nas diretrizes do American College of Medical Genetics and Genomics (ACMG), conforme a atualização de 2020, levando em consideração o modelo Bayesiano. 

<img width="560" alt="Captura de Tela 2025-02-21 às 17 44 01" src="https://github.com/user-attachments/assets/9541e754-54f6-431d-be6c-7adb9b48d8e7" />

Imagem 3. Valores pontuais para as categorias de força de evidência da ACMG/AMP. Disponível em: https://pmc.ncbi.nlm.nih.gov/articles/PMC8011844/

<img width="633" alt="Captura de Tela 2025-02-21 às 17 44 19" src="https://github.com/user-attachments/assets/49af7632-487f-4a2a-a7e8-4170447a19ce" />

Imagem 4. Categorias de classificação de variantes baseadas em pontos. Disponível em: https://pmc.ncbi.nlm.nih.gov/articles/PMC8011844/

A partir dessas informações, obteve-se: 

1) ACMG - Benigna.

- PERM1(NM_001394713.1):c.2330T>C p.(Val777Ala) - chr1:976215-976215.
Frequência em 79% no gnomAD Exome.
- TNFRSF1B(NM_001066.3):c.587T>G p.(Met196Arg) - chr1:12192898-12192898.
Frequência em 23% no gnomAD Exome.
- CDCA8(NM_001256875.2):c.799-8_799-6del p.? - chr1:37708312-37708312.
Frequência em 21% no gnomAD Exome.
- FCGR3A(NM_001127592.2):c.509T>A p.(Leu170His) - chr1: 161548543-161548543.
Frequência em 6% no gnomAD Exome.

BA1 Stand Alone: frequência alta em controles populacionais saudáveis, não sendo consideradas causadoras de uma doença monogênica rara.

2) ACMG - Patogênica. Achados primários.

ABCA4(NM_000350.3):c.6146del p.(Lys2049ArgfsTer12) - chr1: 94005442-94005442.
- Frequência em <0.001% no gnomAD Exome.
- VAF: 48% - variante em heterozigose.
- LOF: Deleçao frameshift que gera um códon de parada prematuro.
- OMIM: degeneração macular relacionada à idade, tipo 2 - (AD)

Assim, PVS1: Perda de função é um mecanismo conhecido da doença (o gene tem 1.029 variantes LOF patogênicas relatadas) (+8), PS4: Relatado em casos afetados no ClinVar (+4) e PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD (+2). 

Esta variante está associada a degeneração macular relacionada à idade (autossômica dominante) e outras doenças oculares autossômicas recessivas.

3) ACMG - Provavelmente Benigna.

HGVS: AK2(NM_001625.4):c.614G>A p.(Gly205Glu) - chr1:33013287-33013287.
- Frequência: não passou nos critérios de qualidade do gnomAD Exome.
- OMIM: Disgenesia reticular (AR)
- VAF: 88% = homozigose.

Assim, BS1: Frequência alélica é maior do que o esperado para a doença (-4) e PP3: Variante em região de splicing (+1).

4) ACMG - VUS.

HGVS: AK2(NM_001625.4):c.602A>T p.(Tyr201Phe) - chr1:33013299-33013299.
- Frequência: 0.0028% no gnomAD Exome.
- OMIM: disgenesia reticular (AR)
- VAF: 88% = homozigose.

Assim, BS1: Frequência alélica é maior do que o esperado para a doença (-4) e PP3: PP3: Evidências computacionais mostram efeito deletério (REVEL SCORE - Pathogenic Moderate) (+1)

O gene AK2 está relacionado à disgenesia reticular (AR), doença genética rara que afeta o sistema imunológico. Forma grave de imunodeficiência primária em crianças, correlacionada ao Caso Clínico.

HSD3B2(NM_000198.4):c.809T>C p.(Ile270Thr) - chr1: 119422310-119422310.
- Frequência: 0.073% no gnomAD Exome.
- OMIM: Adrenal hyperplasia, congenital, due to 3-beta-hydroxysteroid dehydrogenase 2 deficiency (AR).
- VAF: 50% = heterozigose.


5) Achados secundários.

O ACMG (American College of Medical Genetics and Genomics) possui uma lista de achados secundários que contém variantes genéticas consideradas clinicamente significantes, ou seja, podem ter implicações importantes para a saúde do paciente. Quando identificadas em testes genômicos, o ACMG recomenda que os laboratórios reportem as variantes nestes genes listados, uma vez que são associados a doenças genéticas hereditárias, principalmente condições cardiovasculares, oncológicas e metabólicas.

PCSK9(ENST00000302118.5):c.996+44A>G p.? - chr1:55056233-55056233.
- Frequência em 49% no gnomAD Exome.
No entanto, essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara.

CASQ2(NM_001232.4):c.926A>G; p.(Asp309Gly) - chr1: 115705205-115705205.
- Frequência em <0.001% no gnomAD Exome.
- OMIM: Ventricular tachycardia, catecholaminergic polymorphic, 2 (AR).
- VAF: 58% = heterozigose.

Assim, VUS, PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD (+2) e PP3: Variante missense com REVEL Score maior que 0.7 (+1).

O fenótipo associado não é condizente com o caso clínico e a variante está em heterozigose, mas a condição é autossômica recessiva.


6) Outros.

STRIP1(NM_033088.4):c.650+10G>A p.? - chr1: 110040713-110040713.
- Frequência em 49% no gnomAD Exome.
- VAF: 61% = heterozigose.
- 
Variante em região intrônica. Não há informações sobre ela no OMIM e as informações presentes no Clinvar não condizem com a clínica do paciente.

## 3.2 Resultado da análise do chr10

Assim como realizado no chr1, utilizamos a ferramenta Pandas para auxiliar na visualização e no processamento dos dados gerados após o processo de Anotação das Variantes.

Posteriormente, filtramos os resultados obtidos pelo campo "CLNSIG" (classificação no ClinVar), selecionando apenas as variantes classificadas como "Pathogenic", "Likely_pathogenic", "Uncertain_significance" ou “Pathogenic/Likely_pathogenic" (Imagem 4). Dessa forma, observamos 7 variantes.

Como próximo passo, foram encaminhadas para análise manual, focada, principalmente, nos campos: "POS" – posição cromossômica, "CLNDN" – condição associada (OMIM), "CLNSIG" – classificação no ClinVar, "GENEINFO" – nome do gene e "ANN" – anotações funcionais. A partir dessas informações, observou-se:

1) ACMG - Benigna.

CHST3(NM_004273.5):c.*66G>T p.? - chr10:72008537-72008537.
- Frequência: 0.14% gnomAD Exome.
- VAF: 80,71% = heterozigose.

Assim, BP7: variante sinônima ou não codificadora que não está localizada em uma região de splicing e não se prevê que tenha consequências que alterem o splicing (-1).


2) ACMG - Patogênica. Achados primários.


NFKB2(NM_001322934.2):c.2557C>T p.(Arg853Ter) - chr10:102402138-102402138.
- Frequência: sem dados no gnomAD Exome.
- VAF: 51,17% = heterozigose.
- OMIM: Immunodeficiency, common variable, 10 (AD).

Assim, PS2: Relatado no ClinVar , em pelo menos um indivíduo, observada como de novo (+4), PM2: Frequência populacional 0,0%, ariante não encontrada no gnomAD (+2), PVS1: A perda de função é um mecanismo conhecido da doença: 15 variantes nulas patogênicas foram relatadas no ClinVar para este gene, em 7 exons diferentes (+8) e PS3: Estudos experimentais mostraram que esse sinal de parada translacional prematuro afeta a função do NFKB2 (+4).
Dessa forma, esta alteração no gene NFKB2 está altamente associada à imunodeficiência primária, explicando a maioria dos fenótipos listados, especialmente infecções recorrentes, hipogamaglobulinemia e sintomas gastrointestinais.

3) ACMG - VUS.


TMEM254(NM_001270367.1):c.378C>A p.(Phe126Leu) - chr10:80090851-80090851.
- Frequência: 0.019% gnomAD Exome.
- VAF: 41,21% = heterozigose.

Assim, PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD (+2) e BP4: Múltiplas linhas de evidências computacionais sugerem nenhum impacto no gene ou no produto genético (conservação, evolução, impacto de splicing, etc.) (-1).
Não encontramos nenhuma patologia relacionada a esse gene.

4) Outros.


TBATA(NM_001318241.2):c.666T>C p.(Ala222=) - chr10:70777180-70777180.
- Frequência: 29% gnomAD Exome.
- VAF: 54,30% = heterozigose.

NRG3(NM_001370081.1):c.1986C>T p.(Ser662=) - chr10:82985500-82985500.
- Frequência: 32% gnomAD Exome.
- VAF: 44,31% = heterozigose.
Essas variantes têm frequência alta na população para serem causadoras de uma doença monogênica rara.

NT5C2(NM_001351173.2):c.141G>C p.(Lys47Asn) - chr10:103139440-103139440.
- Frequência: 0.017% gnomAD Exome.
- VAF: 41% = heterozigose.
- OMIM: Hereditary Spastic paraplegia 45 (AR).

O fenótipo associado não é condizente com o caso clínico.

TCERG1L(NM_174937.4):c.721G>A p.(Ala241Thr) - chr10:131260394-131260394.
- Frequência: 0.0006% gnomAD Exome.
- VAF: 40% = heterozigose. 
Não encontramos nenhuma patologia relacionada a essa variante.


## 3.3 Resultado da análise do chr11
Seguindo o processo realizado com o chr1 chr10, utilizamos, também, a ferramenta Pandas para auxiliar na visualização e no processamento dos dados gerados após o processo de Anotação das Variantes no chr11 (Imagem 5). 

Então, filtramos os resultados obtidos pelo campo "CLNSIG" (classificação no ClinVar), selecionando apenas as variantes classificadas como "Pathogenic", "Likely_pathogenic", "Uncertain_significance" ou “Pathogenic/Likely_pathogenic" (Imagem 6). Dessa forma, observamos 4 variantes.

Em seguida, foram encaminhadas para análise manual, focada, principalmente, nos campos: "POS" – posição cromossômica, "CLNDN" – condição associada (OMIM), "CLNSIG" – classificação no ClinVar, "GENEINFO" – nome do gene e "ANN" – anotações funcionais. A partir dessas informações, observou-se:


1) ACMG - Benigna

MUC5AC(NM_001304359.2):c.10301C>T p.(Pro3434Leu) - chr11:1188446-1188446.
- Frequência: 15% gnomAD Exome.

Essa variante tem frequência muito alta na população para ser causadora de uma doença monogênica rara e não há evidências funcionais suficientes para classificá-la como patogênica.
Assim, BA1 Stand Alone.

2) ACMG - VUS.

TECTA(NM_005422.4):c.4337C>G p.(Thr1446Arg) - chr11:121157872-121157872.
- Frequência: 0.002% gnomAD Exome.
- VAF: 49% = heterozigose. OMIM: Deafness (AD).
Assim, PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD (+2) e BP4: Variante missense com REVEL score abaixo de 0.4. (-1)
O fenótipo associado não é condizente com o caso clínico.


3) ACMG - Provavelmente Patogênica. Achado secundário.

TPP1(NM_000391.4):c.509-1G>A p.? - chr11:6617154-6617154.
- Frequência: 0.003% gnomAD Exome.
- VAF: 52% = heterozigose.
- OMIM: Ceroid lipofuscinosis, neuronal e Spinocerebellar ataxia (AR).

Assim, PVS1: Perda de função é um mecanismo conhecido da doença (o gene tem 153 variantes LOF patogênicas relatadas (+8), PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD (+2) e PM3: Casos afetados relatados no ClinVar. A interrupção deste local de splicing foi observada em indivíduos com lipofuscinose ceroide neuronal (+2).
Doença relacionada a perda de função do gene TPP1. A lipofuscinose ceróide neuronal (LCN) é um grupo de doenças neurodegenerativas raras, hereditárias e autossômicas recessivas. Também conhecida como doença de Batten, a LCN afeta o desenvolvimento cognitivo e motor. Principais sintomas são: perda de visão, convulsões, declínio das capacidades mentais e motoras, atraso na linguagem.


4) ACMG - Patogênica. Achado secundário.


MEN1(NM_001370259.2):c.1548dup p.(Lys517GlufsTer14) - chr11:64804619-64804619.
- Frequência: <0.001% gnomAD Exome.
- VAF: 52% = heterozigose.
- OMIM: Multiple endocrine neoplasia 1 (AD).
Assim, PVS1: Mudança de quadro, perda de função da Proteina (+8), PM2: Frequência extremamente baixa em bancos de dados populacionais do gnomAD (+2) e PS4: Para doenças raras dominantes, que apareceram em casos afetados enquanto eram extremamente raras na população (+4).

## 5 Conclusões e recomendações
Desenvolvemos um pipeline completo para a análise de variantes germinativas, permitindo a identificação e classificação de variantes com potencial impacto clínico. O fluxo de trabalho estruturado possibilitou a detecção de 22 variantes, distribuídas da seguinte forma: 
- 11 variantes no chr 1;
- 7 variantes no chr 10;
- 4 variantes no chr 11.

Entre essas variantes, destacamos ABCA4 e NFKB2 como candidatas relevantes, associadas a distúrbios oculares e imunodeficiência, respectivamente.

Além disso, identificamos três achados secundários de relevância clínica: AK2, TPP1 e MEN1, que podem ser incluídos no laudo como informações complementares, auxiliando em possíveis direcionamentos clínicos e aconselhamento genético.

A análise realizada permitiu identificar e classificar variantes germinativas, com o auxílio de ferramentas e bases de dados, como, por exemplo, o ClinVar, Varsome e OMIM, para investigar predisposições genéticas e sua relação com o caso clínico. Além disso, ao desenvolver um script para o pipeline germinativo, torna-se possível reutilizá-lo em futuras análises, garantindo maior eficiência e reprodutibilidade.

## 6 Referências

Broad Institute. (2024). Genome Analysis Toolkit (GATK). Disponível em: https://gatk.broadinstitute.org/
DECIPHER. (2024). Database of Chromosomal Imbalance and Phenotype in Humans Using Ensembl Resources. Disponível em: https://www.deciphergenomics.org/
Ensembl Genome Browser. (2024). Variant Effect Predictor (VEP). Disponível em: https://www.ensembl.org/info/docs/tools/vep/
Feng, B. J., Barnard, A. M., & Foulkes, W. D. (2023). Validation and extension of a point-based scoring system for classifying germline variants using the ACMG/AMP guidelines. Human Genomics, 17(1), 18. https://doi.org/10.1186/s40246-023-00549-6
Franklin by Genoox. (2024). Genomic Variant Interpretation Platform. Disponível em: https://franklin.genoox.com/
National Center for Biotechnology Information (NCBI). (2024). ClinVar: Public archive of interpretations of clinically relevant variants. Disponível em: https://www.ncbi.nlm.nih.gov/clinvar/
Online Mendelian Inheritance in Man (OMIM). (2024). An Online Catalog of Human Genes and Genetic Disorders. Disponível em: https://www.omim.org/
Tavtigian, S. V., Greenblatt, M. S., Harrison, S. M., Nussbaum, R. L., Prabhu, S. A., Boucher, K. M., et al. (2020). Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. Genetics in Medicine, 22(1), 77–86. https://doi.org/10.1038/s41436-019-0660-1
UCSC Genome Browser. (2024). Genome Bioinformatics. Disponível em: https://genome.ucsc.edu/
VarSome. (2024). The Human Genomics Community. Disponível em: https://varsome.com/


