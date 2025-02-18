## Atividade Final - chr11

Link pro Colab para execução dos códigos: 

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1KR1l_J46_bDVmNP9RRh-jTfYNmpwVCGt#scrollTo=X4DKx1kv84pc)

## Configuração do ambiente 

Importação do módulo do Google Colab + montagem do Drive. O código permite acessar arquivos do Drive diretamente no ambiente do Colab, permitindo que arquivos gerados durante o processo sejam armazenados. 

```bash
from google.colab import drive
drive.mount('/content/drive', force_remount=True)
```

Criação do diretório principal dentro do Google Drive. 

**Observação:** Sempre utilize `%%bash` no Colab para rodar comandos de shell, pois o ambiente é baseado em Python e o `%%bash` permite a execução de comandos do sistema diretamente.

```bash
%%bash
mkdir /content/drive/MyDrive/AtividadeFinalGerminativo_chr11
```

Criação da estrutura de diretórios no Google Drive, dentro do caminho `/content/drive/MyDrive/AtividadeFinalGerminativo_chr11`, com subpastas para armazenar dados de fastq, bam, vcf, logs e referências.

**Observação:** A partir deste ponto, será necessário especificar a variável `MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"` em cada novo bloco de código, devido a uma peculiaridade do Google Colab que não mantém variáveis entre blocos executados.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"

mkdir $MeuDrive/dados
mkdir $MeuDrive/dados/fastq
mkdir $MeuDrive/dados/bam
mkdir $MeuDrive/dados/vcf
mkdir $MeuDrive/logs
mkdir $MeuDrive/reference
mkdir $MeuDrive/reference/hg38
```

Instalação dos programas necessários para rodar o pipeline (foi quebrado em algumas partes para garantir que todas as etapas estavam rodando normalmente no Colab, além disso, o GATK demora um pouco para rodar, em caso de erros ou desconexão, não precisaríamos rodar tudo novamente). 

`sudo apt install` instala o software e redireciona tanto a saída padrão `(stdout)` quanto os erros `(stderr)` para o arquivo de log em `$MeuDrive/logs/`. 

```bash
%%bash

MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"

echo '1 - Instalação de programas'
mkdir -p logs

echo 'Instalando o bwa'
sudo apt install bwa 1>$MeuDrive/logs/bwa.log 2>$MeuDrive/logs/bwa.log

echo 'Instalando o fastqc'
sudo apt install fastqc 1>$MeuDrive/logs/fastqc.log 2>$MeuDrive/logs/fastqc.log

echo 'Instalando o samtools'
sudo apt install samtools 1>$MeuDrive/logs/samtools.log 2>$MeuDrive/logs/samtools.log

echo 'Instalando o bedtools'
sudo apt install bedtools 1>$MeuDrive/logs/bedtools.log 2>$MeuDrive/logs/bedtools.log

echo 'Instalando o bgzip'
sudo apt install bgzip 1>$MeuDrive/logs/bgzip.log 2>$MeuDrive/logs/bgzip.log

echo 'Instalando o tabix'
sudo apt install tabix 1>$MeuDrive/logs/tabix.log 2>$MeuDrive/logs/tabix.log
```
`wget` baixa arquivos da internet através do link fornecido.

```bash
%%bash

MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"

echo 'Instalando o gatk'
wget https://github.com/broadinstitute/gatk/releases/download/4.1.8.1/gatk-4.1.8.1.zip 1>$MeuDrive/logs/gatk.log 2>$MeuDrive/logs/gatk.log
unzip gatk-4.1.8.1.zip 1>$MeuDrive/logs/gatk.log 2>$MeuDrive/logs/gatk.log
rm gatk-4.1.8.1.zip
```
```bash
%%bash
MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"

echo 'Instalando o picard'
wget https://github.com/broadinstitute/picard/releases/download/2.24.2/picard.jar 1>$MeuDrive/logs/picard.log 2>$MeuDrive/logs/picard.log

echo 'Instalando o snpEff'
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip 1>$MeuDrive/logs/snpEff.log 2>$MeuDrive/logs/snpEff.log
unzip snpEff_latest_core.zip 1>$MeuDrive/logs/snpEff.log 2>$MeuDrive/logs/snpEff.log
rm snpEff_latest_core.zip

echo 'Instalando o multiqc'
sudo apt install multiqc 1>$MeuDrive/logs/multiqc.log 2>$MeuDrive/logs/multiqc.log
```

`curl` realiza requisições HTTP, permitindo baixar arquivos através do link.

```bash
%%bash

MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"

echo '2 - Preparação do Genoma de Referência'

echo 'Baixando o Genoma de Referência'

mkdir -p reference

# baixando o chr11
curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr11.fa.gz" | \
  gunzip -c > $MeuDrive/reference/hg38.fasta

echo 'Indexando o Genoma de Referência'
bwa index \
  -a bwtsw \
  /content/reference/hg38.fasta 1>$MeuDrive/logs/bwa.log 2>$MeuDrive/logs/bwa.log

samtools faidx $MeuDrive/reference/hg38.fasta

java -jar picard.jar CreateSequenceDictionary \
    REFERENCE=$MeuDrive/reference/hg38.fasta \
    OUTPUT=$MeuDrive//reference/hg38.dict 1>$MeuDrive/logs/picard.log 2>$MeuDrive/logs/picard.log
```

## Script para Análise Germinativa

Antes de inciar, é preciso garantir que os arquivos FASTQ R1 e R2 estejam na pasta `dados/fastq/` do Drive e sejam devidamente nomeados nas variáveis "amostra" e "fastq1" e "fastq2".

```bash
%%bash

MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"

amostra="cap-ngse-b-2019-chr11"

fastq1="$MeuDrive/dados/fastq/cap-ngse-b-2019-chr11_S1_L001_R1_001.fastq.gz"
fastq2="$MeuDrive/dados/fastq/cap-ngse-b-2019-chr11_S1_L001_R2_001.fastq.gz"

echo "Iniciando o processamento da $amostra"

# Criando estrutura de diretórios
mkdir -p $MeuDrive/$amostra/input
mkdir -p $MeuDrive/logs

# Movendo arquivos se ainda não estiverem no destino
if [ -e "$fastq1" ]; then
  echo "Movendo '$fastq1' para o diretório input..."
  mv "$fastq1" "$MeuDrive/$amostra/input/"
else
  echo "O arquivo '$fastq1' não existe."
fi

if [ -e "$fastq2" ]; then
  echo "Movendo '$fastq2' para o diretório input..."
  mv "$fastq2" "$MeuDrive/$amostra/input/"
else
  echo "O arquivo '$fastq2' não existe."
fi

echo "==========================="
echo "Início do Pipeline de Análise de Variantes Germinativas"
echo "==========================="

# Preprocessamento de dados
echo "Preprocessamento de dados (ex. QC, trimming)"

input_fastq1="$MeuDrive/$amostra/input/$(basename $fastq1)"
input_fastq2="$MeuDrive/$amostra/input/$(basename $fastq2)"

fastqc "$input_fastq1" >> "$MeuDrive/logs/fastqc.log" 2>&1
fastqc "$input_fastq2" >> "$MeuDrive/logs/fastqc.log" 2>&1
```

A partir deste ponto,  será necessário especificar as variáveis `MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"` e `amostra="cap-ngse-b-2019-chr11"` em cada novo bloco de código, devido a uma peculiaridade do Google Colab que não mantém variáveis entre blocos executados.

```bash
%%bash

MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"
amostra="cap-ngse-b-2019-chr11"
fastq1="cap-ngse-b-2019-chr11_S1_L001_R1_001.fastq.gz"
fastq2="cap-ngse-b-2019-chr11_S1_L001_R2_001.fastq.gz"
input_fastq1="$MeuDrive/$amostra/input/$fastq1"
input_fastq2="$MeuDrive/$amostra/input/$fastq2"
reference="$MeuDrive/reference/hg38.fasta"

mkdir -p "$MeuDrive/$amostra/output"
mkdir -p "$MeuDrive/logs"

# Verificação se o genoma está indexado
if [ ! -f "${reference}.bwt" ]; then
  echo "Indexando o genoma de referência..."
  bwa index "$reference"
fi

# Validação dos arquivos de entrada
if [ ! -f "$input_fastq1" ] || [ ! -f "$input_fastq2" ]; then
  echo "Erro: Arquivo de entrada FASTQ não encontrado!"
  exit 1
fi

# Alinhamento com BWA
echo "Alinhamento de Sequências (ex. BWA)"
bwa mem -R "@RG\tID:$amostra\tSM:$amostra\tLB:$amostra\tPL:ILLUMINA" \
    "$reference" \
    "$input_fastq1" \
    "$input_fastq2" > "$MeuDrive/$amostra/output/$amostra.sam" 2> "$MeuDrive/logs/bwa.log"

if [ $? -eq 0 ]; then
  echo "Alinhamento concluído com sucesso."
else
  echo "Erro durante o alinhamento. Verifique os logs."
fi
```
```bash
%%bash

MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"
amostra="cap-ngse-b-2019-chr11"


# Conversão SAM para BAM e indexação
echo "Conversão para BAM e indexação"

samtools sort -O bam -o $MeuDrive/$amostra/output/$amostra.bam $MeuDrive/$amostra/output/$amostra.sam
samtools index $MeuDrive/$amostra/output/$amostra.bam
```
```bash
%%bash
MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"
amostra="cap-ngse-b-2019-chr11"

bedtools bamtobed -i $MeuDrive/$amostra/output/$amostra.bam >$MeuDrive/$amostra/output/$amostra.bed
bedtools merge -i $MeuDrive/$amostra/output/$amostra.bed >$MeuDrive/$amostra/output/$amostra.merged.bed
bedtools sort -i $MeuDrive/$amostra/output/$amostra.merged.bed >$MeuDrive/$amostra/output/$amostra.sorted.bed
```

Etapa extra: Mark Duplicates (identifica e marca duplicatas no BAM).

```bash
%%bash

MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"
amostra="cap-ngse-b-2019-chr11"

java -jar picard.jar MarkDuplicates \
    I=$MeuDrive/$amostra/output/$amostra.bam \
    O=$MeuDrive/$amostra/output/$amostra.marked.bam \
    M=$MeuDrive/$amostra/output/$amostra.metrics.txt \
    REMOVE_DUPLICATES=true
```

A conversão de `$amostra.bam` para `$amostra.marked.bam` ocorreu devido à etapa de marcação de duplicatas com Picard. Para que o pipeline prossiga corretamente, foi necessário indexar `$amostra.marked.bam` usando `samtools index`, gerando o arquivo `.bai`, que permite acesso rápido às regiões do BAM durante análises posteriores.

```bash
%%bash

MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"
amostra="cap-ngse-b-2019-chr11"

samtools index $MeuDrive/$amostra/output/$amostra.marked.bam

# Chamadas de variantes germinativo
echo "Chamadas de variantes (GATK HaplotypeCaller)"

gatk-4.1.8.1/gatk HaplotypeCaller --verbosity ERROR \
    -R $MeuDrive/reference/hg38.fasta \
    -I $MeuDrive/$amostra/output/$amostra.marked.bam \
    -O $MeuDrive/$amostra/output/$amostra.vcf

bgzip -f $MeuDrive/$amostra/output/$amostra.vcf
tabix -f $MeuDrive/$amostra/output/$amostra.vcf.gz
```
```bash
%%bash

MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"
amostra="cap-ngse-b-2019-chr11"

# Geração de relatórios ou visualizações
echo "Geração do relatório final"

multiqc $MeuDrive/$amostra/

# Fim do pipeline
echo "==========================="
echo "Pipeline de Bioinformática Concluído!"
echo "==========================="
```

## Anotação de Variantes

O SnpEff é uma ferramenta de anotação e previsão de impacto funcional de variantes. Durante os próximos passos, precisaremos baixar a ferramenta e anotar as variantes germinativas do chr10 usando `SnpEff` e `SnpSift`, adicionando informações do ClinVar e dbSNP a um VCF processado.

```bash
%%bash

wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip -o snpEff_latest_core.zip
rm snpEff_latest_core.zip
```

Para rodar sem erros no Colab, foi necessário realizar o update da versão do java. 

```bash
%%bash

sudo apt update
sudo apt install openjdk-21-jre

java -jar snpEff/snpEff.jar download -v GRCh38.p14
```
```bash
%%bash
MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"
amostra="cap-ngse-b-2019-chr11"

java -Xmx8g -jar snpEff/snpEff.jar -v GRCh38.p14 \
    -stats $MeuDrive/$amostra/output/$amostra.html \
    $MeuDrive/$amostra/output/$amostra.vcf.gz > $MeuDrive/$amostra/output/$amostra.ann.vcf
```

```bash
%%bash

mkdir snpEff/./db
mkdir snpEff/./db/GRCh38
mkdir snpEff/./db/GRCh38/clinvar
mkdir snpEff/./db/GRCh38/dbSnp


wget -O snpEff/./db/GRCh38/clinvar/clinvar-latest.vcf.gz \
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -O snpEff/./db/GRCh38/clinvar/clinvar-latest.vcf.gz.tbi \
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

wget -O snpEff/./db/GRCh38/dbSnp/dbSnp.vcf.gz \
    ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz
wget -O snpEff/./db/GRCh38/dbSnp/dbSnp.vcf.gz.tbi \
    ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz.tbi
```

```bash
%%bash
MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"
amostra="cap-ngse-b-2019-chr11"

java -Xmx1g -jar snpEff/SnpSift.jar \
    annotate \
    snpEff/./db/GRCh38/clinvar/clinvar-latest.vcf.gz \
    $MeuDrive/$amostra/output/$amostra.ann.vcf \
    > $MeuDrive/$amostra/output/$amostra.clinvar.ann.vcf
```

Estamos utilizando um arquivo VCF anotado (`$amostra.clinvar.ann.vcf`) e usando o `GATK VariantsToTable` para extrair várias informações importantes das variantes, como posição, hgvs, frequência alélica, entrada no Clinvar (...). O comando vai gerar uma tabela no formato TXT, que será utilizada para análise posterior.

```bash
%%bash

MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr11"
amostra="cap-ngse-b-2019-chr11"

echo "Extraindo informações do VCF"

gatk-4.1.8.1/gatk VariantsToTable \
    -V "$MeuDrive/$amostra/output/$amostra.clinvar.ann.vcf" \
    -F CHROM -F POS -F QUAL -F TYPE -F ID -F ALLELEID \
    -F CLNDN -F CLNSIG -F CLNSIGCONF -F CLNSIGINCL -F CLNVC \
    -F GENEINFO -F AF_EXAC -F CLNHGVS -F ANN -F CLNDNINCL \
    -F ONCDISDBINCL -F LOF \
    -GF AD -GF DP -GF GQ -GF GT -GF AF \
    -O "$MeuDrive/$amostra/output/$amostra.clinvar.ann.txt"
```
