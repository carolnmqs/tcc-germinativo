## Atividade Final - chr1

Link pro Colab para execução dos códigos: 

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1BPivAbGdAe56dTSOIwi09F7N1WoY5aAm#scrollTo=NdujBZr4yyaL)

## Configuração do ambiente 

Importação do módulo do Google Colab + montagem do Drive. O código permite acessar arquivos do Drive diretamente no ambiente do Colab, permitindo que arquivos gerados durante o processo sejam armazenados. 

```bash
from google.colab import drive
drive.mount('/content/drive', force_remount=True)
```

Criação do diretório principal dentro do Google Drive. Obs. sempre usar %%bash no colab, já que o ambiente é baseado em python. 

```bash
%%bash
mkdir /content/drive/MyDrive/AtividadeFinalGerminativo_chr1
```

Criação da estrutura de diretórios no Google Drive, dentro do caminho `/content/drive/MyDrive/AtividadeFinalGerminativo_chr1`, com subpastas para armazenar dados de fastq, bam, vcf, logs e referências.

Observação: A partir deste ponto, será necessário especificar a variável `MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr1"` em cada novo bloco de código, devido a uma peculiaridade do Google Colab que não mantém variáveis entre blocos executados.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/AtividadeFinalGerminativo_chr1"

mkdir $MeuDrive/dados
mkdir $MeuDrive/dados/fastq
mkdir $MeuDrive/dados/bam
mkdir $MeuDrive/dados/vcf
mkdir $MeuDrive/logs
mkdir $MeuDrive/reference
mkdir $MeuDrive/reference/hg38
```

