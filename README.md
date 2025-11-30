## Bioinformatics Learn

## Создание окружения

###  Miniforge3:
```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
~/mambaforge/bin/conda init bash
exec bash
```

```bash
mamba env create -f bioinformatics-learn.yml
conda activate bioinformatics-learn
```

## Обновление окружения при добавлении зависимостей:
```bash
conda env update -f bioinformatics-learn.yml --prune
```
## Тесты:
```bash
python -m pytest -q
```

## Обновление conda base:
```bash
conda update -n base -c conda-forge conda
```

## Installing R:
*Note:* that r-essentials installs a lot of R packages, including ggplot2,
```bash
conda install rpy2 r-essentials r-gridextra
```

### Базовые пакеты R (rpy2)
```bash
conda install r-ggplot2 r-lazyeval r-gridextra rpy2
```

### Dataset
```bash
wget -nd http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.sequence.index -O /datasets/sequence.index
```