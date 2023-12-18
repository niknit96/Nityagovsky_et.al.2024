To repeat the data from this article, follow these steps:

1. Install mamba or conda on Linux:
* mamba: https://github.com/mamba-org/mamba (more faster than conda)
* conda: https://conda.io/projects/conda/en/latest/user-guide/install/linux.html
2. Download code and data from Github:
```
wget https://github.com/niknit96/Nityagovsky_et.al.2024/archive/master.zip
unzip master.zip
```
3. Create and activate conda environment for download SRA files and pre-trained classifiers for QIIME 2 Scikit-learn algorithm:
```
mamba env create --file ./Nityagovsky_et.al.2024-main/sra-tools.yml || conda env create --file ./Nityagovsky_et.al.2024-main/sra-tools.yml
mamba activate sra-tools || conda activate sra-tools
```
4. Run bash script for download SRA files and pre-trained classifiers for QIIME 2 Scikit-learn algorithm:
```
bash ./Nityagovsky_et.al.2024-main/download.sh
```
5. Create and activate conda environment for analysis data:
```
mamba env create --file ./Nityagovsky_et.al.2024-main/Nityagovsky_et.al.2024.yml || conda env create --file ./Nityagovsky_et.al.2024-main/Nityagovsky_et.al.2024.yml
mamba activate Nityagovsky_et.al.2024 || conda activate Nityagovsky_et.al.2024
```
6. Download and install qiime2R, tidyterra, sf and maptiles packages (accessed on 18 December 2023):
```
Rscript -e 'devtools::install_github("jbisanz/qiime2R")'
Rscript -e 'devtools::install_github("dieghernan/tidyterra")'
Rscript -e 'devtools::install_github("r-spatial/sf")'
Rscript -e 'devtools::install_github("riatelab/maptiles")'
```
7. Run bash script for analysis data:
```
bash ./Nityagovsky_et.al.2024-main/analysis.sh
```
8. Results in ./Nityagovsky_et.al.2024-main.
