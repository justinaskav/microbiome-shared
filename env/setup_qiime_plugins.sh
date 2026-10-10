## Check if current environment contains 'qiime2'
if [[ $CONDA_DEFAULT_ENV != "qiime2"* ]]; then
    echo "Please activate the QIIME2 environment before running this script."
    exit 1
fi

if [[ "$OSTYPE" == "darwin"* ]]; then
    conda config --env --set subdir osx-64
fi

mamba config --add channels defaults
mamba config --add channels bioconda
mamba config --add channels conda-forge
# conda config --set channel_priority strict

# echo "Installing q2-krona..."
# conda install -y -c bioconda krona
# pip install git+https://github.com/kaanb93/q2-krona.git

echo "Installing ITSxpress..."
mamba install ITSxpress

mamba install -y \
   -c https://packages.qiime2.org/qiime2/${QIIME_ENV_VERSION:?set QIIME_ENV_VERSION, e.g. 2026.1}/moshpit/released/ \
   -c conda-forge -c bioconda -c defaults \
   q2-fondue

echo "Refreshing QIIME2 plugins..."
qiime dev refresh-cache
