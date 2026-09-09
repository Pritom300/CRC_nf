FROM rocker/tidyverse:4.4.1

RUN apt-get update && apt-get install -y \
    cmake libuv1-dev libcurl4-openssl-dev libssl-dev \
    libxml2-dev libfontconfig1-dev libharfbuzz-dev \
    libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev sqlite3 libsqlite3-dev \
    libcairo2-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# STEP 1: Install CRAN packages FIRST (Using the correct cloud.r-project.org repository)
RUN R -e "install.packages(c( \
    'DBI', 'RSQLite', 'rentrez', 'BiocManager', 'survival', \
    'survminer', 'ggpubr', 'gridExtra', 'glmnet', \
    'randomForest', 'caret', 'reshape2', 'recipes', 'e1071', 'pROC' \
), repos='https://cloud.r-project.org', dependencies=TRUE, Ncpus=4)"

# Copy all pre-downloaded Bioconductor packages
COPY docker_packages/ /tmp/bioc_packages/

# STEP 2: Install Bioconductor packages from local files
RUN R -e "install.packages('/tmp/bioc_packages/BiocGenerics_0.52.0.tar.gz',       repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/Biobase_2.66.0.tar.gz',            repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/S4Vectors_0.44.0.tar.gz',          repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/IRanges_2.40.1.tar.gz',            repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/zlibbioc_1.52.0.tar.gz',           repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/XVector_0.46.0.tar.gz',            repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/GenomeInfoDbData_1.2.13.tar.gz',   repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/UCSC.utils_1.2.0.tar.gz',          repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/GenomeInfoDb_1.42.3.tar.gz',       repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/Biostrings_2.74.1.tar.gz',         repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/KEGGREST_1.46.0.tar.gz',           repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/MatrixGenerics_1.18.1.tar.gz',     repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/S4Arrays_1.6.0.tar.gz',            repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/SparseArray_1.6.2.tar.gz',         repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/DelayedArray_0.32.0.tar.gz',       repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/GenomicRanges_1.58.0.tar.gz',      repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/AnnotationDbi_1.68.0.tar.gz',      repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/org.Hs.eg.db_3.20.0.tar.gz',       repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/hgu133plus2.db_3.13.0.tar.gz',     repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/SummarizedExperiment_1.36.0.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/limma_3.62.2.tar.gz',              repos=NULL, type='source')"
RUN R -e "install.packages('/tmp/bioc_packages/GEOquery_2.74.0.tar.gz',           repos=NULL, type='source')"

# STEP 3: Loud validation health check (Fixed syntax error)
RUN R -e " \
  pkgs <- c('hgu133plus2.db','limma','GEOquery','AnnotationDbi', \
            'Biobase','BiocGenerics','org.Hs.eg.db', \
            'GenomicRanges','SummarizedExperiment','IRanges', \
            'tidyverse','randomForest','caret','survival','pROC'); \
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
  if (length(missing)) stop('MISSING PACKAGES: ', paste(missing, collapse=' ')); \
  cat('ALL PACKAGES OK\n') \
"

WORKDIR /project
