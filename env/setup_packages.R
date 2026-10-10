devtools::install_github("microbiome/mia", upgrade = "never", repos = c("https://cloud.r-project.org/", BiocManager::repositories()))

devtools::install_github("FrederickHuangLin/ANCOMBC", ref = "50131c4", upgrade = "never", repos = c("https://cloud.r-project.org/", BiocManager::repositories())) # Commit before R version bump

devtools::install_github("microbiome/microbiome", upgrade = "never", repos = c("https://cloud.r-project.org/", BiocManager::repositories()))

# Since two of NetCoMi's dependencies are only available on GitHub, 
# it is recommended to install them first:
devtools::install_github("zdk123/SpiecEasi", upgrade = "never", ref = "a92f78c")
devtools::install_github("GraceYoon/SPRING", upgrade = "never")

devtools::install_github("stefpeschel/NetCoMi", repos = c("https://cloud.r-project.org/", BiocManager::repositories()))

installable <- c(
    "jbisanz/qiime2R",
    "benjjneb/dada2",
    "biobakery/maaslin3"
)

for (i in installable) {
    if (!requireNamespace(i, quietly = TRUE)) {
        devtools::install_github(i, upgrade = "never", repos = c("https://cloud.r-project.org/", BiocManager::repositories()))
    }
}

devtools::install_github("mjskay/ggdist", ref = "v3.3.2", upgrade = "never", repos = c("https://cloud.r-project.org/", BiocManager::repositories()))
devtools::install_github("microbiome/miaViz", upgrade = "never", repos = c("https://cloud.r-project.org/", BiocManager::repositories()))

## Reinstall, because Matrix package is from another version
install.packages("lme4", repos = "https://cloud.r-project.org/", dependencies = TRUE)

devtools::install_github("nsgrantham/ggdark", upgrade = "never", repos = c("https://cloud.r-project.org/", BiocManager::repositories()))

devtools::install_github("cafferychen777/ggpicrust2@v2.5.4")
