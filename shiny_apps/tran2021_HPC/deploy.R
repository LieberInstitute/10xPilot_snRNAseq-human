library('rsconnect')
load(here::here("shiny_apps", ".deploy_info.Rdata"))
rsconnect::setAccountInfo(name=deploy_info$name, token=deploy_info$token,
    secret=deploy_info$secret)
options(repos = BiocManager::repositories())
rsconnect::deployApp(appFiles = c('app.R', "sce_hpc_small.rds", "cell_colors_hpc.rds"),
    appName = 'tran2021_HPC', account = 'libd', server = 'shinyapps.io')
