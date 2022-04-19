plpAndromeda <- PatientLevelPrediction::simulatePlpData(plpDataSimulationProfile,
                                                        n=sampleSize)
start <- Sys.time()
arrowDataset <- convertAndromedaToArrow(plpAndromeda$covariateData$covariates, path)
delta <- Sys.time() - start


#' @description
#' converts an andromeda table to arrow dataset and writes it to the specified path
#' @param andromeda andromeda table such as plpData$covariateData$covariates
#' @param path      where to save arrow dataset
convertAndromedaToArrow <- function(andromeda, path) {
  newPath <- file.path(tempdir(), paste(c('arrowDataset', sample(letters, 8)), collapse = ''))
  dir.create(newPath)
  fun <- function(x, path) {
    arrow::write_feather(x,file.path(path, paste(c('file_', sample(letters, 8)), collapse = '')))
      }
  Andromeda::batchApply(andromeda, fun, path=newPath, batchSize = 1e6,
                        progressBar = TRUE)
  
  
  tempDataset <- arrow::open_dataset(newPath, format='feather')
  unlink(path, recursive = TRUE)
  
  arrow::write_dataset(tempDataset, path = path, format='feather')
  
  dataset <- arrow::open_dataset(path, format='feather')  
  return(dataset)
}

