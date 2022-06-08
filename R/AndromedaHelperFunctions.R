# @file ffHelperFunctions.R
#
# Copyright 2021 Observational Health Data Sciences and Informatics
#
# This file is part of PatientLevelPrediction
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

limitCovariatesToPopulation <- function(covariateData, rowIds) {
  ParallelLogger::logTrace(paste0('Starting to limit covariate data to population...'))
  
  metaData <- attr(covariateData, 'metaData')
  
  newCovariateData <- list(covariateRef = covariateData$covariateRef,
                           analysisRef = covariateData$analysisRef)
  
  covariateData$pop <- data.frame(rowId = rowIds)
  
  on.exit(covariateData$pop <- NULL, add = T)
  covariates <- covariateData$covariates %>% dplyr::inner_join(covariateData$pop, by='rowId')
  newCovariateData$covariates <- createArrow(covariates)
 
  metaData$populationSize <- length(rowIds)
  attr(newCovariateData, 'metaData') <- metaData
  class(newCovariateData) <- "CovariateData"
  ParallelLogger::logTrace(paste0('Finished limiting covariate data to population...'))
  return(newCovariateData)
}


# is this used?
# return prev of ffdf 
calculatePrevs <- function(plpData, population){
  #===========================
  # outcome prevs
  #===========================
  
  # add population to sqllite
  population <- tibble::as_tibble(population)
  plpData$covariateData$population <- population %>% dplyr::select(.data$rowId, .data$outcomeCount)
  
  outCount <- nrow(plpData$covariateData$population %>% dplyr::filter(.data$outcomeCount == 1))
  nonOutCount <- nrow(plpData$covariateData$population %>% dplyr::filter(.data$outcomeCount == 0))
  
  # join covariate with label
  prevs <- plpData$covariateData$covariates %>% dplyr::inner_join(plpData$covariateData$population) %>%
    dplyr::group_by(.data$covariateId) %>% 
    dplyr::summarise(prev.out = 1.0*sum(.data$outcomeCount==1, na.rm = TRUE)/outCount,
              prev.noout = 1.0*sum(.data$outcomeCount==0, na.rm = TRUE)/nonOutCount) %>%
    dplyr::select(.data$covariateId, .data$prev.out, .data$prev.noout)
  
  #clear up data
  ##plpData$covariateData$population <- NULL
  
  return(as.data.frame(prevs))
}

# creates an arrow dataset in a temp directory from a query object, unlinks
# old query object files
#' @export
createArrow <- function(object, name=NULL, delete=FALSE) {
  UseMethod("createArrow")
}


#' @export
createArrow.default <- function(object, name=NULL, delete=FALSE) {
  if (is.null(name)) {
    name <- deparse(substitute(object))
  }
  
  getDataObject <- function(object) {
    if (class(object)[[1]] == 'arrow_dplyr_query') {
      data <- object$.data
      data <- getDataObject(data)
    } else {
      return(object)
    }
    return(data)
  }
  
  dataObject <- getDataObject(object)
  oldFiles <- dataObject$files
  oldDir <- dirname(oldFiles[[1]])
  newPath <- file.path(tempdir(), paste(c(name,'_', sample(letters,8)), collapse = ''))
  arrow::write_dataset(object, path=newPath, format='feather')
  dataset <- arrow::open_dataset(newPath, format='feather')
  
  if (delete){
  # clean up old files, use with caution
  for (file in oldFiles) {
    unlink(file)
  }
  unlink(oldDir, recursive = TRUE)
  }
  return(dataset)
}

#' @export
createArrow.data.frame <- function(object, name=NULL) {
  if (is.null(name)) {
    name <- deparse(substitute(object))
  }
  newPath <- file.path(tempdir(), paste(c(name,'_', sample(letters,8)), collapse = ''))
  arrow::write_dataset(object, path=newPath, format='feather')
  dataset <- arrow::open_dataset(newPath, format='feather')
  return(dataset)
}

#' convertCovariateDataToArrow
#' @description Converts a covariateData object from using a rsqlite backend
#' to arrow backend
#' @param covariateData  A covariateData object from featureExtraction
#' @export
convertCovariateDataToArrow <- function(covariateData) {
  newCovariateData <- list(covariateRef=covariateData$covariateRef %>% dplyr::collect(),
                           analysisRef=covariateData$analysisRef %>% dplyr::collect())
  
  newCovariateData$covariates <- convertAndromedaToArrow(plpData$covariateData$covariates)
  
  metaData <- attr(plpData$covariateData, 'metaData') 
  
  attr(newCovariateData, 'metaData') <- metaData
  class(newCovariateData) <- 'CovariateData'
  
  return(newCovariateData)
}


#' convertAndromedaToArrow
#' @description
#' converts an andromeda table to arrow dataset and writes it to the specified path
#' @param andromeda andromeda table such as plpData$covariateData$covariates
#' @param path      where to save arrow dataset
#' @export
convertAndromedaToArrow <- function(andromeda, path) {
  newPath <- file.path(tempdir(), paste(c('arrowDataset', sample(letters, 8)), collapse = ''))
  dir.create(newPath)
  fun <- function(x, path) {
    arrow::write_feather(x,file.path(path, paste(c('file_', sample(letters, 8)), collapse = '')))
  }
  Andromeda::batchApply(andromeda, fun, path=newPath, batchSize = 1e6,
                        progressBar = TRUE)
  
  
  tempDataset <- arrow::open_dataset(newPath, format='feather')
  # unlink(path, recursive = TRUE)
  newPath <- file.path(tempdir(), paste(c('arrowDataset', sample(letters, 8)), collapse = ''))
  arrow::write_dataset(tempDataset, path = newPath, format='feather')
  
  dataset <- arrow::open_dataset(newPath, format='feather')  
  return(dataset)
}