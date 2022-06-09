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
#' convertArrowToCyclopsData
#' @describeIn convertArrowToCyclopsData Convert data from two \code{Dataset} tables
#' @export
convertArrowToCyclopsData <- function(outcomes,
                                      covariates,
                                      modelType = "lr",
                                      addIntercept = TRUE,
                                      checkSorting = NULL,
                                      checkRowIds = TRUE,
                                      normalize = NULL,
                                      quiet = FALSE,
                                      floatingPoint = 64) {
    class(outcomes) <- "data.frame"
    if (!is.null(checkSorting))
      warning("The 'checkSorting' argument has been deprecated. Sorting is now always checked")
    
    if ((modelType == "clr" | modelType == "cpr") & addIntercept) {
      if (!quiet) {
        warning("Intercepts are not allowed in conditional models, removing intercept", call. = FALSE)
      }
      addIntercept = FALSE
    }
    
    if (modelType == "pr" | modelType == "cpr") {
      if (any(outcomes$time <= 0)) {
        stop("time cannot be non-positive", call. = FALSE)
      }
    }
    
    providedNoStrata <- !"stratumId" %in% colnames(outcomes)
    
    if (modelType == "cox" | modelType == "fgr") {
      if (providedNoStrata) {
        outcomes <- outcomes %>%
          mutate(stratumId = 0)
        covariates <- covariates %>%
          mutate(stratumId = 0)
      }
    }
    
    if (checkRowIds) {
      covariateRowIds <- covariates %>%
        distinct(.data$rowId) %>%
        pull()
      outcomeRowIds <- select(outcomes, .data$rowId) %>%
        pull()
      mapping <- match(covariateRowIds, outcomeRowIds)
      if (any(is.na(mapping))) {
        if (!quiet) {
          writeLines("Removing covariate values with rowIds that are not in outcomes")
        }
        covariates <- covariates %>%
          dplyr::filter(.data$rowId %in% outcomeRowIds)
      }
    }
    
    # Sorting should be last, as other operations may change ordering.
    # Also, should always explicitly define sorting, else not guaranteed.]
    if (modelType == "lr" | modelType == "pr") {
      outcomes <- outcomes %>%
        dplyr::arrange(.data$rowId)
      
      covariates <- covariates %>%
        dplyr::arrange(.data$covariateId, .data$rowId)
    }
    
    if (modelType == "clr" | modelType == "cpr") {
      outcomes <- outcomes %>%
        dplyr::arrange(.data$stratumId, .data$rowId)
      
      covariates <- covariates %>%
        dplyr::arrange(.data$covariateId, .data$stratumId, .data$rowId)
    }
    
    if (modelType == "cox" | modelType == "fgr") {
      
      if (modelType == "cox" &
          (select(outcomes, .data$y) %>% dplyr::distinct() %>% dplyr::count() 
           %>% dplyr::collect() > 2)) {
        stop("Cox model only accepts one outcome type")
      }
      
      outcomes <- outcomes %>%
        dplyr::arrange(.data$stratumId, desc(.data$time), .data$y, .data$rowId)
      if (!"time" %in% colnames(covariates)) {
        covariates <- covariates %>%
          dplyr::inner_join(select(outcomes, .data$rowId, .data$time, .data$y), by = "rowId")
      }
      covariates <- covariates %>%
        dplyr::arrange(.data$covariateId, .data$stratumId, desc(.data$time), .data$y, .data$rowId)
    }
    
    dataPtr <- Cyclops:::createSqlCyclopsData(modelType = modelType, floatingPoint = floatingPoint)
    
    outcomes <- dplyr::collect(outcomes)
    if (modelType == "lr" | modelType == "pr") {
      outcomes$stratumId <- NULL
    }
    
    Cyclops:::loadNewSqlCyclopsDataY(object = dataPtr,
                           stratumId = if ("stratumId" %in% colnames(outcomes)) outcomes$stratumId else NULL,
                           rowId = outcomes$rowId,
                           y = outcomes$y,
                           time = if ("time" %in% colnames(outcomes)) outcomes$time else NULL)
    
    if (addIntercept & (modelType != "cox" & modelType != "fgr")) {
      Cyclops:::loadNewSqlCyclopsDataX(dataPtr, 0, NULL, NULL, name = "(Intercept)")
    }
    
    loadCovariates <- function(batch) {
      batch <- as.data.frame(batch)
      covarNames <- unique(batch$covariateId)
      Cyclops:::loadNewSqlCyclopsDataMultipleX(object = dataPtr,
                                     covariateId = batch$covariateId,
                                     rowId = batch$rowId,
                                     covariateValue = batch$covariateValue,
                                     name = covarNames,
                                     append = TRUE)
    }
    arrow::map_batches(covariates, loadCovariates)
    
    
    if (modelType == "pr" || modelType == "cpr")
      Cyclops:::finalizeSqlCyclopsData(dataPtr, useOffsetCovariate = -1)
    
    if (!is.null(normalize)) {
      .normalizeCovariates(dataPtr, normalize)
    }
    
    if ("weights" %in% colnames(outcomes)) {
      dataPtr$weights <- outcomes %>% dplyr::pull(.data$weights)
    } else {
      dataPtr$weights <- NULL
    }
    
    if ("censorWeights" %in% colnames(outcomes)) {
      dataPtr$censorWeights <- outcomes %>% dplyr::pull(.data$censorWeights)
    } else {
      if (modelType == "fgr") {
        dataPtr$censorWeights <- getFineGrayWeights(
          outcomes %>% pull(.data$time),
          outcomes %>% pull(.data$y)
        )$weights
        writeLines("Generating censoring weights")
      } else {
        dataPtr$censorWeights <- NULL
      }
    }
    return(dataPtr)
}
