# Copyright 2021 Observational Health Data Sciences and Informatics
#
# This file is part of FeatureExtraction
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

#' Tidy covariate data in arrow
#'
#' @details
#' Normalize covariate values by dividing by the max and/or remove redundant covariates and/or remove
#' infrequent covariates. For temporal covariates, redundancy is evaluated per time ID.
#'
#' @param covariateData      An object as generated using the \code{\link{getDbCovariateData}}
#'                           function. 
#' @param minFraction        Minimum fraction of the population that should have a non-zero value for a
#'                           covariate for that covariate to be kept. Set to 0 to don't filter on
#'                           frequency.
#' @param normalize          Normalize the covariates? (dividing by the max).
#' @param removeRedundancy   Should redundant covariates be removed?
#'
#' @export
tidyCovariateData <- function(covariateData,
                              minFraction = 0.001,
                              normalize = TRUE,
                              removeRedundancy = TRUE) {
  start <- Sys.time()
  
  newCovariateData <- list(covariateRef = covariateData$covariateRef,
                           analysisRef = covariateData$analysisRef)
  metaData <- attr(covariateData, "metaData")
  populationSize <- metaData$populationSize
  if (covariateData$covariates %>% dplyr::count() %>% dplyr::pull() == 0) {
    newCovariateData$covariates <- covariateData$covariates
  } else {
    newCovariates <- covariateData$covariates
    covariateData$maxValuePerCovariateId <- covariateData$covariates %>% 
      dplyr::group_by(.data$covariateId) %>% 
      dplyr::summarise(maxValue = max(.data$covariateValue, na.rm = TRUE))
    on.exit(covariateData$maxValuePerCovariateId <- NULL)
    
    if (removeRedundancy || minFraction != 0) {
      covariateData$valueCounts <- covariateData$covariates %>% 
        dplyr::group_by(.data$covariateId) %>% 
        dplyr::summarise(n = n(), nDistinct = n_distinct(.data$covariateValue))
      on.exit(covariateData$valueCounts <- NULL, add = TRUE)
    }
    
    ignoreCovariateIds <- c()
    deleteCovariateIds <- c()
    if (removeRedundancy) {
      covariateData$binaryCovariateIds <- covariateData$maxValuePerCovariateId %>% 
        dplyr::inner_join(covariateData$valueCounts, by = "covariateId") %>%
        dplyr::filter(.data$maxValue == 1 & .data$nDistinct == 1) %>%
        dplyr::select(covariateId = .data$covariateId)
      on.exit(covariateData$binaryCovariateIds <- NULL, add = TRUE)
      
      if (covariateData$binaryCovariateIds %>% dplyr::count() %>% dplyr::pull() != 0) {
        if (!is.null(covariateData$covariates$schema$GetFieldByName('timeId'))) { 
          # Temporal
          covariateData$temporalValueCounts <- covariateData$covariates %>% 
            dplyr::inner_join(covariateData$binaryCovariateIds, by = "covariateId") %>% 
            dplyr::group_by(.data$covariateId, .data$timeId) %>% 
            dplyr::count()
          on.exit(covariateData$temporalValueCounts <- NULL, add = TRUE)
          
          # First, find all single covariates that, for every timeId, appear in every row with the same value
          covariateData$deleteCovariateTimeIds <-  covariateData$temporalValueCounts %>% 
            dplyr::filter(n == populationSize) %>% 
            dplyr::select(.data$covariateId, .data$timeId)
          on.exit(covariateData$deleteCovariateTimeIds <- NULL, add = TRUE)
          
          # Next, find groups of covariates (analyses) that together cover everyone:
          analysisIds <- covariateData$temporalValueCounts %>%
            dplyr::anti_join(covariateData$deleteCovariateTimeIds, by = c("covariateId", "timeId")) %>%
            dplyr::inner_join(covariateData$covariateRef, by = "covariateId") %>%
            dplyr::group_by(.data$analysisId) %>%
            dplyr::summarise(n = sum(.data$n, na.rm = TRUE)) %>%
            dplyr::filter(n == populationSize) %>% 
            dplyr::select(.data$analysisId) 
          
          # For those, find most prevalent covariate, and mark it for deletion:
          valueCounts <- analysisIds %>%
            dplyr::inner_join(covariateData$covariateRef, by = "analysisId") %>%
            dplyr::inner_join(covariateData$temporalValueCounts, by = "covariateId") %>%
            dplyr::select(.data$analysisId, .data$covariateId, .data$timeId, .data$n) %>%
          valueCounts <- valueCounts[order(valueCounts$analysisId, -valueCounts$n), ]
          Andromeda::appendToTable(covariateData$deleteCovariateTimeIds, 
                                   valueCounts[!duplicated(valueCounts$analysisId), c("covariateId", "timeId")])
          
          newCovariates <- newCovariates %>%
            dplyr::anti_join(covariateData$deleteCovariateTimeIds, by = c("covariateId", "timeId"))
          
          ParallelLogger::logInfo("Removing ", covariateData$deleteCovariateTimeIds  %>% dplyr::count() %>% dplyr::pull(), " redundant covariate ID - time ID combinations")
        } else {
          # Non-temporal
          
          # First, find all single covariates that appear in every row with the same value
          toDelete <-  covariateData$valueCounts %>% 
            dplyr::inner_join(covariateData$binaryCovariateIds, by = "covariateId") %>% 
            dplyr::filter(n == populationSize) %>% 
            dplyr::select(.data$covariateId) %>% dplyr::collect()
          deleteCovariateIds <- toDelete$covariateId
          
          if (is.null(deleteCovariateIds)) {
            deleteCovariateIds <- c(-1)
          }
          # Next, find groups of covariates (analyses) that together cover everyone:
          analysisIds <- covariateData$valueCounts %>%
            dplyr::inner_join(covariateData$binaryCovariateIds, by = "covariateId") %>% 
            dplyr::filter(!.data$covariateId %in% deleteCovariateIds) %>%
            dplyr::inner_join(covariateData$covariateRef, by = "covariateId") %>%
            dplyr::group_by(.data$analysisId) %>%
            dplyr::summarise(n = sum(.data$n, na.rm = TRUE)) %>%
            dplyr::filter(n == populationSize) %>% 
            dplyr::select(.data$analysisId)
          # For those, find most prevalent covariate, and mark it for deletion:
          valueCounts <- analysisIds %>%
            dplyr::inner_join(covariateData$covariateRef, by = "analysisId") %>%
            dplyr::inner_join(covariateData$valueCounts, by = "covariateId") %>%
            dplyr::select(.data$analysisId, .data$covariateId, .data$n) %>% dplyr::collect()
          valueCounts <- valueCounts[order(valueCounts$analysisId, -valueCounts$n), ]
          deleteCovariateIds <- c(deleteCovariateIds, valueCounts$covariateId[!duplicated(valueCounts$analysisId)])
          ignoreCovariateIds <- valueCounts$covariateId
          ParallelLogger::logInfo("Removing ", length(deleteCovariateIds), " redundant covariates")
        }
      }
      metaData$deletedRedundantCovariateIds <- deleteCovariateIds
    }
    if (minFraction != 0) {
      minCount <- floor(minFraction * populationSize)
      toDelete <- covariateData$valueCounts %>%
        dplyr::filter(.data$n < minCount) %>%
        dplyr::filter(!.data$covariateId %in% ignoreCovariateIds) %>%
        dplyr::select(.data$covariateId) %>%
        dplyr::collect()
      
      metaData$deletedInfrequentCovariateIds <- toDelete$covariateId
      deleteCovariateIds <- c(deleteCovariateIds, toDelete$covariateId)
      ParallelLogger::logInfo("Removing ", nrow(toDelete), " infrequent covariates")
    }
    if (length(deleteCovariateIds) > 0) {
      newCovariates <- newCovariates %>% dplyr::filter(!.data$covariateId %in% deleteCovariateIds)
      newCovariates <- createArrow(newCovariates)
    }
    if (normalize) {
      ParallelLogger::logInfo("Normalizing covariates")
      newCovariates <- newCovariates %>% 
        dplyr::inner_join(covariateData$maxValuePerCovariateId, by = "covariateId") %>%
        dplyr::mutate(covariateValue = .data$covariateValue / .data$maxValue) %>%
        dplyr::select(-.data$maxValue)
      newCovariates <- createArrow(newCovariates)
      metaData$normFactors <- covariateData$maxValuePerCovariateId %>%
        dplyr::collect()
    } 
    newCovariateData$covariates <- newCovariates
  }
  
  class(newCovariateData) <- "CovariateData"
  attr(class(newCovariateData), "package") <- "FeatureExtraction"
  attr(newCovariateData, "metaData") <- metaData
  
  delta <- Sys.time() - start
  ParallelLogger::logInfo("Tidying covariates took ", signif(delta, 3), " ", attr(delta, "units"))
  
  return(newCovariateData)
}

# Copyright 2021 Observational Health Data Sciences and Informatics
#
# This file is part of FeatureExtraction
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

#' Tidy covariate data in memory
#'
#' @details
#' Normalize covariate values by dividing by the max and/or remove redundant covariates and/or remove
#' infrequent covariates. For temporal covariates, redundancy is evaluated per time ID.
#'
#' @param covariateData      An object as generated using the \code{\link{getDbCovariateData}}
#'                           function. 
#' @param minFraction        Minimum fraction of the population that should have a non-zero value for a
#'                           covariate for that covariate to be kept. Set to 0 to don't filter on
#'                           frequency.
#' @param normalize          Normalize the covariates? (dividing by the max).
#' @param removeRedundancy   Should redundant covariates be removed?
#'
#' @export
tidyCovariateDataMemory <- function(covariateData,
                                   minFraction = 0.001,
                                   normalize = TRUE,
                                   removeRedundancy = TRUE) {
  start <- Sys.time()
  
  newCovariateData <- list(covariateRef = covariateData$covariateRef,
                           analysisRef = covariateData$analysisRef)
  metaData <- attr(covariateData, "metaData")
  populationSize <- metaData$populationSize
  if (covariateData$covariates %>% dplyr::count() %>% dplyr::pull() == 0) {
    newCovariateData$covariates <- covariateData$covariates
  } else {
    newCovariates <- covariateData$covariates
    covariateData$maxValuePerCovariateId <- covariateData$covariates %>% 
      dplyr::group_by(.data$covariateId) %>% 
      dplyr::summarise(maxValue = max(.data$covariateValue, na.rm = TRUE)) %>% 
      dplyr::collect()
    on.exit(covariateData$maxValuePerCovariateId <- NULL)
    
    if (removeRedundancy || minFraction != 0) {
      covariateData$valueCounts <- covariateData$covariates %>% 
        dplyr::group_by(.data$covariateId) %>% 
        dplyr::summarise(n = n(), nDistinct = n_distinct(.data$covariateValue)) %>%
        dplyr::collect()
      on.exit(covariateData$valueCounts <- NULL, add = TRUE)
    }
    
    ignoreCovariateIds <- c()
    deleteCovariateIds <- c()
    if (removeRedundancy) {
      covariateData$binaryCovariateIds <- covariateData$maxValuePerCovariateId %>% 
        dplyr::inner_join(covariateData$valueCounts, by = "covariateId") %>%
        dplyr::filter(.data$maxValue == 1 & .data$nDistinct == 1) %>%
        dplyr::select(covariateId = .data$covariateId)
      on.exit(covariateData$binaryCovariateIds <- NULL, add = TRUE)
      
      if (covariateData$binaryCovariateIds %>% dplyr::count() %>% dplyr::pull() != 0) {
        if ('timeId' %in% names(covariateData$covariates)) { 
          # Temporal
          covariateData$temporalValueCounts <- covariateData$covariates %>% 
            dplyr::inner_join(covariateData$binaryCovariateIds, by = "covariateId") %>% 
            dplyr::group_by(.data$covariateId, .data$timeId) %>% 
            dplyr::count()
          on.exit(covariateData$temporalValueCounts <- NULL, add = TRUE)
          
          # First, find all single covariates that, for every timeId, appear in every row with the same value
          covariateData$deleteCovariateTimeIds <-  covariateData$temporalValueCounts %>% 
            dplyr::filter(n == populationSize) %>% 
            dplyr::select(.data$covariateId, .data$timeId)
          on.exit(covariateData$deleteCovariateTimeIds <- NULL, add = TRUE)
          
          # Next, find groups of covariates (analyses) that together cover everyone:
          analysisIds <- covariateData$temporalValueCounts %>%
            dplyr::anti_join(covariateData$deleteCovariateTimeIds, by = c("covariateId", "timeId")) %>%
            dplyr::inner_join(covariateData$covariateRef, by = "covariateId") %>%
            dplyr::group_by(.data$analysisId) %>%
            dplyr::summarise(n = sum(.data$n, na.rm = TRUE)) %>%
            dplyr::filter(n == populationSize) %>% 
            dplyr::select(.data$analysisId) 
          
          # For those, find most prevalent covariate, and mark it for deletion:
          valueCounts <- analysisIds %>%
            dplyr::inner_join(covariateData$covariateRef, by = "analysisId") %>%
            dplyr::inner_join(covariateData$temporalValueCounts, by = "covariateId") %>%
            dplyr::select(.data$analysisId, .data$covariateId, .data$timeId, .data$n) %>%
            dplyr::collect()
          valueCounts <- valueCounts[order(valueCounts$analysisId, -valueCounts$n), ]
          Andromeda::appendToTable(covariateData$deleteCovariateTimeIds, 
                                   valueCounts[!duplicated(valueCounts$analysisId), c("covariateId", "timeId")])
          
          newCovariates <- newCovariates %>%
            dplyr::anti_join(covariateData$deleteCovariateTimeIds, by = c("covariateId", "timeId"))
          
          ParallelLogger::logInfo("Removing ", covariateData$deleteCovariateTimeIds  
                                  %>% dplyr::count() %>% dplyr::pull(), " redundant covariate ID - time ID combinations")
        } else {
          # Non-temporal
          
          # First, find all single covariates that appear in every row with the same value
          toDelete <-  covariateData$valueCounts %>% 
            dplyr::inner_join(covariateData$binaryCovariateIds, by = "covariateId") %>% 
            dplyr::filter(n == populationSize) %>% 
            dplyr::select(.data$covariateId) %>% 
            dplyr::collect()
          deleteCovariateIds <- toDelete$covariateId
          
          # Next, find groups of covariates (analyses) that together cover everyone:
          analysisIds <- covariateData$valueCounts %>%
            dplyr::inner_join(covariateData$binaryCovariateIds, by = "covariateId") %>% 
            dplyr::filter(!.data$covariateId %in% deleteCovariateIds) %>%
            dplyr::inner_join(covariateData$covariateRef, by = "covariateId") %>%
            dplyr::group_by(.data$analysisId) %>%
            dplyr::summarise(n = sum(.data$n, na.rm = TRUE)) %>%
            dplyr::filter(n == populationSize) %>% 
            dplyr::select(.data$analysisId) 
          # For those, find most prevalent covariate, and mark it for deletion:
          valueCounts <- analysisIds %>%
            dplyr::inner_join(covariateData$covariateRef, by = "analysisId") %>%
            dplyr::inner_join(covariateData$valueCounts, by = "covariateId") %>%
            dplyr::select(.data$analysisId, .data$covariateId, .data$n) %>%
            dplyr::collect()
          valueCounts <- valueCounts[order(valueCounts$analysisId, -valueCounts$n), ]
          deleteCovariateIds <- c(deleteCovariateIds, valueCounts$covariateId[!duplicated(valueCounts$analysisId)])
          ignoreCovariateIds <- valueCounts$covariateId
          ParallelLogger::logInfo("Removing ", length(deleteCovariateIds), " redundant covariates")
        }
      }
      metaData$deletedRedundantCovariateIds <- deleteCovariateIds
    }
    if (minFraction != 0) {
      minCount <- floor(minFraction * populationSize)
      toDelete <- covariateData$valueCounts %>%
        dplyr::filter(.data$n < minCount) %>%
        dplyr::filter(!.data$covariateId %in% ignoreCovariateIds) %>%
        dplyr::select(.data$covariateId) %>%
        dplyr::collect()
      
      metaData$deletedInfrequentCovariateIds <- toDelete$covariateId
      deleteCovariateIds <- c(deleteCovariateIds, toDelete$covariateId)
      ParallelLogger::logInfo("Removing ", nrow(toDelete), " infrequent covariates")
    }
    if (length(deleteCovariateIds) > 0) {
      newCovariates <- newCovariates %>% dplyr::filter(!.data$covariateId %in% deleteCovariateIds) 
    }
    if (normalize) {
      ParallelLogger::logInfo("Normalizing covariates")
      newCovariates <- newCovariates %>% 
        dplyr::inner_join(covariateData$maxValuePerCovariateId, by = "covariateId") %>%
        dplyr::mutate(covariateValue = .data$covariateValue / .data$maxValue) %>%
        dplyr::select(-.data$maxValue)
      metaData$normFactors <- covariateData$maxValuePerCovariateId %>%
        dplyr::collect()
    } 
    newCovariateData$covariates <- newCovariates
  }
  
  class(newCovariateData) <- "CovariateData"
  attr(class(newCovariateData), "package") <- "FeatureExtraction"
  attr(newCovariateData, "metaData") <- metaData
  
  delta <- Sys.time() - start
  ParallelLogger::logInfo("Tidying covariates took ", signif(delta, 3), " ", attr(delta, "units"))
  return(newCovariateData)
}

