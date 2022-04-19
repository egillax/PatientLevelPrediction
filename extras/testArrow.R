library(PatientLevelPredictionArrow)

arrow <- T
data(plpDataSimulationProfile)
sampleSize <- 1e5
plpData <- simulatePlpData(
  plpDataSimulationProfile,
  n = sampleSize
)

populationSet <- PatientLevelPredictionArrow::createStudyPopulationSettings(
  requireTimeAtRisk = F, 
  riskWindowStart = 1, 
  riskWindowEnd = 365)


modelSettings <- PatientLevelPredictionArrow::setLassoLogisticRegression()

if (arrow) {
  res2 <- runPlp(
    plpData = plpData,
    outcomeId = 3,
    modelSettings = modelSettings,
    analysisId = 'Test',
    analysisName = 'Testing Arrow',
    populationSettings = populationSet,
    splitSettings = createDefaultSplitSetting(),
    sampleSettings = createSampleSettings(),  # none
    featureEngineeringSettings = createFeatureEngineeringSettings(), # none
    preprocessSettings = createPreprocessSettings(),
    logSettings = createLogSettings(),
    executeSettings = createExecuteSettings(
      runSplitData = T,
      runSampleData = F,
      runfeatureEngineering = F,
      runPreprocessData = T,
      runModelDevelopment = T,
      runCovariateSummary = T
    ),
    saveDirectory = '~/test/arrow_new_plp/'
  )
} else {
  
library(PatientLevelPrediction)
res2 <- PatientLevelPrediction::runPlp(
    plpData = plpData,
    outcomeId = 3,
    modelSettings = modelSettings,
    analysisId = 'Test',
    analysisName = 'Testing Original',
    populationSettings = populationSet,
    splitSettings = createDefaultSplitSetting(),
    sampleSettings = createSampleSettings(),  # none
    featureEngineeringSettings = createFeatureEngineeringSettings(), # none
    preprocessSettings = createPreprocessSettings(),
    logSettings = createLogSettings(),
    executeSettings = createExecuteSettings(
      runSplitData = T,
      runSampleData = F,
      runfeatureEngineering = F,
      runPreprocessData = T,
      runModelDevelopment = T,
      runCovariateSummary = T
    ),
    saveDirectory = '~/test/original_plp/'
  )
}

