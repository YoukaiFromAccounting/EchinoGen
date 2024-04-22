### R code from vignette source 'ASEB.Rnw'

###################################################
### code chunk number 1: ASEB.Rnw:78-81
###################################################
  library(ASEB)
  ff <- system.file("extdata", "background_sites.fa", package="ASEB")
  readSequence(ff)


###################################################
### code chunk number 2: ASEB.Rnw:89-94
###################################################
  backgroundSites <- readSequence(system.file("extdata", "background_sites.fa", package="ASEB")) 
  prodefinedSites <- readSequence(system.file("extdata", "predefined_sites.fa", package="ASEB"))
  testSites <- readSequence(system.file("extdata", "sites_to_test.fa", package="ASEB"))
  resultList <- asebSites(backgroundSites, prodefinedSites, testSites, permutationTimes=100)
  resultList$results[1:2,]


###################################################
### code chunk number 3: ASEB.Rnw:101-105
###################################################
  backgroundSitesFile <- system.file("extdata", "background_sites.fa", package="ASEB")
  prodefinedSitesFile <- system.file("extdata", "predefined_sites.fa", package="ASEB")
  testSitesFile <- system.file("extdata", "sites_to_test.fa", package="ASEB")
  asebSites(backgroundSitesFile, prodefinedSitesFile, testSitesFile, permutationTimes=100)


###################################################
### code chunk number 4: ASEB.Rnw:113-114
###################################################
  drawEScurve(resultList$curveInfo, max_p_value=0.1, min_es=0.1)


###################################################
### code chunk number 5: ASEB.Rnw:124-129
###################################################
  backgroundSites <- readSequence(system.file("extdata", "background_sites.fa", package="ASEB")) 
  prodefinedSites <- readSequence(system.file("extdata", "predefined_sites.fa", package="ASEB"))
  testProteins <- readSequence(system.file("extdata", "proteins_to_test.fa", package="ASEB"))
  resultList <- asebProteins(backgroundSites, prodefinedSites, testProteins, permutationTimes=100)
  resultList$results[1:2,]


###################################################
### code chunk number 6: ASEB.Rnw:134-138
###################################################
  backgroundSitesFile <- system.file("extdata", "background_sites.fa", package="ASEB")
  prodefinedSitesFile <- system.file("extdata", "predefined_sites.fa", package="ASEB")
  testProteinsFile <- system.file("extdata", "sites_to_test.fa", package="ASEB")
  asebProteins(backgroundSitesFile, prodefinedSitesFile, testProteinsFile, permutationTimes=100)


###################################################
### code chunk number 7: ASEB.Rnw:146-147
###################################################
  drawStat(curveInfoDataFrame=resultList$curveInfo);


