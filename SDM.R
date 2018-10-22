workdir <- "E:/Lab Jobs/Conservation Biology&Ecology/Pheasant/CROSS BREEDING"
setwd(workdir)
allPA <- read.csv('all_PA_proj_thinned.csv')
redPA <- allPA[allPA$species=='红腹锦鸡',]
whitePA <- allPA[allPA$species=='白腹锦鸡',]

require(raster)
require(biomod2)

pathfac <- paste0(workdir,'/SDM/rasters_factor')
pathcon <- paste0(workdir,'/SDM/rasters_con')

rasters_factor <- list.files(path = pathfac)
rasters_con <- list.files(path = pathcon)

stack1 <- stack((paste(pathfac,rasters_factor,sep = '/')))
stack1[[2]] <- stack1[[2]] * (stack1[[1]]<=2)
stack1[[1]] <- asFactor(stack1[[1]])
stack1[[2]] <- asFactor(stack1[[2]])
stack1[[1]] <- ratify(stack1[[1]])
stack1[[2]] <- ratify(stack1[[2]])

stack2 <- stack(paste(pathcon,rasters_con,sep='/'))
my.stack <- stack(stack1,stack2)

rm(stack1)
rm(stack2)

distributiondata.species1 <- 
  BIOMOD_FormatingData(
    resp.var = as.numeric(redPA[,4]),
    expl.var = my.stack,
    resp.xy = (redPA[, 2:3]),
    resp.name = 'red'
  )

distributiondata.species2 <- 
  BIOMOD_FormatingData(
    resp.var = as.numeric(whitePA[,4]),
    expl.var = my.stack,
    resp.xy = (whitePA[, 2:3]),
    resp.name = 'white'
  )

SDM.1 <-
  BIOMOD_Modeling(
    data = distributiondata.species1,
    models = 'RF',
    DataSplit=85,
    models.eval.meth = c('TSS', 'ROC') ,
    VarImport = 3,
    SaveObj = T
  )

SDM.2 <-
  BIOMOD_Modeling(
    data = distributiondata.species2,
    models = 'RF',
    DataSplit=85,
    models.eval.meth = c('TSS', 'ROC') ,
    VarImport = 3,
    SaveObj = T
  )
distribution.est.1.full <- BIOMOD_Projection(SDM.1, my.stack, 'distribution1.est',build.clamping.mask = F)
plot(distribution.est.1.full)
distribution.est.1 <- distribution.est.1.full@proj@val[[2]]
plot(distribution.est.1)
#writeRaster(distribution.est.1,'red.tif',overwrite=T)

distribution.est.2.full <- BIOMOD_Projection(SDM.2, my.stack, 'distribution2.est',build.clamping.mask = F)
plot(distribution.est.2.full)
distribution.est.2 <- distribution.est.2.full@proj@val[[2]]
plot(distribution.est.2)
#writeRaster(distribution.est.2,'white.tif',overwrite=T)

