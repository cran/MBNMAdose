## ----setup, include = FALSE---------------------------------------------------
library(MBNMAdose)
#devtools::load_all()
library(rmarkdown)
library(knitr)
library(dplyr)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  include=TRUE,
  tidy.opts=list(width.cutoff=80),
  tidy=TRUE
)

## ---- echo=FALSE--------------------------------------------------------------
kable(head(triptans), digits=2) 

## ---- echo=TRUE---------------------------------------------------------------
kable(head(ssri), digits=2) 

## ---- echo=FALSE--------------------------------------------------------------
kable(head(gout), digits=2) 

## ---- echo=FALSE--------------------------------------------------------------
kable(head(osteopain), digits=2) 

## ---- echo=FALSE--------------------------------------------------------------
kable(head(alog_pcfb), digits=2) 

## -----------------------------------------------------------------------------
# Using the triptans dataset
network <- mbnma.network(triptans)
summary(network)

## ---- message=FALSE, warning=FALSE--------------------------------------------
# Prepare data using the gout dataset
goutnet <- mbnma.network(gout)
summary(goutnet)

## -----------------------------------------------------------------------------
plot(goutnet, label.distance = 5)

## -----------------------------------------------------------------------------
# Plot at the agent-level
plot(goutnet, level="agent", label.distance = 6)

## -----------------------------------------------------------------------------
# Plot connections to placebo via a two-parameter dose-response function (e.g. Emax)
plot(goutnet, level="agent", doselink = 2, remove.loops = TRUE, label.distance = 6)

## ---- results="hide"----------------------------------------------------------
# Colour vertices by agent
plot(goutnet, v.color = "agent", label.distance = 5)

## ---- results="hide", message=FALSE, warning=FALSE----------------------------
# Run a random effect split NMA using the alogliptin dataset
alognet <- mbnma.network(alog_pcfb)
nma.alog <- nma.run(alognet, method="random")

## -----------------------------------------------------------------------------
print(nma.alog)

# Draw plot of NMA estimates plotted by dose
plot(nma.alog)

## ---- results="hide", warning=FALSE-------------------------------------------
# Run an Emax dose-response MBNMA
mbnma <- mbnma.run(tripnet, fun=demax(emax="rel", ed50="rel"), 
                   method="random")

## ---- message=FALSE, warning=FALSE--------------------------------------------
# Print neat summary of output
summary(mbnma)

## ---- results="hide"----------------------------------------------------------
# Emax model with single parameter estimated for Emax
emax <- mbnma.run(tripnet, fun=demax(emax="rel", ed50="common"), 
                  method="random")

## -----------------------------------------------------------------------------
summary(emax)

## ---- results="hide", message=FALSE, warning=FALSE----------------------------
# Using the osteoarthritis dataset
pain.df <- osteopain

# Set class equal to agent for all agents
pain.df$class <- pain.df$class

# Set a shared class (NSAID) only for Naproxcinod and Naproxen
pain.df$class[pain.df$agent %in% c("Naproxcinod", "Naproxen")] <- 
  "NSAID"

# Run a restricted cubic spline MBNMA with a random class effect on beta.1
classnet <- mbnma.network(pain.df)
splines <- mbnma.run(classnet, fun=dspline(type="bs", knots=2), class.effect = list(beta.1="random"))

## ---- eval=FALSE--------------------------------------------------------------
#  # Using the depression SSRI dataset
#  depnet <- mbnma.network(ssri)
#  
#  # An example specifying a quadratic dose-response function
#  quadfun <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
#  
#  quad <- mbnma.run(depnet, fun=duser(fun=quadfun, beta.1 = "rel", beta.2 = "rel"))|

## ---- eval=FALSE--------------------------------------------------------------
#  # Using the depression SSRI dataset
#  depnet <- mbnma.network(ssri)
#  
#  dr.funs <- dmulti(list(
#    "Placebo"=dpoly(degree=2),
#    "citalopram"=dpoly(degree=2),
#    "escitalopram"=dpoly(degree=2),
#    "fluoxetine"=dspline(type="ns",knots=2),
#    "paroxetine"=dpoly(degree=2),
#    "sertraline"=dspline(type="ns",knots=2)
#    ))
#  
#  multifun <- mbnma.run(depnet, fun=dr.funs, method="common", n.iter=50000)
#  summary(multifun)

## ---- eval=FALSE--------------------------------------------------------------
#  dspline(type="bs", knots=3)
#  # ...is equivalent to
#  dspline(type="bs", knots=c(0.25,0.5,0.75))
#  
#  # Using a natural cubic spline on the SSRI dataset
#  depnet <- mbnma.network(ssri)
#  ns <- mbnma.run(depnet, fun=dspline(type="ns", knots=c(0.25,0.5,0.75)))

## -----------------------------------------------------------------------------
print(mbnma$model.arg$priors)

## ---- eval=FALSE--------------------------------------------------------------
#  # Define replacement prior
#  new.priors <- list(
#    sd = "dnorm(0, 1) T(0,)"
#    )
#  
#  # Run an MBNMA model with new priors
#  emax <- mbnma.run(alognet, fun=demax(), method="random",
#                     priors=new.priors)

## ---- message=FALSE, warning=FALSE--------------------------------------------
# Generate dataset without placebo
noplac.gout <- 
  gout[!gout$studyID %in% c(2001, 3102),] # Drop two-arm placebo studies
noplac.gout <- 
  noplac.gout[noplac.gout$agent!="Plac",] # Drop placebo arm from multi-arm studies

# Create mbnma.network object
noplac.net <- mbnma.network(noplac.gout)

## -----------------------------------------------------------------------------
# Plot network
plot(noplac.net, label.distance=5)

## -----------------------------------------------------------------------------
# Network plot at the agent level illustrates how doses can connect using MBNMA
plot(noplac.net, level="agent", remove.loops = TRUE, label.distance = 4)

## -----------------------------------------------------------------------------
# Network plot assuming connectivity via two doses
# Allows estimation of a single-parameter dose-response function
plot(noplac.net, level="agent", remove.loops = TRUE, label.distance = 4,
     doselink=2)

# Network plot assuming connectivity via three doses
# Allows estimation of a two-parameter dose-response function
plot(noplac.net, level="agent", remove.loops = TRUE, label.distance = 4,
     doselink=3)

## ---- nonparam, results="hide"------------------------------------------------
nonparam <- mbnma.run(tripnet, fun=dnonparam(direction="increasing"), method="random")

## -----------------------------------------------------------------------------
print(nonparam)

## ---- results="hide"----------------------------------------------------------
tripnet <- mbnma.network(triptans)
trip.emax <- mbnma.run(tripnet, fun=demax(emax="rel", ed50="rel")) 

## ---- results="hide"----------------------------------------------------------
# Plot boxplots of residual deviance contributions (scatterplot is the default)
devplot(trip.emax, plot.type = "box")

## ---- results="hide", warning=FALSE-------------------------------------------
# Plot fitted and observed values with treatment labels
fitplot(trip.emax)

## ---- results="hide"----------------------------------------------------------
plot(trip.emax)

## -----------------------------------------------------------------------------
# Specify treatments (agents and doses) for which to estimate relative effects
treats <- list("Placebo"=0,
               "eletriptan"= 1,
               "sumatriptan"=2,
               "almotriptan"=1)

# Print relative effects on the natural scale
rels <- get.relative(trip.emax, treatments = treats, eform=TRUE)
print(rels)

# Rank relative effects
rank(rels)

## -----------------------------------------------------------------------------
ranks <- rank(trip.emax, lower_better = FALSE)
print(ranks)
summary(ranks)

## -----------------------------------------------------------------------------
# Ranking histograms for Emax
plot(ranks, params = "emax")

# Ranking histograms for ED50
plot(ranks, params = "ed50")

## -----------------------------------------------------------------------------
# Cumulative ranking plot for both dose-response parameters
cumrank(ranks, sucra=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  E0 <- triptans[triptans$dose==0,]

## -----------------------------------------------------------------------------
# Predict 20 doses for each agent, with a stochastic distribution for E0
doses <- list("Placebo"=0, 
                  "eletriptan"=3,
                  "sumatriptan"=3,
                  "almotriptan"=3,
                  "zolmitriptan"=3,
                  "naratriptan"=3,
                  "rizatriptan"=3)

pred <- predict(trip.emax, E0="rbeta(n, shape1=2, shape2=10)",
                      max.doses=doses, n.dose=20)


# Predict exact doses for two agents, and estimate E0 from the data
E0.data <- triptans[triptans$dose==0,]
doses <- list("eletriptan"=c(0,1,3),
                  "sumatriptan"=c(0,3))

## ---- results="hide"----------------------------------------------------------
pred <- predict(trip.emax, E0=E0.data,
                      exact.doses=doses)

## -----------------------------------------------------------------------------
summary(pred)

## -----------------------------------------------------------------------------
# Predict responses using default doses up to the maximum of each agent in the dataset
pred <- predict(trip.emax, E0=0.2, n.dose=20)

plot(pred)

## -----------------------------------------------------------------------------
plot(pred, disp.obs = TRUE)

## ---- results="hide", warning=FALSE, message=FALSE----------------------------
alognet <- mbnma.network(alog_pcfb)
alog.emax <- mbnma.run(alognet, fun=demax(), method="random")
pred <- predict(alog.emax, E0=0, n.dose=20)
plot(pred, overlay.split = TRUE, method="random")

## -----------------------------------------------------------------------------
pred <- predict(trip.emax, E0=0.2, n.doses=4,
                max.doses = list("eletriptan"=5, "sumatriptan"=5, 
                              "frovatriptan"=5, "zolmitriptan"=5))

ranks <- rank(pred)
plot(ranks)

## ---- results="hide", warning=FALSE, message=FALSE----------------------------
# Using the alogliptin dataset
alognet <- mbnma.network(alog_pcfb)
nma <- nma.run(alognet, method="random")
ume <- nma.run(alognet, method="random", UME = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
kable(data.frame(
  "Model"=c("NMA", "UME"),
  "resdev"=c(nma$jagsresult$BUGSoutput$median$totresdev,
                        ume$jagsresult$BUGSoutput$median$totresdev),
  "sd"=c(MBNMAdose:::neatCrI(nma$jagsresult$BUGSoutput$summary[rownames(nma$jagsresult$BUGSoutput$summary)=="sd", c(3,5,7)], digits = 2),
                       MBNMAdose:::neatCrI(ume$jagsresult$BUGSoutput$summary[rownames(ume$jagsresult$BUGSoutput$summary)=="sd", c(3,5,7)], digits=2))
), digits=2,
col.names=c("Model", "Residual Deviance", "Betwen-study SD"))

## ---- results="hide", warning=FALSE, message=FALSE----------------------------
# Compares residual deviance contributions from NMA and UME models
devdev(nma, ume, dev.type="resdev")

## ---- results="hide", warning=FALSE, message=FALSE, fig.show = "hide", eval=FALSE----
#  # Using the psoriasis dataset (>75% improvement in PASI score)
#  psoriasis$r <- psoriasis$r75
#  psorinet <- mbnma.network(psoriasis)
#  
#  # Identify comparisons on which node-splitting is possible
#  splitcomps <- inconsistency.loops(psorinet$data.ab, incldr=TRUE)
#  print(splitcomps)
#  
#  # If we want to fit an Emax dose-response function, there is insufficient
#  #indirect evidence in all but the first 6 comparisons
#  nodesplit <- mbnma.nodesplit(psorinet, fun=demax(), comparisons=splitcomps[1:6,], method="common")

## ---- eval=FALSE--------------------------------------------------------------
#  print(nodesplit)

## ---- eval=FALSE--------------------------------------------------------------
#  # Plot forest plots of direct, indirect and pooled (MBNMA) results for each comparison
#  plot(nodesplit, plot.type="forest")
#  
#  # Plot posterior densities of direct and indirect results for each nodesplit comparisons
#  plot(nodesplit, plot.type="density")

