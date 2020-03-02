## ----setup, include = FALSE----------------------------------------------
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

## ---- echo=FALSE---------------------------------------------------------
kable(head(HF2PPITT), digits=2) 

## ---- echo=FALSE---------------------------------------------------------
kable(head(GoutSUA_2wkCFB), digits=2) 

## ---- echo=FALSE---------------------------------------------------------
kable(head(osteopain_2wkabs), digits=2) 

## ---- echo=FALSE---------------------------------------------------------
kable(head(alog_pcfb), digits=2) 

## ------------------------------------------------------------------------
# Using the triptans dataset
network <- mbnma.network(HF2PPITT)
print(network)

## ---- message=FALSE, warning=FALSE---------------------------------------
# Prepare data using the gout dataset
goutnet <- mbnma.network(GoutSUA_2wkCFB)

## ------------------------------------------------------------------------
plot(goutnet, label.distance = 5)

## ------------------------------------------------------------------------
# Plot at the agent-level
plot(goutnet, level="agent", label.distance = 6)

## ------------------------------------------------------------------------
# Plot connections to placebo via a two-parameter dose-response function (e.g. Emax)
plot(goutnet, level="agent", doselink = 2, remove.loops = TRUE, label.distance = 6)

## ---- results="hide"-----------------------------------------------------
# Colour vertices by agent
plot(goutnet, v.color = "agent", label.distance = 5)

## ---- results="hide", message=FALSE, warning=FALSE-----------------------
# Run a random effect split NMA using the alogliptin dataset
alognet <- mbnma.network(alog_pcfb)
nma.alog <- nma.run(alognet, method="random")

## ------------------------------------------------------------------------
print(nma.alog)

# Draw plot of NMA estimates plotted by dose
plot(nma.alog)

## ---- results="hide", warning=FALSE--------------------------------------
# Run an Emax dose-response MBNMA
mbnma <- mbnma.run(tripnet, fun="emax", 
                   beta.1="rel", beta.2="rel", method="common")

## ---- message=FALSE, warning=FALSE---------------------------------------
summary(mbnma)

## ---- eval=FALSE---------------------------------------------------------
#  # An alternative would be to use an Emax wrapper for mbnma.run() which would give the
#  #same result but with more easily interpretable parameter names
#  mbnma.emax(tripnet, emax="rel", ed50="rel", method="common")

## ---- results="hide"-----------------------------------------------------
# Emax model with single parameter estimated for Emax
emax <- mbnma.emax(tripnet, emax="rel", ed50="common", method="random")

## ------------------------------------------------------------------------
summary(emax)

## ---- results="hide", message=FALSE, warning=FALSE-----------------------
# Using the osteoarthritis dataset (contains info on classes)
painnet <- mbnma.network(osteopain_2wkabs)

## ---- results="hide"-----------------------------------------------------
# Emax MBNMA with a single absolute parameter for Emax and random relative effects on ED50
# Common class effects on ED50
emax <- mbnma.emax(painnet, emax="common", ed50="rel", method="random",
                   class.effect=list(ed50="common"))

## ---- eval=FALSE---------------------------------------------------------
#  summary(emax)

## ---- eval=FALSE---------------------------------------------------------
#  # Define relative magnitudes of Emax and ED50
#  rel.size <- c(4,1)
#  
#  mbnma.emax(tripnet, emax="rel", ed50="rel", method="random",
#             var.scale=rel.size)

## ---- eval=FALSE---------------------------------------------------------
#  # An example specifying a quadratic dose-response function
#  quad <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
#  
#  mbnma.run(alognet, fun="user", user.fun=quad,
#            beta.1 = "rel", beta.2 = "rel")

## ---- eval=FALSE---------------------------------------------------------
#  # Frovatriptan and rizatriptan modelled using an exponential function as only two active
#  #doses are available for these agents
#  # All other agents are modelled using an Emax function
#  multifun <- mbnma.run(painnet, fun=c(rep("emax", 3), "exponential", rep("emax", 2), "exponential", "emax"),
#              beta.1="rel", beta.2="rel", beta.3="random")
#  summary(multifun)

## ------------------------------------------------------------------------
print(mbnma$model.arg$priors)

## ---- eval=FALSE---------------------------------------------------------
#  # Define replacement prior
#  new.priors <- list(
#    sd = "dnorm(0, 1) T(0,)"
#    )
#  
#  # Run an MBNMA model with new priors
#  emax <- mbnma.emax(alognet, emax="rel", ed50="rel", method="random",
#                     priors=new.priors)

## ---- result="hide", message=FALSE, warning=FALSE, error=TRUE, eval=FALSE----
#  # Exponential model currently returns an error on gout dataset
#  exponetial <- mbnma.exponential(goutnet, lambda="rel")

## ---- eval=FALSE---------------------------------------------------------
#  # Standardise gout dataset doses relative to maximum investigated dose
#  std.gout <- GoutSUA_2wkCFB
#  std.gout <- std.gout %>% group_by(agent) %>% mutate(dose=dose/max(dose))
#  std.gout$dose[std.gout$agent=="Plac"] <- 0
#  std.gout <- mbnma.network(std.gout)
#  
#  exponential <- mbnma.exponential(std.gout, lambda="rel")

## ---- message=FALSE, warning=FALSE---------------------------------------
# Generate dataset without placebo
noplac.gout <- 
  GoutSUA_2wkCFB[!GoutSUA_2wkCFB$studyID %in% c(2001, 3102),] # Drop two-arm placebo studies
noplac.gout <- 
  noplac.gout[noplac.gout$agent!="Plac",] # Drop placebo arm from multi-arm studies

# Create mbnma.network object
noplac.net <- mbnma.network(noplac.gout)

## ------------------------------------------------------------------------
# Plot network
plot(noplac.net, label.distance=5)

## ------------------------------------------------------------------------
# Network plot at the agent level illustrates how doses can connect using MBNMA
plot(noplac.net, level="agent", remove.loops = TRUE, label.distance = 4)

## ------------------------------------------------------------------------
# Network plot assuming connectivity via two doses
# Allows estimation of a single-parameter dose-response function
plot(noplac.net, level="agent", remove.loops = TRUE, label.distance = 4,
     doselink=2)

# Network plot assuming connectivity via three doses
# Allows estimation of a two-parameter dose-response function
plot(noplac.net, level="agent", remove.loops = TRUE, label.distance = 4,
     doselink=3)

## ---- results="hide"-----------------------------------------------------
nonparam <- mbnma.run(tripnet, fun="nonparam.up", method="random")

## ------------------------------------------------------------------------
print(nonparam)

## ---- results="hide"-----------------------------------------------------
tripnet <- mbnma.network(HF2PPITT)
trip.emax <- mbnma.emax(tripnet, emax="rel", ed50="rel") 

## ---- results="hide"-----------------------------------------------------
# Plot boxplots of residual deviance contributions (scatterplot is the default)
devplot(trip.emax, plot.type = "box")

## ---- results="hide", warning=FALSE--------------------------------------
# Plot fitted and observed values with treatment labels
fitplot(trip.emax)

## ---- results="hide"-----------------------------------------------------
plot(trip.emax)

## ------------------------------------------------------------------------
ranks <- rank(trip.emax, direction = 1)
summary(ranks)

## ------------------------------------------------------------------------
# Ranking histograms for Emax
plot(ranks, params = "d.emax")

# Ranking histograms for ED50
plot(ranks, params = "d.ed50")

## ---- eval=FALSE---------------------------------------------------------
#  E0 <- HF2PPITT[HF2PPITT$dose==0,]

## ------------------------------------------------------------------------
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
E0.data <- HF2PPITT[HF2PPITT$dose==0,]
doses <- list("eletriptan"=c(0,1,3),
                  "sumatriptan"=c(0,3))

## ---- results="hide"-----------------------------------------------------
pred <- predict(trip.emax, E0=E0.data,
                      exact.doses=doses)

## ------------------------------------------------------------------------
summary(pred)

## ------------------------------------------------------------------------
# Predict responses using default doses up to the maximum of each agent in the dataset
pred <- predict(trip.emax, E0=0.2, n.dose=20)

plot(pred)

## ------------------------------------------------------------------------
plot(pred, disp.obs = TRUE)

## ---- results="hide", warning=FALSE, message=FALSE-----------------------
alognet <- mbnma.network(alog_pcfb)
alog.emax <- mbnma.emax(alognet, method="random")
pred <- predict(alog.emax, E0=0, n.dose=20)
plot(pred, overlay.split = TRUE, method="random")

## ------------------------------------------------------------------------
pred <- predict(trip.emax, E0=0.2, n.doses=4,
                max.doses = list("eletriptan"=5, "sumatriptan"=5, 
                              "frovatriptan"=5, "zolmitriptan"=5))

ranks <- rank(pred)
plot(ranks)

## ---- results="hide", warning=FALSE, message=FALSE-----------------------
# Using the alogliptin dataset
alognet <- mbnma.network(alog_pcfb)
nma <- nma.run(alognet, method="random")
ume <- nma.run(alognet, method="random", UME = TRUE)

## ---- echo=FALSE---------------------------------------------------------
kable(data.frame(
  "Model"=c("NMA", "UME"),
  "resdev"=c(nma$jagsresult$BUGSoutput$median$totresdev,
                        ume$jagsresult$BUGSoutput$median$totresdev),
  "sd"=c(MBNMAdose:::neatCrI(nma$jagsresult$BUGSoutput$summary[rownames(nma$jagsresult$BUGSoutput$summary)=="sd", c(3,5,7)], digits = 2),
                       MBNMAdose:::neatCrI(ume$jagsresult$BUGSoutput$summary[rownames(ume$jagsresult$BUGSoutput$summary)=="sd", c(3,5,7)], digits=2))
), digits=2,
col.names=c("Model", "Residual Deviance", "Betwen-study SD"))

## ---- results="hide", warning=FALSE, message=FALSE, fig.show = "hide", eval=FALSE----
#  # Using the triptans dataset
#  tripnet <- mbnma.network(HF2PPITT)
#  
#  # Nodesplit using a random effects model
#  nodesplit <- nma.nodesplit(tripnet, method="random")

## ---- eval=FALSE---------------------------------------------------------
#  print(nodesplit)

## ---- eval=FALSE---------------------------------------------------------
#  # Plot forest plots of direct, indirect and pooled (NMA) results for each comparison
#  plot(nodesplit, plot.type="forest")
#  
#  # Plot posterior densities of direct and indirect results for each nodesplit comparisons
#  plot(nodesplit, plot.type="density")

