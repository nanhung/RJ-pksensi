# Required packages -------------------------------------------------------
library(pksensi)
library(httk)
library(deSolve)

# GNU MCSim installation --------------------------------------------------
mcsim_install()

# Example 1: One-compartment pbtk model ------------------------------------
# Construct 1-cpt pbtk model for deSolve package
pbtk1cpt <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dAgutlument = - kgutabs * Agutlument
    dAcompartment = kgutabs * Agutlument - ke * Acompartment
    dAmetabolized = ke * Acompartment
    Ccompartment = Acompartment / vdist * BW;
    list(c(dAgutlument, dAcompartment, dAmetabolized), 
         "Ccompartment" = Ccompartment) 
  })
}

# Extract parameter value from httk package
pars1comp <- (parameterize_1comp(chem.name = "acetaminophen"))

# Define parameters that will be used in 1-cpt pbtk model
parms <- c(vdist = pars1comp$Vdist, 
           ke = pars1comp$kelim, 
           kgutabs = pars1comp$kgutabs, 
           BW = pars1comp$BW)

# Define initial condition of each variable
initState <- c(Agutlument = 10, Acompartment = 0, Ametabolized = 0)

# Set output time steps
t <- seq(from = 0.01, to = 24.01, by = 1)

# Run deSolve::ode
y <- ode(y = initState, times = t, func = pbtk1cpt, parms = parms)

# Figure 2 ----------------------------------------------------------------
par(mar=c(4,2,2,1))
plot(y)

# GSA for 1-cpt pbtk model -------------------------------------------------
# Set parameter range 
q <- c("qunif", "qunif", "qunif", "qnorm")
q.arg <- list(list(min = pars1comp$Vdist / 2, max = pars1comp$Vdist * 2),
              list(min = pars1comp$kelim / 2, max = pars1comp$kelim * 2),
              list(min = pars1comp$kgutabs / 2, max = pars1comp$kgutabs * 2),
              list(mean = pars1comp$BW, sd = 5))
params <- c("vdist", "ke", "kgutabs", "BW")


# Parameter matrix generation
set.seed(1234)
x <- rfast99(params, n = 200, q = q, q.arg = q.arg, replicate = 10)

# Review data structure
dim(x$a)

# Figure 3 ----------------------------------------------------------------
par(mfrow=c(4,4),mar=c(0.8,0.8,0.8,0),oma=c(4,4,2,1), pch =".")
for (j in c("vdist", "ke", "kgutabs", "BW")) {
  if ( j == "BW") {
    plot(x$a[,1,j], ylab = "BW")
  } else plot(x$a[,1,j], xaxt="n", ylab = "")
  for (i in 2:3) {
    if ( j == "BW") {
      plot(x$a[,i,j], ylab = "", yaxt="n")  
    } else plot(x$a[,i,j], xaxt="n", yaxt="n", ylab = "")
  } 
  hist <- hist(x$a[,,j], plot=FALSE, 
               breaks=seq(from=min(x$a[,,j]), to=max(x$a[,,j]), length.out=20))
  barplot(hist$density, axes=FALSE, space=0, horiz = T, main = j) 
}
mtext("Model evaluation", SOUTH<-1, line=2, outer=TRUE)

# Figure 4 ----------------------------------------------------------------
outputs <- c("Ccompartment", "Ametabolized")
out <- solve_fun(x, time = t, func = pbtk1cpt, 
                 initState = initState, outnames = outputs)

# Figure 5 ----------------------------------------------------------------
plot(out)
plot(out, vars = "Ametabolized")

# Figure 6 ----------------------------------------------------------------
par(mfcol=c(4,4),mar=c(0.8,0.8,0,0),oma=c(4,4,2,1), pch = ".")
plot(x$a[,1,"vdist"], out$y[,1,"0.01",1], xaxt="n", main = "\nvdist")
plot(x$a[,1,"vdist"], out$y[,1,"2.01",1], xaxt="n")
plot(x$a[,1,"vdist"], out$y[,1,"6.01",1], xaxt="n")
plot(x$a[,1,"vdist"], out$y[,1,"24.01",1])
for (j in c("ke", "kgutabs", "BW")){
  for (k in c("0.01", "2.01", "6.01", "24.01")){
    if (k == "0.01") {
      plot(x$a[,1,j], out$y[,1,k,1], yaxt = "n", xaxt="n", main = paste0("\n", j))
    } else if (k == "24.01") {
      plot(x$a[,1,j], out$y[,1,k,1], yaxt = "n")
    } else plot(x$a[,1,j], out$y[,1,k,1], xaxt = "n", yaxt = "n")
  }
}
mtext("Parameter", SOUTH<-1, line=2, outer=TRUE)
mtext("Ccompartment", WEST<-2, line=2, outer=TRUE)

# Review data structure
dim(out$y)

# Check result in console
check(out, SI.cutoff = 0.05)


# Model implementation with GNU MCSim ------------------------------------
# Benchmark time from deSolve
system.time(out <- solve_fun(x, time = t, 
                             func = pbtk1cpt, initState = initState, 
                             outnames = outputs))

# Benchmark time from GNU MCSim
pbtk1cpt_model()
mName <- "pbtk1cpt"
compile_model(mName, application = "mcsim")
conditions <- c("Agutlument = 10") 
system.time(out <- solve_mcsim(x, mName = mName, params = params, 
                               vars = outputs, time = t, 
                               condition = conditions))

# Example 2: Acetaminophen-PBPK model ------------------------------------
# Download and compile model
mName <- "pbpk_apap"
pbpk_apap_model()
compile_model(mName, application = "mcsim")

# GSA for Acetaminophen-PBPK model ---------------------------------------
# Set parameter info and range 
Tg <- log(0.23)
Tp <- log(0.033)
CYP_Km <- log(130)
SULT_Km_apap <- log(300)
SULT_Ki <- log(526)
SULT_Km_paps <- log(0.5)
UGT_Km <- log(6.0e3)
UGT_Ki <- log(5.8e4)
UGT_Km_GA <-log(0.5)
Km_AG <- log(1.99e4)
Km_AS <- log(2.29e4)
rng <- 1.96 

params <- c("lnTg", "lnTp", "lnCYP_Km","lnCYP_VmaxC",
            "lnSULT_Km_apap","lnSULT_Ki","lnSULT_Km_paps","lnSULT_VmaxC",
            "lnUGT_Km","lnUGT_Ki","lnUGT_Km_GA","lnUGT_VmaxC",
            "lnKm_AG","lnVmax_AG","lnKm_AS","lnVmax_AS",
            "lnkGA_syn","lnkPAPS_syn", "lnCLC_APAP","lnCLC_AG","lnCLC_AS")
dist <- rep("Uniform", 21)
q <- rep("qunif", 21)
q.arg <-list(list(Tg-rng, Tg+rng), list(Tp-rng, Tp+rng), 
             list(CYP_Km-rng, CYP_Km+rng), list(-2., 5.),
             list(SULT_Km_apap-rng, SULT_Km_apap+rng),
             list(SULT_Ki-rng, SULT_Ki+rng),
             list(SULT_Km_paps-rng, SULT_Km_paps+rng),
             list(0, 10), list(UGT_Km-rng, UGT_Km+rng),
             list(UGT_Ki-rng, UGT_Ki+rng),
             list(UGT_Km_GA-rng, UGT_Km_GA+rng),
             list(0, 10), list(Km_AG-rng, Km_AG+rng),
             list(7., 15), list(Km_AS-rng, Km_AS+rng),
             list(7., 15), list(0., 13), list(0., 13),
             list(-6., 1), list(-6., 1), list(-6., 1))

# Define the exposure condition
conditions <- c("mgkg_flag = 1",
                "OralExp_APAP = NDoses(2, 1, 0, 0, 0.001)",
                "OralDose_APAP_mgkg = 20.0")

# Set the examined variables
vars <- c("lnCPL_APAP_mcgL", "lnCPL_AG_mcgL", "lnCPL_AS_mcgL")

# Set output time steps
times <- seq(from = 0.1, to = 12.1, by = 0.2)

set.seed(1111)
out <- solve_mcsim(mName = mName, params = params, vars = vars,
                   monte_carlo = 1000, dist = dist, q.arg = q.arg, 
                   time = times, condition = conditions, 
                   rtol = 1e-7, atol = 1e-9)

# Acetaminophen PK data
head(APAP)

# Figure 7 ----------------------------------------------------------------
par(mfrow = c(1,3), mar = c(4,4,1,1))
pksim(out, xlab = "Time (h)", ylab = "Conc. (ug/L)", main = "APAP")
points(APAP$Time, log(APAP$APAP * 1000))
pksim(out, vars = "lnCPL_AG_mcgL", xlab = "Time (h)", main = "APAP-G", 
      ylab = " ", legend = FALSE)
points(APAP$Time, log(APAP$AG * 1000))
pksim(out, vars = "lnCPL_AS_mcgL", xlab = "Time (h)", main = "APAP-S", 
      ylab = " ", legend = FALSE)
points(APAP$Time, log(APAP$AS * 1000))


# Generate parameter matrix
set.seed(1234)
x <- rfast99(params = params, n = 512, q = q, q.arg = q.arg, replicate = 10) 

# Solve ODE through GNU MCSim
out <- solve_mcsim(x, mName = mName,
                   params = params, 
                   time = times, 
                   vars = vars,
                   condition = conditions, 
                   rtol = 1e-7, atol = 1e-9)

# Figure 8 ----------------------------------------------------------------
plot(out, vars = "lnCPL_APAP_mcgL")

# Figure 9 ----------------------------------------------------------------
heat_check(out, order = "total order")

# Figure 10 ----------------------------------------------------------------
heat_check(out, index = "CI", order = "total order")
