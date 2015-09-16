###############################################################################
sum.grid <- function(x,rn,grid.size) {
  y <- NULL
  for (i in 1:dim(x)[2]) y <- cbind(y,rowSums(matrix(x[,i],ncol=grid.size)))
  dimnames(y) <- list(rn,colnames(x))  
  y
}

###############################################################################
expand <- function(x,grid.size) {
  matrix(rep(t(x),grid.size),ncol=dim(x)[2],byrow=T,dimnames=list(rep(rownames(x),grid.size),colnames(x)))
}

###################################### Plot grid concentrations; grid is a vector with x changing fastest
plot.grid <- function(grid,x,zlim=NULL,title="") {
  grid.mat <- matrix(grid,byrow=T,ncol=x)
  if(!is.null(zlim)) zlim<-c(0,zlim)
  myImagePlot(grid.mat,zlim=zlim,title=title)
}

# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  # Data Map
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  #rect(0.5,0.5,length(xLabels)+0.5,length(yLabels)+0.5)
}
############################################################################

MakePlots <- function(out) {
  par(cex.axis=2,cex.lab=2.5,cex.main=3,font.lab=2,font.axis=2,lwd=2,mar=c(5,7,4,2)+0.1,tcl=0.4,las=1,mgp=c(3,0.75,0),mfrow=c(2,5))
  
  RGB.palette <- colorRampPalette(c("red","yellow","blue"),space="rgb")
  
  n_taxa <- dim(out$MicrobesSeries)[2]
  colors.RGB <- RGB.palette(out$params["Enz_per_taxon_max",]+1)
  plot(0:out$end_time,out$Mic_Sum[,"C"],ylim=c(0,max(out$MicrobesSeries/out$grid.size,na.rm=T)),ylab="",main="Microbial C",xlab="Day",type="n",axes=F,frame.plot=T)
  for (i in 1:n_taxa){
    lines(0:out$end_time,out$MicrobesSeries[,i]/out$grid.size,col=colors.RGB[(rowSums(out$EnzGenes)+1)[i]],lwd=2)
  }
  axis(1,lwd=2); axis(2,lwd=2)
  
  n_subs <- length(dimnames(out$SubstratesSeries)[2][[1]])
  color.sub <- c("brown","darkblue","goldenrod","salmon","orange","pink","gray","purple","blue","turquoise","forestgreen","green")
  plot(0:out$end_time,out$Substrate_Sum,ylim=c(0,max(out$SubstratesSeries/out$grid.size,na.rm=T)),main=expression(bold(paste("Substrate C (mg ",cm^-3,")"))),ylab=NA,xlab="Day",type="n",axes=F,frame.plot=T)
  axis(1,lwd=0,lwd.ticks=2); axis(2,lwd=0,lwd.ticks=2)
  for (i in 1:n_subs){
    lines(0:out$end_time,out$SubstratesSeries[,i]/out$grid.size,col=color.sub[i],lwd=2)
  }
  legend("topright",legend=c("Dead microbe","Inactive enzyme","Cellulose","Hemicellulose","Starch","Chitin","Lignin","Protein 1","Protein 2","Protein 3","Phospholipid","Nucleic acid"),lty=1,col=color.sub,box.lty=0,cex=1.1)

  plot(0:out$end_time,out$Mic_Sum[,"C"]/out$Mic_Sum[,"P"],ylim=c(0,max(out$MicrobesPSeries,na.rm=T)),ylab=NA,main="Microbial C:P",xlab="Day",type="l",axes=F,frame.plot=T)
  for (i in 1:n_taxa){
    lines(0:out$end_time,out$MicrobesPSeries[,i],col=colors.RGB[(rowSums(out$EnzGenes)+1)[i]])
  }
  axis(1,lwd=2); axis(2,lwd=2)
  
  plot(0:out$end_time,out$Mic_Sum[,"C"]/out$Mic_Sum[,"N"],ylim=c(0,max(out$MicrobesNSeries,na.rm=T)),ylab=NA,main="Microbial C:N",xlab="Day",type="l",axes=F,frame.plot=T)
  for (i in 1:n_taxa){
    lines(0:out$end_time,out$MicrobesNSeries[,i],col=colors.RGB[(rowSums(out$EnzGenes)+1)[i]])
  }
  axis(1,lwd=2); axis(2,lwd=2)
  
  plot(0:out$end_time,out$Mic_Sum[,"N"]/out$Mic_Sum[,"P"],ylim=c(0,max(out$MicrobesNPSeries,na.rm=T)),ylab=NA,main="Microbial N:P",xlab="Day",type="l",axes=F,frame.plot=T)
  for (i in 1:n_taxa){
    lines(0:out$end_time,out$MicrobesNPSeries[,i],col=colors.RGB[(rowSums(out$EnzGenes)+1)[i]])
  }
  axis(1,lwd=2); axis(2,lwd=2)
  
  # Color code enzymes according to the main substrate they degrade
  enz.potential <- out$ReqEnz[[1]]*out$Vmax
  colors <- color.sub[apply(enz.potential,2,which.max)]
  plot(0:out$end_time,out$Enzyme_Sum,ylim=c(0,max(out$EnzymesSeries/out$grid.size,na.rm=T)),ylab=NA,main="Enzyme",xlab="Day",type="n",axes=F,frame.plot=T)
  for (i in 1:length(colors)){
    lines(0:out$end_time,out$EnzymesSeries[,i]/out$grid.size,col=colors[i],lwd=2)
  }
  axis(1,lwd=2); axis(2,lwd=2)
  
  plot(0:out$end_time,out$NH4Series,ylim=c(0,max(cbind(out$NH4Series/out$grid.size,out$PO4Series/out$grid.size))),
       ylab=NA,main="NH4 | PO4",xlab="Day",type="n",axes=F,frame.plot=T)
  lines(0:out$end_time,out$NH4Series/out$grid.size,col="orange")
  lines(0:out$end_time,out$PO4Series/out$grid.size,col="purple")
  axis(1,lwd=2); axis(2,lwd=2)
  
  n_monomers <- length(dimnames(out$MonomersSeries)[2][[1]])
  colors <- c("black","black",color.sub)
  plot(0:out$end_time,out$Monomer_Sum,ylim=c(0,max(out$MonomersSeries/out$grid.size,na.rm=T)),ylab=NA,main="Monomer C",xlab="Day",type="n",axes=F,frame.plot=T)
  for (i in 1:n_monomers){
    lines(0:out$end_time,out$MonomersSeries[,i]/out$grid.size,col=colors[i])
  }
  axis(1,lwd=2); axis(2,lwd=2)
  
  plot(0:out$end_time,out$RespSeries/out$grid.size,ylim=c(0,max(out$RespSeries/out$grid.size)),ylab=NA,main="Respiration",xlab="Day",type="l",axes=F,frame.plot=T)
  axis(1,lwd=2); axis(2,lwd=2)
    
  index <- is.finite(log10(apply(out$MicrobesSeries,2,max)))
  plot(rowSums(out$EnzGenes)[index],log10(apply(out$MicrobesSeries,2,max)/out$grid.size)[index],main=expression(bold(paste(Log[10]," max density"))),
       ylab=NA,xlab="Number of enzyme genes",pch=19,cex=2,axes=F,frame.plot=T,#ylim=c(-2,1),
       col=colors.RGB[1+as.vector(rowSums(out$EnzGenes)[index])])
  axis(1,lwd=2);  axis(2,lwd=2)
  abline(lm(log10(apply(out$MicrobesSeries,2,max)/out$grid.size)[index]~rowSums(out$EnzGenes)[index]),lwd=1.75)
  
} # End of MakePlots() ########################################################


###############################################################################
RunPulse <- function(
  params,
  timestamp,
  rng.seed,
  grid.size,
  Microbes,
  Substrates,
  SubInput,
  Enzymes,
  Monomers,
  MonInput,
  MonomersProduced,
  ReprodNew,
  Colonization.reset,
  Ea,
  Vmax0,
  Km0,
  ReqEnz,
  EnzGenes,
  EnzProdInduce,
  EnzProdConstit,
  UptakeGenes,
  UptakeGenesForEnz,
  Uptake_ReqEnz,
  EnzAttrib,
  Uptake_Ea,
  Uptake_Vmax0,
  Uptake_Km0,
  CUE.ref,
  OptimalRatios,
  MinRatios,
  RangeRatios,
  fb,
  Temp,
  Psi,
  Tolerance
) {
  # Initialize names and indicies
  n_taxa <- params["n_taxa",]
  n_substrates <- params["n_substrates",]
  n_enzymes <- params["n_enzymes",]
  Mon.names <- c("NH4","PO4",rownames(Substrates[1:n_substrates,]))
  n_monomers <- length(Mon.names)
  Mon.names.G <- rownames(Monomers)
  Enz.names <- colnames(Vmax0)
  Sub.names <- rownames(Substrates[1:n_substrates,])
  Mic.names <- rownames(Microbes[1:n_taxa,])
  is.NH4 <- which(rownames(Monomers)=="NH4")
  is.PO4 <- which(rownames(Monomers)=="PO4")
  org <- which(rownames(Monomers) != c("NH4","PO4"))
  mineral <- which(rownames(Monomers) == c("NH4","PO4"))
  is.deadMic <- which(rownames(Substrates)=="DeadMic")
  is.deadEnz <- which(rownames(Substrates)=="DeadEnz")
  is.cellulose <- which(rownames(Substrates)=="Cellulose")
  is.lignin <- which(rownames(Substrates)=="Lignin")
  MonomerRatios <- Monomers
  MonomerRatios[,] <- 0
  MonomerRatios[is.NH4,"N"] <- 1
  MonomerRatios[is.PO4,"P"] <- 1
  Tref <- 293
  # Create matrices to be filled in
  tev <- matrix(rep(0,n_substrates*grid.size*n_enzymes),ncol=n_enzymes)
  Max_Uptake1 <- matrix(rep(0,grid.size*n_taxa*n_monomers),nrow=n_taxa)
  Max_Uptake2 <- matrix(Max_Uptake1,ncol=n_taxa,dimnames=list(Mon.names.G,Mic.names))
  TUC.mat <- TUN.mat <- TUP.mat <- matrix(rep(0,grid.size*n_taxa*n_monomers),nrow=n_monomers)
  EP.mat <- matrix(rep(0,grid.size*n_taxa*n_enzymes),nrow=n_taxa,dimnames=list(Mic.names,rownames(Enzymes)))
  Death.mat <- matrix(rep(0,grid.size*n_taxa*3),nrow=n_taxa)
  DeadEnz.mat <- matrix(rep(0,grid.size*n_enzymes*3),nrow=n_enzymes)
  # Create index for matrix manipulation of taxon uptake
  # index1 is the shape of the starting matrix
  TU.index1 <- matrix(1:(grid.size*n_taxa*n_monomers),ncol=n_taxa)
  TU.index <- matrix(TU.index1,nrow=n_monomers)
  for(i in 1:grid.size) {
    rows <- ((i-1)*n_monomers+1):(i*n_monomers)
    cols <- ((i-1)*n_taxa+1):(i*n_taxa)
    TU.index[,cols] <- TU.index1[rows,]
  }
  # Create index for matrix manipulation of taxon enzyme production
  EP.index1 <- matrix(1:(grid.size*n_taxa*n_enzymes),ncol=n_enzymes)
  EP.index <- matrix(EP.index1,nrow=n_taxa)
  for(i in 1:grid.size) {
    rows <- ((i-1)*n_taxa+1):(i*n_taxa)
    cols <- ((i-1)*n_enzymes+1):(i*n_enzymes)
    EP.index[,cols] <- EP.index1[rows,]
  }
  # Index to manipulate Vmax values
  vm.index1 <- matrix(1:(grid.size*n_substrates*n_enzymes),ncol=n_substrates)
  vm.index <- matrix(vm.index1,ncol=n_enzymes)
  for(i in 1:grid.size) {
    vm.index[((i-1)*n_substrates+1):(i*n_substrates),] <- t(vm.index1[((i-1)*n_enzymes+1):(i*n_enzymes),])
  }
  # Set grid sums to initial values
  Monomers.grid <- sum.grid(Monomers,Mon.names,grid.size)
  Enzymes.grid <- sum.grid(Enzymes,Enz.names,grid.size)
  Substrates.grid <- sum.grid(Substrates,Sub.names,grid.size)
  Microbes.grid <- sum.grid(Microbes,Mic.names,grid.size)
  n_fungi <- sum(fb)/grid.size
  Fungi <- which(fb[1:n_taxa]==1)
  is.fungi <- fb==1
  Fungi.count.zero <- Fungi.count <- rep(0,n_taxa)
  Fungi.vec <- Fungi.vec.zero <- matrix(rep(0,n_taxa*grid.size),ncol=1,dimnames=list(rownames(Microbes),"Count"))
  Fungi.count[Fungi] <- Microbes.grid[Fungi,"C"]/(0.5*params["max_size_f",])
  z <- rep(1:n_taxa,grid.size)
  kill.reset <- rep(0,grid.size*n_taxa)
  
  # Initialize the time series of data to hold
  RespSeries <- 0
  Cum_Leaching_N <- 0
  Cum_Leaching_P <- 0
  EnzymesSeries <- t(Enzymes.grid[,"C"])
  SubstratesSeries <- Substrates.grid[,"C"]
  MonomersSeries <- Monomers.grid[,"C"]
  NH4Series <- Monomers.grid["NH4","N"]
  PO4Series <- Monomers.grid["PO4","P"]
  MicrobesSeries <- Microbes.grid[,"C"]
  MicrobesNSeries <- Microbes.grid[,"C"]/Microbes.grid[,"N"]
  MicrobesPSeries <- Microbes.grid[,"C"]/Microbes.grid[,"P"]
  MicrobesNPSeries <- Microbes.grid[,"N"]/Microbes.grid[,"P"]
  
  # Sum microbe, substrate, monomer, enzyme pools
  Mic_Sum <- Cum_Microbe <- colSums(Microbes.grid)
  Substrate_Sum <- sum(Substrates.grid[,"C"])
  Cum_Substrate <- colSums(Substrates.grid)
  Monomer_Sum <- sum(Monomers.grid[,"C"])
  Cum_Monomer <- colSums(Monomers.grid)
  Enzyme_Sum <- Cum_Enzyme_C <- sum(Enzymes.grid)
  Cum_Enzyme_N <- sum(Enzymes.grid*EnzAttrib[,"N_cost"])
  Cum_Enzyme_P <- sum(Enzymes.grid*EnzAttrib[,"P_cost"])
  
  # Relative enzyme carbon cost for each enzyme gene
  Induce_Enzyme_C <- t(EnzProdInduce)*EnzAttrib[,"C_cost"]
  Constit_Enzyme_C <- t(EnzProdConstit)*EnzAttrib[,"C_cost"]
  # Relative uptake enzyme carbon cost for each enzyme gene
  Uptake_Cost <- t(UptakeGenes)
  rownames(Uptake_Cost) <- sprintf("%s%03d","Upt",1:params["n_uptake",])
  
  # Expand variables to the size of the grid
  Induce_Enzyme_C <- expand(t(Induce_Enzyme_C),grid.size)
  Constit_Enzyme_C <- expand(t(Constit_Enzyme_C),grid.size)
  Uptake_Cost <- expand(t(Uptake_Cost),grid.size)
  tUptakeGenes <- expand(t(UptakeGenes),grid.size)
  UptakeGenes <- expand(UptakeGenes,grid.size)
  
  # Create variable to hold grid at maximum enzyme biomass
  grid.max <- NULL

  for(i_t in 1:params["end_time",]) { # Loop through time
    
    # Mortality rate increases linearly with -Psi times slope beta
    # which can differ for bacteria vs fungi; tolerance of 1 eliminates
    # dependence on Psi. Fungi and bacteria can also differ in baseline
    # mortality according to Death_Ratio
    r.death <- (1-fb)*params["Death_Rate",]*(1 - (params["beta.bac",]*Psi[i_t])*(1 - Tolerance)) +
      fb*params["Death_Rate",]*((1 - (params["beta.fungi",]*Psi[i_t])*(1 - Tolerance))*params["Death_Ratio",])
   
    # Calculate Vmax values with Arrhenius equation
    # Gas constant = 0.008314 kJ/(mol K)
    # Arrhenius equation for Vmax multiplied by exponential decay for Psi sensitivity
    Vmax <- Vmax0*exp((-Ea/0.008314)*(1/(Temp[i_t]+273)-1/Tref))*exp(params["Psi.slope.Vmax",]*Psi[i_t])
    Uptake_Vmax <- Uptake_Vmax0*exp((-Uptake_Ea/0.008314)*(1/(Temp[i_t]+273)-1/Tref))*exp(params["Psi.slope.uptake",]*Psi[i_t])
    # Expand to the size of the grid after calculating
    etVmax <- expand(t(Vmax),grid.size)
    Uptake_Vmax <- expand(Uptake_Vmax,grid.size)

    # Recalculate Km values with Arrhenius equation
    Km <- Km0*exp((-params["Km_Ea",]/0.008314)*(1/(Temp[i_t]+273)-1/Tref))
    Uptake_Km <- Uptake_Km0*exp((-params["Km_Ea",]/0.008314)*(1/(Temp[i_t]+273)-1/Tref))
    
    # Carbon use efficiency dependent on temperature and number of enzyme genes
    CUE <- CUE.ref + (Temp[i_t] - (Tref-273))*params["CUE_temp",]
    
    # Reset the reproduction matrix
    Reprod <- ReprodNew
    # Reset the colonization matrix
    Colonization <- Colonization.reset
    # Fungal translocation: calculate average biomass within fungal taxa
    Mean.fungi <- Microbes.grid/as.vector(Fungi.count)
    Mean.fungi[!is.finite(Mean.fungi)] <- 0
    # Expand the fungal average across the grid
    eMF <- expand(Mean.fungi,grid.size)
    # Reset the fungal count to zero
    Fungi.count <- Fungi.count.zero
    
    # Select the daughter cells that are fungi versus bacteria
    daughters.b <- which(Reprod[,"C"]>0 & fb==0)
    daughters.f <- which(Reprod[,"C"]>0 & fb==1)
    num.b <- length(daughters.b)
    num.f <- length(daughters.f)
    shift_x <- shift_y <- rep(0,grid.size*n_taxa)
    shift_x[daughters.b] <- sample(c(-params["dist",]:params["dist",]),num.b,replace=T) # vector of dispersal movements in x direction
    shift_y[daughters.b] <- sample(c(-params["dist",]:params["dist",]),num.b,replace=T) # vector of dispersal movements in y direction
    shift_x[daughters.f] <- 1 # Fungi always move positively in x direction
    # vector of dispersal movements in y direction; constrained to one box away determined by probability "direct"
    shift_y[daughters.f] <- sample(c(-1:1),num.f,replace=T,prob=c(0.5*(1-params["direct",]),params["direct",],0.5*(1-params["direct",])))
    new_x <- (rep(1:params["x",],each=n_taxa,times=params["y",]) + shift_x + params["x",])%%params["x",] # calculate x coordinates of dispersal destinations
    new_x[new_x==0] <- params["x",] # Substitute coordinates when there is no shift
    new_y <- (rep(1:params["y",],each=params["x",]*n_taxa) + shift_y + params["y",])%%params["y",] # calculate y coordinates of dispersal destinations
    new_y[new_y==0] <- params["y",]
    Reprod[daughters.f,] <- eMF[daughters.f,] # set all fungi equal to their grid averages for translocation before colonization
    index_col <- n_taxa*((new_y-1)*params["x",]+(new_x-1))+z # convert x,y coordinates to a vector of destination locations
    # Transfer cells to new locations and sum when two or more of the same taxa go to same location
    for(i in c(daughters.b,daughters.f)) {
      Colonization[index_col[i],] <- Colonization[index_col[i],] + Reprod[i,]
    }
    
    # Set Monomers to the grid average
    Monomers <- expand(Monomers.grid/grid.size,grid.size)

    # Translocate nutrients within fungal taxa
    i <- Microbes[fb==1,"C"]>0
    Microbes[fb==1,][i,] <- eMF[fb==1,][i,]
    
    # Colonization of dispersing microbes
    Microbes <- Microbes + Colonization
    
    # Multiply Vmax values for each substrate by the quantity of each enzyme
    tev[,] <- (as.vector(Enzymes)*etVmax)[vm.index]
    
    # Equation for Michaelis-Menten enzyme catalysis
    rss <- rowSums(Substrates)
    Decay <- tev*rss/(Km+rss)
    
    # Loop that pulls out each batch of required enzymes and sums across redundant enzymes
    DecaySums <- NULL
    for(i in 1:dim(ReqEnz)[3]) {
      DecaySums <- cbind(rowSums(as.matrix(ReqEnz[,,i])*Decay),DecaySums)
    }
    
    # Assess the rate-limiting enzyme and set decay to that rate
    # Compare to substrate available and take the min, allowing for a tolerance of 1e-9
    # Link cellulose degradation to lignin concentration (LCI)
    DecayRates <- pmin(pmin(DecaySums[,1],DecaySums[,2],na.rm=T),rss-1e-9*rss,na.rm=T)
    ss7 <- rowSums(Substrates[is.lignin,])
    DecayRates[is.cellulose] <- DecayRates[is.cellulose]*(1+params["LCI_slope",]*ss7/(Substrates[is.cellulose,1]+ss7))
    
    # Monomer stoichiometry: organic monomers follow substrate
    SubstrateRatios <- Substrates/rss
    SubstrateRatios[rss==0,] <- 0
    MonomerRatios[org,] <- SubstrateRatios
    
    # Update monomer pools
    # Decayed substrate stays in the grid box for the whole iteration
    # Preferential access by resident microbes
    # Any leftovers diffuse (get averaged) across grid at next iteration
    Monomers[org,] <- Monomers[org,] + (DecayRates+MonInput[org])*MonomerRatios[org,]
    Monomers[mineral,] <- Monomers[mineral,] + MonInput[mineral]*MonomerRatios[mineral,]
    # Keep track of mass balance for inputs
    Cum_Monomer <- Cum_Monomer + colSums(as.vector(MonInput)*MonomerRatios)
    
    # Recalculate Monomer stoichiometry after changes due to dead microbial biomass
    rsm <- rowSums(Monomers)
    MonomerRatios[org,] <- Monomers[org,]/rsm[org]
    MonomerRatios[org,][rsm[org]==0,] <- 0
    
    # Multiply microbial biomass by each taxon's uptake allocation to get biomass of each uptake enzyme by taxon
    rsi <- rowSums(Microbes)
    MicCXGenes <- rsi*UptakeGenes
    
    # Equation for hypothetical potential uptake (per unit of compatible uptake protein)
    Potential_Uptake <- Uptake_ReqEnz*Uptake_Vmax*rsm/(Uptake_Km+rsm)
    
    # Matrix multiplication to get max possible uptake by taxon and monomer
    # Must extract each grid point separately for operation
    for(i in 1:grid.size) {
      i.a <- ((i-1)*n_monomers+1):(i*n_monomers)
      i.b <- ((i-1)*n_taxa+1):(i*n_taxa)
      Max_Uptake <- MicCXGenes[i.b,]%*%t(Potential_Uptake[i.a,])
      Max_Uptake1[,i.a] <- Max_Uptake # wide format; mon*grid cols
      Max_Uptake2[i.a,] <- t(Max_Uptake) # long format; mon*grid rows
    }
    # Sum the total potential uptake of each monomer
    csmu <- colSums(Max_Uptake1)
    
    # Take the min of the monomer available and the max potential uptake
    Min_Uptake <- pmin(csmu,rsm,na.rm=T)
    
    # Scale the uptake to what's available
    Uptake <- Max_Uptake2*Min_Uptake/csmu
    Uptake[csmu==0,] <- 0
    # Prevent total uptake from getting too close to zero
    Uptake <- Uptake - 1e-9*Uptake
    
    # Calculate total uptake by monomer
    Monomer_Uptake <- rowSums(Uptake)*MonomerRatios
    
    # Calculate total uptake by taxon
    TUC.mat[,] <- (Uptake*MonomerRatios[,"C"])[TU.index]
    TUN.mat[,] <- (Uptake*MonomerRatios[,"N"])[TU.index]
    TUP.mat[,] <- (Uptake*MonomerRatios[,"P"])[TU.index]
    # Sum across monomers
    Taxon_Uptake_C <- colSums(TUC.mat)
    Taxon_Uptake_N <- colSums(TUN.mat)
    Taxon_Uptake_P <- colSums(TUP.mat)
    
    # Enzyme production: rows are taxa, cols are enzyme genes
    # Two alternatives: proportional to biomass or proportional to uptake C
    Taxon_Enzyme_Production <- Taxon_Uptake_C*Induce_Enzyme_C + rsi*Constit_Enzyme_C
    
    # Constrain enzyme production if it will cost more N than currently available in microbial biomass
    Enzyme_Cost_N <- colSums(t(Taxon_Enzyme_Production)*EnzAttrib[,"N_cost"])
    # Only select entries with microbes present
    i <- which(Enzyme_Cost_N>0)
    Taxon_Enzyme_Production[i,] <- Taxon_Enzyme_Production[i,]*pmin(Microbes[i,"N"],Enzyme_Cost_N[i])/Enzyme_Cost_N[i]
    
    # Total enzyme production across all taxa
    EP.mat[,] <- Taxon_Enzyme_Production[EP.index]
    Enzyme_Production <- colSums(EP.mat)
    
    # Total enzyme carbon cost for each taxon
    ttep <- t(Taxon_Enzyme_Production)
    Enzyme_Maint <- colSums(ttep*EnzAttrib[,"Maint_cost"])
    Enzyme_Cost <- colSums(ttep) + Enzyme_Maint
    Enzyme_Cost_N <- colSums(ttep*EnzAttrib[,"N_cost"])
    Enzyme_Cost_P <- colSums(ttep*EnzAttrib[,"P_cost"])
    
    # Uptake enzyme production: cols are taxa, rows are enzymes
    # Total uptake enzyme carbon cost for each taxon
    Uptake_Maint <- colSums(t(rsi*Uptake_Cost))*params["Uptake_Maint_cost",]
    
    # Calculate loss rates and random death
    Enzyme_Loss <- params["Enzyme_Loss_Rate",]*Enzymes
    
    # Update enzyme pools
    Enzymes <- Enzymes + Enzyme_Production - Enzyme_Loss
    
    # Update microbe pools
    Microbes[,"C"] <- Microbes[,"C"] + Taxon_Uptake_C*CUE - Enzyme_Cost - Uptake_Maint
    Microbes[,"N"] <- Microbes[,"N"] + Taxon_Uptake_N - Enzyme_Cost_N
    Microbes[,"P"] <- Microbes[,"P"] + Taxon_Uptake_P - Enzyme_Cost_P
    
    # Kill microbes that are starving and transfer biomass to substrate
    Death <- Colonization.reset
    if (is.na(sum(Microbes))) {break}
    starve_index <- (Microbes[,"C"]<params["C_min",] | Microbes[,"N"]<params["N_min",] | Microbes[,"P"]<params["P_min",]) & Microbes[,"C"]>0
    Death[starve_index,] <- Microbes[starve_index,]
    Microbes[starve_index,] <- 0
    
    # Determine which locations have microbes
    is.m <- which(Microbes[,"C"]>0)
    # Zero the death index
    kill <- kill.reset
    # Sample for death events (=1)
    kill[is.m] <- runif(length(is.m)) < r.death[is.m]
    # Kill microbes
    killed <- kill*Microbes
    Death <- Death + killed
    Microbes <- Microbes - killed
    
    # Adjust stoichiometry of microbial biomass
    # Index locations of microbial cells
    mic.index <- as.numeric(which(Microbes[,"C"]>0))
    # drop=F retains matrix structure if only one row is selected
    Mic.subset <- Microbes[mic.index,,drop=F]
    # Calculate microbial ratios
    MicrobeRatios <- (Mic.subset/rowSums(Mic.subset))
    # Select the corresponding minimum quotas
    MinRat <- MinRatios[mic.index,,drop=F]
    # Index only microbes that have below-minimum quotas
    rat.index <- mic.index[MicrobeRatios[,"C"] < MinRat[,"C"] | MicrobeRatios[,"N"] < MinRat[,"N"] | MicrobeRatios[,"P"] < MinRat[,"P"] ]
    Mic.subset <- Microbes[rat.index,,drop=F]
    StartMicrobes <- Mic.subset
    MicrobeRatios <- Excess <- (Mic.subset/rowSums(Mic.subset))
    MinRat <- MinRatios[rat.index,,drop=F]
    
    # Calculate difference between min and actual ratios    
    Deficit <- Deficit.0 <- MicrobeRatios - MinRat
    # Determine the limiting nutrient that will be conserved
    Limiting <- max.col(-Deficit/MinRat,ties.method="first")
    # Set all deficient ratios to their minima
    MicrobeRatios[Deficit<0] <- MinRat[Deficit<0]
    # Reduce the mass fractions for non-deficient elements in proportion to the distance from the minimum
    # Calculate how far above the minimum each non-deficient element is
    Excess[Deficit>0] <- Deficit[Deficit>0]
    # Set deficient element fractions to zero
    Excess[Deficit<0] <- 0
    # Partition the total deficit to the excess element(s) in proportion to their distances from their minima
    Deficit.0[Deficit>0] <- 0
    MicrobeRatios[Deficit>0] <- MicrobeRatios[Deficit>0] + (rowSums(Deficit.0)*Excess/rowSums(Excess))[Deficit>0]
    
    # Construct hypothetical nutrient quotas for each possible minimum nutrient
    MC <- Mic.subset[,"C"]; MN <- Mic.subset[,"N"]; MP <- Mic.subset[,"P"]
    MRC <- MicrobeRatios[,"C"]; MRN <- MicrobeRatios[,"N"]; MRP <- MicrobeRatios[,"P"]
    new.C <- cbind(MC,MN*MRC/MRN,MP*MRC/MRP)
    new.N <- cbind(MC*MRN/MRC,MN,MP*MRN/MRP)
    new.P <- cbind(MC*MRP/MRC,MN*MRP/MRN,MP)
    # Insert the appropriate set of nutrient quotas scaled to the minimum nutrient
    select <- cbind(1:length(rat.index),Limiting)
    Microbes[rat.index,] <- cbind(new.C[select],new.N[select],new.P[select])
    # Sum up the element losses from biomass across whole grid and calculate average loss
    MicLoss <- StartMicrobes - Microbes[rat.index,]
    
    # Update substrate pools
    # Substrate stoichiometry
    rss <- rowSums(Substrates)
    SubstrateRatios <- Substrates/rss
    SubstrateRatios[rss==0,] <- 0
    # Add inputs and remove decay
    Substrates <- Substrates + SubInput - DecayRates*SubstrateRatios
    # Account for inputs in mass balance
    Cum_Substrate <- Cum_Substrate + colSums(SubInput)
    # Add dead microbial biomass
    Death.mat[,] <- Death
    Substrates[is.deadMic,] <- Substrates[is.deadMic,] + colSums(Death.mat)
    # Add dead enzymes
    DeadEnz.mat[,] <- c(Enzyme_Loss,Enzyme_Loss*EnzAttrib[,"N_cost"],Enzyme_Loss*EnzAttrib[,"P_cost"])
    Substrates[is.deadEnz,] <- Substrates[is.deadEnz,] + colSums(DeadEnz.mat)
    
    # Update monomer pools
    Monomers[is.NH4,"N"] <- Monomers[is.NH4,"N"] + sum(MicLoss[,"N"])/grid.size
    Monomers[is.PO4,"P"] <- Monomers[is.PO4,"P"] + sum(MicLoss[,"P"])/grid.size
    Monomers <- Monomers - Monomer_Uptake
    
    # Sum microbes prior to reproduction
    Microbes.grid <- sum.grid(Microbes,Mic.names,grid.size)
    # Reset the vector of fungal locations
    Fungi.vec <- Fungi.vec.zero
    # Add one or two fungi to the count vector based on size
    Fungi.vec[fb==1][Microbes[fb==1,"C"]>0] <- 1
    Fungi.vec[fb==1][Microbes[fb==1,"C"]>params["max_size_f",]] <- 2
    # Sum up the count vector
    Fungi.count <- sum.grid(Fungi.vec,Mic.names,grid.size)
    
    # Cell division
    MicrobesBeforeDivision <- Microbes
    Microbes[fb==0,][Microbes[fb==0,"C"]>params["max_size_b",],] <- Microbes[fb==0,][Microbes[fb==0,"C"]>params["max_size_b",],]/2
    Microbes[fb==1,][Microbes[fb==1,"C"]>params["max_size_f",],] <- Microbes[fb==1,][Microbes[fb==1,"C"]>params["max_size_f",],]/2
    
    # Add daughter cells to matrix of new reproduction
    ReprodNew <- MicrobesBeforeDivision-Microbes
    
    # Sum pools across grid
    Respiration <- sum(Taxon_Uptake_C*(1-CUE)) + sum(Enzyme_Maint) + sum(Uptake_Maint) + sum(MicLoss[,"C"])
    Monomers.grid <- sum.grid(Monomers,Mon.names,grid.size)
    Enzymes.grid <- sum.grid(Enzymes,Enz.names,grid.size)
    Substrates.grid <- sum.grid(Substrates,Sub.names,grid.size)
    
    # Update monomer leaching
    Leaching_N <- Monomers.grid["NH4","N"]*params["Leaching",]*exp(params["Psi.slope.leach",]*Psi[i_t])
    Leaching_P <- Monomers.grid["PO4","P"]*params["Leaching",]*exp(params["Psi.slope.leach",]*Psi[i_t])
    # Keep track of mass balance
    Cum_Leaching_N <- Cum_Leaching_N + Leaching_N
    Cum_Leaching_P <- Cum_Leaching_P + Leaching_P
    Monomers.grid["NH4","N"] <- Monomers.grid["NH4","N"] - Leaching_N
    Monomers.grid["PO4","P"] <- Monomers.grid["PO4","P"] - Leaching_P
    
    # Record total pools for each time step
    RespSeries <- c(RespSeries,Respiration)
    MonomersSeries <- rbind(MonomersSeries,Monomers.grid[,"C"])
    NH4Series <- rbind(NH4Series,Monomers.grid["NH4","N"])
    PO4Series <- rbind(PO4Series,Monomers.grid["PO4","P"])
    EnzymesSeries <- rbind(EnzymesSeries,t(Enzymes.grid))
    SubstratesSeries <- rbind(SubstratesSeries,Substrates.grid[,"C"])
    MicrobesSeries <- rbind(MicrobesSeries,Microbes.grid[,"C"])
    MicrobesNSeries <- rbind(MicrobesNSeries,Microbes.grid[,"C"]/Microbes.grid[,"N"])
    MicrobesPSeries <- rbind(MicrobesPSeries,Microbes.grid[,"C"]/Microbes.grid[,"P"])
    MicrobesNPSeries <- rbind(MicrobesNPSeries,Microbes.grid[,"N"]/Microbes.grid[,"P"])
    
    # Sum microbial pools
    Mic_Sum <- rbind(Mic_Sum,colSums(Microbes.grid))
    
    # Output grid if enzyme mass has peaked
    if (i_t>3 & params["output.mid",]==1) {
      if (Enzyme_Sum[i_t]<Enzyme_Sum[i_t-1] & Enzyme_Sum[i_t-1]>Enzyme_Sum[i_t-2]) {
        grid.max <- list("Substrates"=Substrates,"Monomers"=Monomers,"Enzymes"=Enzymes,"Microbes"=Microbes)
        grid.time <- i_t
      }
    }
    
    # Output grid images
    if (i_t == 1 || i_t%%params["print.grid",]==0){
      apply(matrix(1:n_taxa),1,function(i)plot.grid(Microbes[n_taxa*(0:(grid.size-1))+i,"C"],params["x",],zlim=72*fb[1:n_taxa][i]+3,paste(i,":",fb[1:n_taxa][i],":",i_t,sep="")))
    }
    
    # Sum substrate, monomer, enzyme pools
    Substrate_Sum <- c(Substrate_Sum,sum(Substrates.grid[,"C"]))
    Monomer_Sum <- c(Monomer_Sum,sum(Monomers.grid[,"C"]))
    Enzyme_Sum <- c(Enzyme_Sum,sum(Enzymes.grid))
    
    print(i_t)
  }
  
  # Calculate recoveries
  Recovered_C <- (sum(RespSeries) + sum(Substrates.grid[,"C"]) + sum(Monomers.grid[,"C"]) + sum(Enzymes.grid) + sum(Microbes.grid[,"C"]))-(Cum_Microbe["C"] + Cum_Substrate["C"] + Cum_Monomer["C"] + Cum_Enzyme_C)
  Recovered_N <- (Cum_Leaching_N + sum(Substrates.grid[,"N"]) + sum(Monomers.grid[,"N"]) + sum(Enzymes.grid*EnzAttrib[,"N_cost"]) + sum(Microbes.grid[,"N"]))-(Cum_Microbe["N"] + Cum_Substrate["N"] + Cum_Monomer["N"] + Cum_Enzyme_N)
  Recovered_P <- (Cum_Leaching_P + sum(Substrates.grid[,"P"]) + sum(Monomers.grid[,"P"]) + sum(Enzymes.grid*EnzAttrib[,"P_cost"]) + sum(Microbes.grid[,"P"]))-(Cum_Microbe["P"] + Cum_Substrate["P"] + Cum_Monomer["P"] + Cum_Enzyme_P)
  
  index <- is.finite(log10(apply(MicrobesSeries,2,max)))
  correl <- cor.test(rowSums(EnzGenes)[index],log(apply(MicrobesSeries,2,max)/grid.size)[index])  
  
  # Output grid at end if necessary
  if (is.null(grid.max)) {
    grid.max <- list("Substrates"=Substrates,"Monomers"=Monomers,"Enzymes"=Enzymes,"Microbes"=Microbes)
    grid.time <- params["end_time",]
  }
  
  # Correlations of taxon abundance from the grid
  corr.vec <- matrix(NA,ncol=2,nrow=n_taxa*n_taxa)
  count <- 1
  taxon.mat <- matrix(grid.max$Microbes[,"C"],ncol=n_taxa,byrow=T)
  for (i in 1:n_taxa) {
    for (j in 1:n_taxa) {
      taxon.a <- taxon.mat[,i]
      taxon.b <- taxon.mat[,j]
      overlap <- which(taxon.a>0 & taxon.b>0)
      taxon.1 <- taxon.a[overlap]
      taxon.2 <- taxon.b[overlap]
      if (length(taxon.1)>2) {
        corr.list <- cor.test(taxon.1,taxon.2)
      } else {
        corr.list <- list("estimate"=NA,"p.value"=NA)
      }
      corr.vec[count,1] <- corr.list$estimate
      corr.vec[count,2] <- corr.list$p.value*((params["n_taxa",]*params["n_taxa",]-params["n_taxa",])/2) # Bonferroni corrected p-values
      count <- count+1
    }
  }
  
  # Output pulse results
  out.pulse <- list(
    "params"=params,
    "seed"=rng.seed,
    "Ea"=Ea,
    "Uptake_Ea"=Uptake_Ea,
    "Vmax"=Vmax0,
    "Uptake_Vmax"=Uptake_Vmax0,
    "Km"=Km0[1:n_substrates,],
    "Uptake_Km"=Uptake_Km0[1:n_monomers,],
    "ReqEnz"=list(ReqEnz[,,1][1:n_substrates,],ReqEnz[,,2][1:n_substrates,]),
    "Uptake_ReqEnz"=Uptake_ReqEnz[1:n_monomers,],
    "OptimalRatios"=OptimalRatios,
    "RangeRatios"=RangeRatios,
    "EnzGenes"=EnzGenes,
    "EnzProdInduce"=EnzProdInduce,
    "EnzProdConstit"=EnzProdConstit,
    "UptakeGenes"=UptakeGenes[1:n_taxa,],
    "UptakeGenesForEnz"=UptakeGenesForEnz,
    "EnzAttrib"=EnzAttrib,
    "MonInput"=MonInput[1:n_monomers],
    "SubInput"=SubInput[1:n_substrates,],
    "CUE"=CUE.ref,
    "RespSeries"=RespSeries,
    "EnzymesSeries"=EnzymesSeries,
    "SubstratesSeries"=SubstratesSeries,
    "MonomersSeries"=MonomersSeries,
    "NH4Series"=NH4Series,
    "PO4Series"=PO4Series,
    "MicrobesSeries"=MicrobesSeries,
    "MicrobesNSeries"=MicrobesNSeries,
    "MicrobesPSeries"=MicrobesPSeries,
    "MicrobesNPSeries"=MicrobesNPSeries,
    "grid"=grid.max,
    "grid.time"=grid.time,
    "grid.size"=grid.size,
    "corr.vec"=corr.vec,
    "cor.test"=correl,
    "end_time"=params["end_time",],
    "Mic_Sum"=Mic_Sum,
    "Substrate_Sum"=Substrate_Sum,
    "Monomer_Sum"=Monomer_Sum,
    "Enzyme_Sum"=Enzyme_Sum,
    "timestamp"=timestamp,
    "fb"=fb[1:n_taxa],
    "Temp"=Temp,
    "Psi"=Psi,
    "Tolerance"=Tolerance,
    "RecoveredC_Resp_Sub_Mon_Enz_Mic"=
      c(Recovered_C,sum(RespSeries),sum(Substrates.grid[,"C"]),sum(Monomers.grid[,"C"]),sum(Enzymes.grid),sum(Microbes.grid[,"C"])),
    "RecoveredN_Leach_Sub_Mon_Enz_Mic"=
      c(Recovered_N,Cum_Leaching_N,sum(Substrates.grid[,"N"]),sum(Monomers.grid[,"N"]),sum(Enzymes.grid*EnzAttrib[,"N_cost"]),sum(Microbes.grid[,"N"])),
    "RecoveredP_Leach_Sub_Mon_Enz_Mic"=
      c(Recovered_P,Cum_Leaching_P,sum(Substrates.grid[,"P"]),sum(Monomers.grid[,"P"]),sum(Enzymes.grid*EnzAttrib[,"P_cost"]),sum(Microbes.grid[,"P"]))
  )
  out.pulse
} # End of RunPulse() #########################################################

###############################################################################
TraitModel <- function(job.time,task.ID){
  
  # Read in parameters
  params.file <- paste("params/params",task.ID,".txt",sep="")
  params1 <- read.table(params.file, header=F, row.names=1, sep="\t",stringsAsFactors=F)
  params <- params1
  grid.size <- params["x",]*params["y",]

  seed.val <- as.POSIXct(job.time,format="%y%m%d%H%M%S") + task.ID # convert to date integer, add task id, and use as seed
  set.seed(seed.val) 
  timestamp <- format(seed.val, "%y%m%d%H%M%S") # convert date seed to 12 digit number
  rng.seed <- as.numeric(timestamp)
  if(params["set.seed",]==1) {
    set.seed(as.POSIXct(as.character(params["seed",]),format="%y%m%d%H%M%S"))
    rng.seed <- params["seed",]
  }
  
  n_monomers <- params["n_substrates",]+2
  n_genes <- params["n_enzymes",]
  n_upgenes <- params["n_uptake",]
  
  # Read in activation energies for substrates
  Ea.frame <- read.table("Ea.txt", header=TRUE, row.names=1, sep="\t",stringsAsFactors=F)
  
  # Choose taxa with "fungal" strategy, specifically ability to translocate C and nutrients
  fb <- sample(c(1,0),params["n_taxa",],replace=T,prob=c(params["fb",],(1-params["fb",])))

  # Set stress tolerance allocation for each taxon
  Tolerance <- c(runif(params["n_taxa",],params["Tol_min",],params["Tol_max",]))
  if(params["Tol_min",]==params["Tol_max",]) {dummy <- runif(params["n_taxa",])}
  
  # Pre-exponential constants for enzymes
  # Rows are substrates; cols are enzymes
  Vmax0 <- c(runif(params["n_substrates",]*params["n_enzymes",],params["Vmax0_min",],params["Vmax0_max",]))
  Vmax0 <- matrix(Vmax0,nrow=params["n_substrates",],ncol=params["n_enzymes",],byrow=T,dimnames=list(sprintf("%s%03d","Sub",1:params["n_substrates",]),sprintf("%s%03d","Enz",1:params["n_enzymes",])))

  # Pre-exponential constants for uptake
  # Rows are monomers; cols are uptake enzymes
  Uptake_Vmax0 <- c(runif(params["n_uptake",]*n_monomers,params["Uptake_Vmax0_min",],params["Uptake_Vmax0_max",]))
  Uptake_Vmax0 <- matrix(Uptake_Vmax0,nrow=n_monomers,ncol=params["n_uptake",],byrow=T,dimnames=list(sprintf("%s%03d","Mon",1:n_monomers),sprintf("%s%03d","Upt",1:params["n_uptake",])))
  
  # Enzyme specificity matrix of activation energies
  # Rows are substrates; cols are enzymes
  Ea <- sapply(seq(len=dim(Ea.frame)[1]),function(i)runif(params["n_enzymes",],Ea.frame$Ea_min[i],Ea.frame$Ea_max[i]))
  Ea <- matrix(Ea,nrow=params["n_substrates",],ncol=params["n_enzymes",],byrow=T,dimnames=list(sprintf("%s%03d","Sub",1:params["n_substrates",]),sprintf("%s%03d","Enz",1:params["n_enzymes",])))
  
  # Uptake specificity matrix of activation energies
  # Rows are monomers; cols are uptake enzymes
  Uptake_Ea <- c(runif(params["n_uptake",]*n_monomers,params["Uptake_Ea_min",],params["Uptake_Ea_max",]))
  Uptake_Ea <- matrix(Uptake_Ea,nrow=n_monomers,ncol=params["n_uptake",],byrow=T,dimnames=list(sprintf("%s%03d","Mon",1:n_monomers),sprintf("%s%03d","Upt",1:params["n_uptake",])))
  
  # Enzymes required for substrate degradation
  # Rows are substrates; cols are enzymes
  # Same number within row implies redundancy
  # Ensures each substrate is degraded by at least 1 enzyme and every enzyme degrades at least 1 substrate
  probability_vector <- rep(0,params["n_enzymes",])
  probability_vector[1:params["Enzymes_per_sub",]] <- 1
  ReqEnz1 <- NULL
  for(i in 1:params["n_substrates",]) {
    ReqEnz1 <- c(ReqEnz1,sample(probability_vector,params["n_enzymes",],replace=F))
  }
  ReqEnz1 <- matrix(ReqEnz1,nrow=params["n_substrates",],ncol=params["n_enzymes",],byrow=T)
  probability_vector <- rep(0,params["n_substrates",])
  probability_vector[1] <- 1
  for (i in 1:params["n_enzymes",]) {
    if(sum(ReqEnz1[,i])==0) ReqEnz1[,i] <- sample(probability_vector,params["n_substrates",],replace=F)
  }
  
  # Choose some substrates that require multiple enzymes
  probability_vector <- rep(0,params["n_enzymes",])
  if(params["Avg_extra_req_enz",]>0) probability_vector[1:params["Avg_extra_req_enz",]] <- 1
  ReqEnz2 <- sample(probability_vector,params["n_substrates",]*params["n_enzymes",],replace=T)
  ReqEnz2[ReqEnz1==1] <- 0
  ReqEnz2 <- matrix(ReqEnz2,nrow=params["n_substrates",],ncol=params["n_enzymes",],byrow=T)
  for(j in 1:dim(ReqEnz2)[1]) {
    # Put in NAs if the substrate does not require a second enzyme
    if(rowSums(ReqEnz2)[j]==0) ReqEnz2[j,] <- NA
  }
  
  # Generate the correct structure for the required enzyme matrices
  ReqEnz <- c(expand(ReqEnz1,grid.size),expand(ReqEnz2,grid.size))
  # rows, cols, stacks; default assemble by cols
  ReqEnz <- array(ReqEnz,dim=c(grid.size*params["n_substrates",],params["n_enzymes",],2),dimnames=list(rep(sprintf("%s%03d","Sub",1:params["n_substrates",]),grid.size),sprintf("%s%03d","Enz",1:params["n_enzymes",]),c("set1","set2")))
    
  # Enzymes used for uptake
  # Rows are monomers; cols are uptake enzymes
  # Same number within row implies redundancy
  # Make sure each monomer is taken up by at least one transporter and every transporter takes up at least one monomer
  probability_vector <- rep(0,params["n_uptake",])
  probability_vector[1:params["Uptake_per_monomer",]] <- 1
  Uptake_ReqEnz1 <- NULL
  for(i in 1:n_monomers) {
    Uptake_ReqEnz1 <- c(Uptake_ReqEnz1,sample(probability_vector,params["n_uptake",],replace=F))
  }
  Uptake_ReqEnz <- matrix(Uptake_ReqEnz1,nrow=n_monomers,ncol=params["n_uptake",],byrow=T,dimnames=list(sprintf("%s%03d","Mon",1:n_monomers),sprintf("%s%03d","Upt",1:params["n_uptake",])))
  probability_vector <- rep(0,n_monomers)
  probability_vector[1] <- 1
  for (i in 1:params["n_uptake",]) {
    if(sum(Uptake_ReqEnz[,i])==0) Uptake_ReqEnz[,i] <- sample(probability_vector,n_monomers,replace=F)
  }
  
  # Substrate concentrations
  substrates.frame <- read.table("substrates.txt", header=TRUE, sep="\t",stringsAsFactors=F,row.names=1)
  Substrates <- data.matrix(substrates.frame)
 
  # Substrate input rates
  SubInputC <- read.table("inputs.txt", header=TRUE, sep="\t",stringsAsFactors=F,row.names=1)[,1]
  SubInputN <- SubInputC*Substrates[,"N"]/Substrates[,"C"]
  SubInputP <- SubInputC*Substrates[,"P"]/Substrates[,"C"]
  SubInput <- c(SubInputC,SubInputN,SubInputP)
  SubInput <- matrix(SubInput,nrow=params["n_substrates",],ncol=3,byrow=F,dimnames=list(sprintf("%s%03d","SubIn",1:params["n_substrates",]),c("C","N","P")))
  SubInput["SubIn001",] <- SubInput["SubIn002",] <- 0
  Substrates["DeadMic",] <- Substrates["DeadEnz",] <- 0
  
  # Monomers produced by each substrate (binary)
  # Rows are substrates; cols are monomers
  # All rows should sum to 1
  MonomersProduced <- cbind(rep(0,params["n_substrates",]),rep(0,params["n_substrates",]),diag(params["n_substrates",]))
  MonomersProduced <- as.matrix(MonomersProduced)
  dimnames(MonomersProduced) <- list(sprintf("%s%03d","Sub",1:params["n_substrates",]),sprintf("%s%03d","Mon",-1:(n_monomers-2)))
  dimnames(MonomersProduced)[[2]][1:4]<- c("NH4","PO4","DeadMic","DeadEnz")
  
  # Optimal stoichiometry of bacterial taxa
  OptimalRatios <- rep(c(params["Cfrac_b",],params["Nfrac_b",],params["Pfrac_b",]),params["n_taxa",])
  OptimalRatios <- matrix(OptimalRatios,nrow=params["n_taxa",],ncol=3,byrow=T,dimnames=list(sprintf("%s%03d","Tax",1:params["n_taxa",]),c("C","N","P")))
  # Substitute fungal stoichiometry
  OptimalRatios[fb==1,] <- matrix(rep(c(params["Cfrac_f",],params["Nfrac_f",],params["Pfrac_f",]),sum(fb)),ncol=3,byrow=T)
  
  RangeRatios <- rep(c(params["Crange",],params["Nrange",],params["Prange",]),params["n_taxa",])
  RangeRatios <- matrix(RangeRatios,nrow=params["n_taxa",],ncol=3,byrow=T,dimnames=list(sprintf("%s%03d","Tax",1:params["n_taxa",]),c("C","N","P")))
  RangeWeight <- RangeRatios/OptimalRatios
  
  # Calculate minimum cell quotas
  MinRatios <- OptimalRatios-RangeRatios
  
  # Enzyme genes possessed by taxa (binary)
  # Rows are taxa; cols are genes
  # Not all taxa need have an enzyme gene
  # This method draws the number of genes per taxon from a uniform distribution
  # and allows production rates between Enz_Prod_min and max
  genes_per_taxon <- sample(params["Enz_per_taxon_min",]:params["Enz_per_taxon_max",],params["n_taxa",],replace=T)
  EnzGenes <- sapply(seq(params["n_taxa",]),function(i){
    probability_vector <- rep(0,n_genes)
    probability_vector[0:genes_per_taxon[i]] <- 1
    sample(probability_vector,n_genes,replace=F)})
  EnzGenes <- matrix(EnzGenes,nrow=params["n_taxa",],ncol=n_genes,byrow=T,dimnames=list(sprintf("%s%03d","Tax",1:params["n_taxa",]),sprintf("%s%03d","Enz",1:n_genes)))
  Normalize <- 1
  if (params["NormalizeProd",]==1) {
    Normalize <- rowSums(EnzGenes)/params["Enz_per_taxon_max",]
  }
  EnzProdInduce <- apply(EnzGenes,1,function(x)runif(1,params["Enz_Prod_min",],params["Enz_Prod_max",]))
  EnzProdInduce <- EnzProdInduce*EnzGenes/Normalize
  EnzProdInduce[is.na(EnzProdInduce)] <- 0
  EnzProdConstit <- apply(EnzGenes,1,function(x)runif(1,params["Constit_Prod_min",],params["Constit_Prod_max",]))
  EnzProdConstit <- EnzProdConstit*EnzGenes/Normalize
  EnzProdConstit[is.na(EnzProdConstit)] <- 0
  
  # Uptake genes possessed by taxa (binary)
  # Rows are taxa; cols are genes
  # All taxa must have at least one uptake gene
  # For each taxon, check if it has the enz gene associated with the monomer(s) to be targeted by the uptake gene
  # If so, assign that uptake gene to the taxon with probability p
  # Enz-Taxon%*%Enz-Sub%*%Sub-Monomer to get Taxon-Monomer matrix
  # Matrix multiply this by the required uptake enzyme matrix, then assign uptake genes
  RE2 <- ReqEnz[,,2][1:params["n_substrates",],]
  RE2[is.na(RE2)] <- 0
  enz_sub <- ReqEnz[,,1][1:params["n_substrates",],] + RE2
  # Matrix multiplication to relate taxa to the monomers they can generate with their enzymes
  UptakeGenes <- EnzGenes%*%t(enz_sub)%*%MonomersProduced
  UptakeGenes[,1:2] <- 1
  UptakeGenes[UptakeGenes>0] <- 1
  # Make sure every taxon is likely to have an uptake enzyme for at least 1 organic monomer
  # Not guaranteed unless uptake_prob = 1
  probability_vector <- rep(0,n_monomers-2)
  probability_vector[1] <- 1
  for (i in 1:params["n_taxa",]) {
    if(sum(UptakeGenes[i,3:n_monomers])==0) UptakeGenes[i,3:n_monomers] <- sample(probability_vector,n_monomers-2,replace=F)
  }

  # Give taxa random number of additional uptake genes between the number they have and n_upgenes
  for (i in 1:params["n_taxa",]) {
    n.zero <- length(UptakeGenes[i,][UptakeGenes[i,]==0])
    probability_vector <- rep(0,n.zero)
    probability_vector[1:sample(0:n.zero,1)] <- 1
    UptakeGenes[i,][UptakeGenes[i,]==0] <- sample(probability_vector,n.zero,replace=F)
  }
  UptakeGenesForEnz <- UptakeGenes
  # If true then the uptake potential is normalized to the number of uptake genes
  Normalize <- 1
  if (params["NormalizeUptake",]==1) {
    Normalize <- rowSums(UptakeGenes)/params["n_uptake",]
  }
  UptakeGenes <- UptakeGenes/Normalize
  # Choose total amount of uptake allocation at random
  UptakeProd <- c(runif(params["n_taxa",],params["Uptake_C_cost_min",],params["Uptake_C_cost_max",]))
  UptakeGenes <- UptakeGenes*UptakeProd
    
  # Calculate CUE as a function of the number of enzyme genes, uptake genes, and reference CUE
  CUE.ref <- params["CUE_enz",]*rowSums(EnzGenes)/(params["Enz_per_taxon_max",]) + 
    params["CUE_uptake",]*rowSums(UptakeGenesForEnz)/(params["n_uptake",]) + params["CUE_ref",] + params["Tol_CUE",]*Tolerance
  
  # Enzyme attributes
  # Rows are enzymes; cols are attributes
  EnzAttrib_C <- rep(params["Enz_C_cost",],params["n_enzymes",])
  EnzAttrib <- rbind(EnzAttrib_C,params["Enz_N_cost",],params["Enz_P_cost",],params["Enz_Maint_cost",])
  EnzAttrib <- matrix(EnzAttrib,nrow=params["n_enzymes",],ncol=4,byrow=T,dimnames=list(sprintf("%s%03d","Enz",1:params["n_enzymes",]),c("C_cost","N_cost","P_cost","Maint_cost")))
  
  # Monomer pool sizes for all elements
  Monomers <- rbind(c(0,params["Init_NH4",],0),c(0,0,params["Init_PO4",]),Substrates*(params["Monomer_Substrate_Ratio",]))
  Monomers[is.na(Monomers)] <- 0
  Monomers <- as.matrix(Monomers)
  dimnames(Monomers) <- list(c("NH4","PO4","DeadMic","DeadEnz",sprintf("%s%03d","Mon",3:(n_monomers-2))),c("C","N","P"))
  
  # Monomer input rates
  MonInput <- c(params["Input_NH4",],params["Input_PO4",],read.table("inputs.txt", header=TRUE, sep="\t",stringsAsFactors=F,row.names=1)[,2])
  MonInput <- as.matrix(MonInput)
  rownames(MonInput) <- rownames(Monomers)
  
  # Enzyme pool sizes for all elements
  Enzymes <- matrix(c(runif(params["n_enzymes",],params["Enz_min",],params["Enz_max",])),nrow=params["n_enzymes",],ncol=1,dimnames=list(sprintf("%s%03d","Enz",1:params["n_enzymes",]),"C"))
  
  # Microbial pool sizes for all elements
  BacC <- 0.5*params["max_size_b",]
  BacN <- BacC*params["Nfrac_b",]/params["Cfrac_b",]
  BacP <- BacC*params["Pfrac_b",]/params["Cfrac_b",]
  FunC <- 0.5*params["max_size_f",]
  FunN <- FunC*params["Nfrac_f",]/params["Cfrac_f",]
  FunP <- FunC*params["Pfrac_f",]/params["Cfrac_f",]
  Microbes <- rep(c(BacC,BacN,BacP),params["n_taxa",]*grid.size)
  Microbes <- matrix(Microbes,nrow=params["n_taxa",]*grid.size,ncol=3,byrow=T,dimnames=list(rep(sprintf("%s%03d","Tax",1:params["n_taxa",]),grid.size),c("C","N","P")))
  fb <- rep(fb,grid.size)
  Microbes[fb==1,"C"] <- FunC; Microbes[fb==1,"N"] <- FunN; Microbes[fb==1,"P"] <- FunP
  # Randomly assign microbes to each grid box
  p.b <- params["taxa_per_box",]
  p.f <- p.b*params["max_size_b",]/params["max_size_f",]
  choose_taxa <- sample(c(1,0),params["n_taxa",]*grid.size,replace=T,c(p.b,(1-p.b)))
  choose_taxa[fb==1] <- sample(c(1,0),sum(fb),replace=T,c(p.f,(1-p.f)))
  Microbes[choose_taxa==0,] <- 0
  # Initialize the reproduction list
  ReprodNew <- Microbes
  ReprodNew[ReprodNew>0] <- 0
  Colonization.reset <- ReprodNew
  
  # Account for efficiency-specificity tradeoff by dividing Vmax_0 by the number of substrates (or monomers) targeted
  # and multiplied by a specificity factor
  # Leave Km unchanged, effectively increasing it by the factor that Vmax_0 is reduced
  total_substrates <- colSums(ReqEnz[,,1][1:params["n_substrates",],]) + colSums(RE2)
  if (params["Specif_factor",]==0) {total_substrates[total_substrates>1] <- 1} else
  {total_substrates[total_substrates>1] <- total_substrates[total_substrates>1]*params["Specif_factor",]}
  Vmax0 <- t(t(Vmax0)/total_substrates)
  Vmax0[!is.finite(Vmax0)] <- 0
  
  total_monomers <- colSums(Uptake_ReqEnz)
  if (params["Specif_factor",]==0) {total_monomers[total_monomers>1] <- 1} else
  {total_monomers[total_monomers>1] <- total_monomers[total_monomers>1]*params["Specif_factor",]}
  Uptake_Vmax0 <- t(t(Uptake_Vmax0)/total_monomers)
  Uptake_Vmax0[!is.finite(Uptake_Vmax0)] <- 0
  
  # Implement Vmax-Km tradeoff as a direct correlation with slope Vmax_Km, 
  # and error term normally distributed with magnitude Km_error. Minimum Km constrained to Km_min
  Km <- abs(rnorm(Vmax0, mean=Vmax0*params["Vmax_Km",], sd = mean(Vmax0)*params["Km_error",])+params["Vmax_Km_int",])
  Km[Km < params["Km_min",]] <- params["Km_min",]
  Km <- matrix(Km,nrow=params["n_substrates",],ncol=params["n_enzymes",],byrow=F,dimnames=list(sprintf("%s%03d","Sub",1:params["n_substrates",]),sprintf("%s%03d","Enz",1:params["n_enzymes",])))
  Uptake_Km <- abs(rnorm(Uptake_Vmax0, mean=Uptake_Vmax0*params["Uptake_Vmax_Km",], sd = mean(Uptake_Vmax0)*params["Km_error",])+params["Uptake_Vmax_Km_int",])
  Uptake_Km[Uptake_Km < params["Uptake_Km_min",]] <- params["Uptake_Km_min",]
  Uptake_Km <- matrix(Uptake_Km,nrow=n_monomers,ncol=params["n_uptake",],byrow=F,dimnames=list(sprintf("%s%03d","Mon",1:n_monomers),sprintf("%s%03d","Upt",1:params["n_uptake",])))
  
  # Expand matricies to cover the whole grid
  Substrates <- expand(Substrates,grid.size)
  SubInput <- expand(SubInput,grid.size)
  Monomers <- expand(Monomers,grid.size)
  MonInput <- expand(MonInput,grid.size)
  Enzymes <- expand(Enzymes,grid.size)
  Km0 <- expand(Km,grid.size)
  Uptake_Km0 <- expand(Uptake_Km,grid.size)
  Uptake_ReqEnz <- expand(Uptake_ReqEnz,grid.size)
  MinRatios <- expand(MinRatios,grid.size)
  
  climate <- read.table("climate2.txt", header=TRUE, sep="\t",stringsAsFactors=F)
  
  # Set up output list containing each pulse
  out <- vector("list",params["pulses",])
  for (i_p in 1:params["pulses",]) {
    zeros <- ifelse(i_p<10,"0","")
    filename <- paste("outputs7/", timestamp, "_", zeros, i_p, "G.png", sep = "")
    rows <- params["end_time",]%/%params["print.grid",] + 1
    png(file=filename,width=300*params["n_taxa",],height=300*rows)
    par(mfrow=c(rows,params["n_taxa",]),mar=c(3,3,3,3))
    Temp <- climate$Temp[(1+params["end_time",]*(i_p-1)):(params["end_time",]*i_p)]
    Psi <- climate$Psi[(1+params["end_time",]*(i_p-1)):(params["end_time",]*i_p)]
    out[[i_p]] <- RunPulse(
      params,
      timestamp,
      rng.seed,
      grid.size,
      Microbes,
      Substrates,
      SubInput,
      Enzymes,
      Monomers,
      MonInput,
      MonomersProduced,
      ReprodNew,
      Colonization.reset,
      Ea,
      Vmax0,
      Km0,
      ReqEnz,
      EnzGenes,
      EnzProdInduce,
      EnzProdConstit,
      UptakeGenes,
      UptakeGenesForEnz,
      Uptake_ReqEnz,
      EnzAttrib,
      Uptake_Ea,
      Uptake_Vmax0,
      Uptake_Km0,
      CUE.ref,
      OptimalRatios,
      MinRatios,
      RangeRatios,
      fb,
      Temp,
      Psi,
      Tolerance
    )
    dev.off()
    filename <- paste("outputs7/", timestamp, "_", zeros, i_p, ".png", sep = "")
    png(file=filename,width=2000,height=800)
    MakePlots(out[[i_p]])
    dev.off()
    # Recalculate microbial frequencies and repopulate the grid with microbes based on frequencies
    cum.abundance <- apply(out[[i_p]]$MicrobesSeries, 2, sum)
    # Account for different cell sizes of bacteria and fungi
    cum.abundance[out[[i_p]]$fb==1] <- cum.abundance[out[[i_p]]$fb==1]*params["max_size_b",]/params["max_size_f",]
    frequencies <- cum.abundance/sum(cum.abundance)
    
    # Microbial pool sizes for all elements
    BacC <- 0.5*params["max_size_b",]
    BacN <- BacC*params["Nfrac_b",]/params["Cfrac_b",]
    BacP <- BacC*params["Pfrac_b",]/params["Cfrac_b",]
    FunC <- 0.5*params["max_size_f",]
    FunN <- FunC*params["Nfrac_f",]/params["Cfrac_f",]
    FunP <- FunC*params["Pfrac_f",]/params["Cfrac_f",]
    Microbes <- rep(c(BacC,BacN,BacP),params["n_taxa",]*grid.size)
    Microbes <- matrix(Microbes,nrow=params["n_taxa",]*grid.size,ncol=3,byrow=T,dimnames=list(rep(sprintf("%s%03d","Tax",1:params["n_taxa",]),grid.size),c("C","N","P")))
    Microbes[fb==1,"C"] <- FunC; Microbes[fb==1,"N"] <- FunN; Microbes[fb==1,"P"] <- FunP
    # Randomly assign microbes to each grid box based on prior densities
    probs <- matrix(c(frequencies,(1-frequencies)),ncol=2)
    choose_taxa <- matrix(rep(0,grid.size*params["n_taxa",]),nrow=params["n_taxa",])
    for (i in 1:params["n_taxa",]) {
      choose_taxa[i,] <- sample(c(1,0),grid.size,replace=T,probs[i,])
    }
    Microbes[as.vector(choose_taxa)==0,] <- 0
    # Reset reproduction
    ReprodNew[ReprodNew>0] <- 0  
  }
  out
  
} # End of TraitModel() #######################################################

