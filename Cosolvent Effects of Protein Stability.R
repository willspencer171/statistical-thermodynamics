# This Script is intended for the use of automating the production of publish-quality graphs
# from a series of datasets collected from various sources

if(!("tidyverse" %in% installed.packages())){
  install.packages("tidyverse")
}

#These are the packages that we will be using through this analysis
library(ggplot2)
library(scales)
library(rbenchmark)
options(warn = -1)

# Set the working directory
setwd("D:/University of York/Biochemistry/Final Year Project/Data Files") 
rm(list = ls())

# Prediction resolution
Res <- 30

# Constants
rConst <- 8.314
Kelvin <- 273

# Set target folder with data files in
targetFolder <- list.files(path = paste(getwd(), "/Data Files", sep = ""), pattern = "*.csv", full.names = F)


# This loop generates the graph and completes the statistical analysis of data
for (item in targetFolder) {
  data <- read.csv(T, file = paste("./Data Files/", item, sep = ""), fileEncoding = 'UTF-8-BOM')
  ratio <- data$Ratio[1]
  item <- list(substr(item, 1, nchar(item)-4))
  names(item) <- c("Name of Paper")
  combiPresent <- 1
  deltaHPresent <- 1

  ### Predicting the stabiliser
  mod <- na.exclude(lm(data = data, SY~SX));summary(mod)
  pred <- predict(mod);pred
  item$StabPred <- seq(from = pred[1], to = tail(pred, 1), by = abs(pred[1]-tail(pred ,1))/(Res-1))
  item$StabR2 <- summary(mod)$r.squared # R squared value for the raw data
  
  ### Predicting the denaturant
  mod <- na.exclude(lm(data = data, DY~DX));summary(mod)
  pred <- predict(mod);pred
  item$DenatPred <- seq(from = pred[1], to = tail(pred, 1), by = -(pred[1]-tail(pred ,1))/(Res-1))
  item$DenatR2 <- summary(mod)$r.squared
  
  ### Predicting the Experimental Combination
  tryCatch({
    mod <- na.exclude(lm(data = data, CY~CX));summary(mod)
    pred <- predict(mod);pred
    item$ExpC <- seq(from = pred[1], to = tail(pred, 1), by = -(pred[1]-tail(pred ,1))/(Res-1))
    item$ExpCR2 <- summary(mod)$r.squared
    },
    error=function(cond) {
      combiPresent <<- 0
      message(paste(item$`Name of Paper`, "No experimental combination found", sep = " - "))
    }
    )

  ### Producing Estimate of Combination (dilute assumption)
  item$ComboPred <- (item$StabPred + item$DenatPred)/2
  if (is.na(data$DeltaH[1])){
    deltaHPresent <<- 0
    message(paste(item$`Name of Paper`, "No DeltaH found", sep = " - "))
  }else{
    item$deltaH <- data$DeltaH[1]
    item$deltaG <- 0
    for (i in 1:length(data$chem[data$chem != ""])){

      order <- c("DY", "SY", "DX", "SX")
      dTm <- tail(na.exclude(data[order[i]]), 1) - data[[order[i]]][1];
      item$deltaG[i] <- -(data$DeltaH * dTm/(rConst * (data[[order[i]]][1]) * (1/tail(na.exclude(data[order[i+2]]), 1))))
      
      
      names(item$deltaG)[i] <- paste(data$chem[i], "DeltaG")
      i <- i + 1
    }
  }

  # lower resolution defined just for graph aesthetics
  CPLowRes <- c()
  for (i in 1:length(data$DX)) {
    CPLowRes <- c(CPLowRes, item$ComboPred[(Res/length(data$DX))*i])
    i <- i+1
  }
  
  ### The graph
  graph <- ggplot(data = data)+
      
    # Denaturant
    geom_point(aes(x = DX, y = DY+Kelvin), shape = 1, size = 3)+
    stat_smooth(method = "lm", aes(x = DX, y = DY+Kelvin), se = F, 
                colour = "red", 
                size = 0.5)+
    
    # Stabiliser
    geom_point(aes(x = SX*ratio, y = SY+Kelvin), shape = 2, size = 3)+
    stat_smooth(method = "lm", aes(x = SX*ratio, y = SY+Kelvin), se =F, 
                colour = "green", 
                size = 0.5)+
    
    # Combination (predictions)
    geom_point(aes(x = DX, y = CPLowRes+Kelvin), shape = 3, size = 3)+
    stat_smooth(method = "lm", aes(x = DX, y = CPLowRes+Kelvin), se = F, 
                colour = "purple", 
                size = 0.5, linetype = "dashed")+
    
    # Adding R Squared values as annotations
    # For the Stabiliser
    annotate("text",
             x = data$SX[ceiling(length(na.exclude(data$SX))-1)]*ratio,
             y = data$SY[ceiling(length(na.exclude(data$SY))-1)]+2+Kelvin,
             label = substitute(R^2 == rVal, list(rVal = signif(item$StabR2, digits = 3))),
             colour = "green",
             size = 5)+
    
    # For the Denaturant
    annotate("text",
             x = data$DX[length(na.exclude(data$DX))-1],
             y = data$DY[floor(length(na.exclude(data$DY))-1)]-2+Kelvin,
             label = substitute(R^2 == rVal, list(rVal = signif(item$DenatR2, digits = 3))),
             colour = "red",
             size = 5)+
    
    annotate("text",
             x = 0,
             y = min(c(data$SY, data$DY), na.rm = T)+Kelvin,
             label = paste(item$`Name of Paper`, ": ", data$chem[1], " vs ", data$chem[2], sep = ""),
             hjust = 0,
             size = 6)+
    
    # Graph Aesthetics
    scale_x_continuous(name = paste("[", data$chem[1], "]", " // M", sep = ""), 
                       labels = number_format(accuracy = 0.01), 
                       sec.axis = sec_axis(~. /2, name = paste("[", data$chem[2], "]", " // M", sep = ""), labels = number_format(accuracy = 0.01)))+
    scale_y_continuous(name = expression("Melting Temperature (T"[" m"]~")" ~ "// K"), 
                       breaks = seq(signif(min(c(data$SY, data$DY), na.rm = T), digits = 2)+Kelvin, 
                                    signif(max(c(data$SY, data$DY), na.rm = T), digits = 2)+Kelvin, 
                                    4))+
    theme_bw()+
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 14))
  
  if(combiPresent == 1){
    graph <- graph + geom_point(aes(x = CX, y = CY+Kelvin), shape = 0, size = 3)+
    stat_smooth(method = "lm", aes(x = CX, y = CY+Kelvin), se = F, colour = "black", size = 0.5)+
    
    annotate("text",
             x = data$CX[ceiling(length(na.exclude(data$CX))-1)],
             y = data$CY[ceiling(length(na.exclude(data$CY))/2)]-1+Kelvin,
             label = substitute(R^2 == rVal, list(rVal = signif(item$ExpCR2, digits = 3))),
             colour = "black",
             size = 5)
  }
 

  if(deltaHPresent == 1 & combiPresent == 0){
    names(item) <- c("Name of Paper",
                     "Stabiliser Model",
                     "Stabiliser Model R Squared",
                     "Denaturant Model",
                     "Denaturant Model R Squared",
                     "Combination Model",
                     "DeltaH")
  }else if(deltaHPresent == 1 & combiPresent == 1){
    names(item) <- c("Name of Paper",
                     "Stabiliser Model",
                     "Stabiliser Model R Squared",
                     "Denaturant Model",
                     "Denaturant Model R Squared",
                     "Experimental Combination Fit",
                     "Experimental Combination Fit R Squared",
                     "Combination Model",
                     "DeltaH")
  }else if(deltaHPresent == 0 & combiPresent == 1){
    names(item) <- c("Name of Paper",
                     "Stabiliser Model",
                     "Stabiliser Model R Squared",
                     "Denaturant Model",
                     "Denaturant Model R Squared",
                     "Experimental Combination Fit",
                     "Experimental Combination Fit R Squared",
                     "Combination Model")
  }
  # Compiles all the predictive data along with the name of the file, and the R squared values associated
  write.csv(item[-c(1)], paste("./Outputs/", item$`Name of Paper`, "_Data Output.csv", sep = ""))
  
  # Graph output
  png(paste("./Outputs/Images/", item, "_colour.png", sep=""), width = 1080, height = 1080, bg = NA, res = 150)
  print(graph)
  dev.off()

}

