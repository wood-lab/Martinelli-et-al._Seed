
### WORKING DIRECTORY AND PACKAGES
#################################################################################
setwd('/Users/jmartine/Dropbox/My Mac (Julietas-MacBook-Air.local)/Desktop/')

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(RColorBrewer) #display.brewer.all()
library(rcompanion)
library(wesanderson)
library(lme4)
library(car)
library(ggeffects)

### LOADING DATA
data <- read.csv('seed_data_tweaked.csv', header=TRUE, row.names = NULL, stringsAsFactors=FALSE,  sep=',') 
str(data)

# subsetting to remove HI
data1 <- subset(data, data$State == "WA"|data$State == "AK") # removing HI
dataHI <- subset(data, data$State == "HI") # ONLY KEEPING HI

# subsetting to remove HI-2: "data2"
data2 <- subset(data, data$FLUPSY == "WA-1"|data$FLUPSY == "WA-2"|
                  data$FLUPSY == "WA-3"| data$FLUPSY == "WA-4"|
                  data$FLUPSY == "HI-1"|data$FLUPSY == "AK-1") # removing HI-2

# subsetting to remove AK and HI: "data3"
data3 <- subset(data, data$State == "WA") # selecting only WA

### ABUNDANCE PER FLUPSY
abund_wa <- data3 %>%
  group_by(FLUPSY) %>% #
  summarize_at(vars(Count), list(sum = sum), na.rm=TRUE)

### CALCULATING PREVALENCE PER STATE AND FARM/FLUPSY
data_state <- data %>%
        group_by(State) %>% #
        summarize_at(vars(Infested), list(mean = mean), na.rm=TRUE)

data_farm <- data %>%
        group_by(State, FLUPSY) %>% #
        summarize_at(vars(Infested), list(mean = mean), na.rm=TRUE)
farm_prev <- as.data.frame(data_farm)
Prev <- (farm_prev$mean)*100
farm_prev <- cbind(farm_prev, Prev)

### just for WA and AK
data1_prev <- data1 %>%
  group_by(State, FLUPSY) %>% #
  summarize_at(vars(Infested), list(mean = mean), na.rm=TRUE)

data1_prev <- as.data.frame(data1_prev)
Prev <- (data1_prev$mean)*100
akwa_prev <- cbind(data1_prev, Prev)

### just for HI
dataHI_prev <- dataHI %>%
  group_by(State, FLUPSY) %>% #
  summarize_at(vars(Infested), list(mean = mean), na.rm=TRUE)

dataHI_prev <- as.data.frame(dataHI_prev)
Prev <- (dataHI_prev$mean)*100
hi_prev <- cbind(dataHI_prev, Prev)

### PLOTTING PREVALENCE
akwa_plot <- ggplot(akwa_prev, aes(x= FLUPSY, y=Prev, fill= State, show.legend = FALSE))  +
  # scale_fill_manual(values=wes_palette("GrandBudapest1")) + 
  scale_fill_brewer(palette = "YlGnBu") + 
  geom_bar(stat = "identity", col='black', show.legend = FALSE) +
  ylim(0,20) +
  labs(x = 'ALL FLUPSY', y = 'Prevalence (%)', size=16) 
akwa_plot + theme_classic(base_size = 18) 


hi_plot <- ggplot(hi_prev, aes(x= FLUPSY, y=Prev, fill= State, show.legend = FALSE))  +
       scale_fill_manual(values=wes_palette("GrandBudapest1")) + 
        #scale_fill_brewer(palette = "Spectral") +
        geom_bar(stat = "identity", col='black', show.legend = FALSE) +
        ylim(0,100) +
        labs(x = 'FLUPSY and bottle culture', y = 'Prevalence (%)', size=16) 
hi_plot + theme_classic(base_size = 18) 


### TESTING DIFF HAWAI'I
# SUBSETTING Treated (H1-2) vs Untreated seed (H1-1)
hi_flupsy <- subset(data, data$State == "HI") # selecting only HI
hi_flupsy1 <- subset(hi_flupsy, hi_flupsy$FLUPSY == "HI-1") # selecting only HI
hi_flupsy2 <- subset(hi_flupsy, hi_flupsy$FLUPSY == "HI-2") # selecting only HI

### NON-PARAMETRIC WILCOX TEST
wilcox.test(hi_flupsy1$Infested, y=hi_flupsy2$Infested) # W = 5000, p-value < 2.2e-16


### MODEL TO TEST OTHER STATES
# This model excludes the treated FLUPSY from Hawai'i
# A second model only for WA will also be tested

### MODEL: WITHOUT HI
model <- glm(Infested ~ FLUPSY + Shell.height..mm., family="binomial", data = data1)
summary(model)
anova(model)
car::Anova(model, type=3) # getting p-values 

### MODEL 1: WITHOUT HI-2
model1 <- glm(Infested ~ FLUPSY + Shell.height..mm., family="binomial", data = data2)
summary(model1)
anova(model1)
car::Anova(model1, type=3) # getting p-values 

### MODEL 2: ONLY WA
model2 <- glm(Infested ~ FLUPSY + Shell.height..mm., family="binomial", data = data3)
summary(model2)
anova(model2)
vif(model2)
car::Anova(model2, type=3) # getting p-values 

## FLUPSY SIG, HEIGHT NOT
plotmod <- ggpredict(model2,c("Shell.height..mm.","FLUPSY"))
plotmod2 <- ggpredict(model2,c("FLUPSY"))
plotmod3 <- ggpredict(model1,c("Shell.height..mm.","FLUPSY"))

plotwa <- ggplot(plotmod,aes(x,predicted,color=group)) +
  #scale_color_manual(values=wes_palette("GrandBudapest1")) + 
  scale_colour_brewer(palette = "BrBG") +
  geom_point(size=4) +
  geom_errorbar(data=plotmod, mapping=aes(x=x, ymin=conf.low, ymax=conf.high), width=0.1) +
  geom_line(aes(group=group)) +
  xlab("Shell height (mm)") +
  ylab(expression(paste("Predicted infestation"))) +
  ylim(0,0.5) +
  theme_classic() +
  theme(plot.title=element_text(size=16,hjust=0.5,face="plain"), axis.text.y=element_text(size=16), 
        axis.title.y=element_text(size=16), axis.text.x=element_text(size=16),legend.text = element_text(size=14), 
        axis.title.x=element_text(size=16), legend.title = element_text(size=14),panel.grid.minor=element_line(color=NA))
plotwa

plotall <- ggplot(plotmod3,aes(x,predicted,color=group)) +
  #scale_color_manual(values=wes_palette("GrandBudapest1")) + 
  scale_colour_brewer(palette = "BrBG") + #BrBG
  geom_point(size=4) +
  geom_errorbar(data=plotmod3, mapping=aes(x=x, ymin=conf.low, ymax=conf.high), width=0.1) +
  geom_line(aes(group=group)) +
  xlab("Shell height (cm)") +
  ylab(expression(paste("Predicted infestation"))) +
  #ylim(0,0.5) +
  theme_classic() +
  theme(plot.title=element_text(size=16,hjust=0.5,face="plain"), axis.text.y=element_text(size=16), 
        axis.title.y=element_text(size=16), axis.text.x=element_text(size=16),legend.text = element_text(size=14), 
        axis.title.x=element_text(size=16), legend.title = element_text(size=14),panel.grid.minor=element_line(color=NA))
plotall
