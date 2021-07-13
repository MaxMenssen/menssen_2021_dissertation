#-------------------------------------------------------------------------------
#------------------ HCD for female Harlan Spraque Dawley rats ------------------
#-------------------------------------------------------------------------------

# The HCD was obtained from the NTP historical control data reports 2020

#-------------------------------------------------------------------------------

library(tidyverse)

# Data import
hcd <- read.csv2("HCD_spraque_dawley\\hcd_spraque_dawley_female.csv") %>% 
  
  mutate(# Study_number=factor(Study_number),
    # Start=factor(Start),
    # Pathway=factor(Pathway),
    
    p_total=as.double(Tumors_Total/Total),
    p_mammary=as.double(FFCA/Total)) %>% 
  
  select(!c(Lab,
            Study_number,
            Start,
            Pathway)) %>% 
  
  pivot_longer(cols = starts_with("p"),
               names_to = "Prop_Type",
               # names_prefix = "p",
               values_to = "Proportion") %>% 
  
  arrange(Prop_Type)


hcd <- data.frame(hcd[-(22:24),])
str(hcd)
hcd$Source_Type <- factor(factor(hcd$Source) : factor(hcd$Prop_Type),
                          levels=c("HCD:p_total", "HCD:p_mammary", "Seralini:p_mammary"))


#-------------------------------------------------------------------------------

# Graphic without PI

ggplot(hcd, aes(x=Source_Type , y=Proportion)) + 
  
  theme_bw()+
  
  geom_point(data=filter(hcd,
                         Source_Type=="HCD:p_total" | Source_Type=="HCD:p_mammary"),
             aes(color=factor(Total)), size=3,
             position=position_jitter(height=0, width=0.1))+
  
  geom_point(data=filter(hcd,
                         Source_Type=="Seralini:p_mammary"),
             aes(color=factor(Total)), size=3)+
  
  ylab("Tumor rate (y/n)")+
  xlab("")+
  
  scale_colour_discrete(name  ="Numbers of rats\nper group (n)")+
  scale_x_discrete(name ="", 
                   limits=c("HCD:p_total",
                            "HCD:p_mammary", 
                            "Seralini:p_mammary"),
                   labels=c("HCD total",
                            "HCD mammary",
                            "Seralini"))+
  geom_text(aes(x=3.27, y=0.302, label="Control"), size=4.5)+
  geom_text(aes(x=3.27, y=0.502, label="Treat. min"), alpha=0.05, size=4.5)+
  geom_text(aes(x=3.27, y=0.802, label="Treat. max"), alpha=0.05, size=4.5)

ggsave("HCD_Seralini.png", width=26, height=18, units="cm")




#-------------------------------------------------------------------------------
#------------------------- PI for the HCD --------------------------------------
#-------------------------------------------------------------------------------

library(predint)

# Total numbers of tumors
total <- mutate(hcd,
                No_Tumors_Total=Total-Tumors_Total) %>% 
  select(Tumors_Total,
         No_Tumors_Total)

total <- total[1:9,]


# Quasibinomial assumption
qb_pi_total <- quasi_bin_pi(histdat=total, 
                            newsize=c(10, 10, 10))


qb_pi_total$lower_pi <- qb_pi_total$lower/10
qb_pi_total$upper_pi <- qb_pi_total$upper/10
qb_pi_total

# Betabinomial assumption
bb_pi_total <- beta_bin_pi(histdat=total, 
                           newsize=c(10, 10, 10))

bb_pi_total$lower_pi <- bb_pi_total$lower/10
bb_pi_total$upper_pi <- bb_pi_total$upper/10
bb_pi_total$Source_Type <- "HCD:p_total"

#-------------------------------------------------------------------------------

# Mammary gland
FFCA <-  mutate(hcd,
                No_FFCA=Total-FFCA) %>% 
  select(FFCA,
         No_FFCA)

FFCA <- FFCA[1:9,]

# Quasibinomial assumption
qb_pi_FFCA <- quasi_bin_pi(histdat=FFCA, 
                           newsize=c(10, 10, 10))


qb_pi_FFCA$lower_pi <- qb_pi_FFCA$lower/10
qb_pi_FFCA$upper_pi <- qb_pi_FFCA$upper/10
qb_pi_FFCA$Source_Type <- "HCD:p_mammary"

# Betabinomial assumption
bb_pi_FFCA <- beta_bin_pi(histdat=FFCA, 
                          newsize=c(10, 10, 10))

bb_pi_FFCA$lower_pi <- bb_pi_FFCA$lower/10
bb_pi_FFCA$upper_pi <- bb_pi_FFCA$upper/10
bb_pi_FFCA$Source_Type <- "HCD:p_mammary"

#-------------------------------------------------------------------------------



ggplot(hcd, aes(x=Source_Type , y=Proportion)) + 
  
  theme_bw()+
  
  geom_linerange(aes(x=1.1,
                     ymin=qb_pi_total$lower_pi[1],
                     ymax=qb_pi_total$upper_pi[1]),
                 lty="dashed")+
  geom_linerange(aes(x=0.9,
                     ymin=bb_pi_total$lower_pi[1],
                     ymax=bb_pi_total$upper_pi[1]))+
  
  geom_linerange(aes(x=2.1,
                     ymin=qb_pi_FFCA$lower_pi[1],
                     ymax=qb_pi_FFCA$upper_pi[1]),
                 lty="dashed")+
  geom_linerange(aes(x=1.9,
                     ymin=bb_pi_FFCA$lower_pi[1],
                     ymax=bb_pi_FFCA$upper_pi[1]))+
  
  geom_point(data=filter(hcd,
                         Source_Type=="HCD:p_total" | Source_Type=="HCD:p_mammary"),
             aes(color=factor(Total)), size=3,
             position=position_jitter(height=0, width=0.1))+
  
  geom_point(data=filter(hcd,
                         Source_Type=="Seralini:p_mammary"),
             aes(color=factor(Total)), size=3)+
  
  
  
  ylab("Tumor rate (y/n)")+
  xlab("")+
  
  ylim(c(0,1))+
  scale_colour_discrete(name  ="Numbers of rats\nper group (n)")+
  scale_x_discrete(name ="", 
                   limits=c("HCD:p_total",
                            "HCD:p_mammary", 
                            "Seralini:p_mammary"),
                   labels=c("HCD total",
                            "HCD mammary",
                            "Seralini"))+
  geom_text(aes(x=3.27, y=0.302, label="Control"), size=4.5)+
  geom_text(aes(x=3.27, y=0.502, label="Treat. min"), alpha=0.05, size=4.5)+
  geom_text(aes(x=3.27, y=0.802, label="Treat. max"), alpha=0.05, size=4.5)


# Beta binomial PI: Solid line
# Quasi binomial PI: Dashed line
ggsave("HCD_Seralini_prediction_intervals_m3.png", width=26, height=18, units="cm")








