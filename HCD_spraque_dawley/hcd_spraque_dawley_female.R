#-------------------------------------------------------------------------------
#------------------ HCD for female Harlan Spraque Dawley rats ------------------
#-------------------------------------------------------------------------------

# The HCD was obtained from the NTP historical control data reports 2020

#-------------------------------------------------------------------------------

library(tidyverse)

# Data import
hcd <- read.csv2("HCD_spraque_dawley\\hcd_spraque_dawley_female.csv") %>% 
  
  mutate(
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
  
  scale_color_manual(values=c("#7CAE00", "#E69F00", "#56B4E9"),
                     name  ="Cluster size (n)")+
  scale_x_discrete(name ="", 
                   limits=c("HCD:p_total",
                            "HCD:p_mammary", 
                            "Seralini:p_mammary"),
                   labels=c("HCD total",
                            "HCD mammary",
                            "Seralini"))+
  geom_text(aes(x=3.27, y=0.302, label="Control"), 
            size=4.5)+
  geom_text(aes(x=3.27, y=0.502, label="Treat. min"),
            size=4.5)+
  geom_text(aes(x=3.27, y=0.802, label="Treat. max"), 
            size=4.5)+
  
  theme(axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold",
                                   size=14),
        axis.title.x = element_text(size=14, 
                                    face="bold"),
        axis.title.y = element_text(size=14, 
                                    face="bold"),
        legend.title = element_text(size=10, 
                                    face="bold"))

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

###

# Quasi binomial glm
fit_qb <- glm(cbind(Tumors_Total, No_Tumors_Total) ~ 1,
              family=quasibinomial(),
              data=total)

summary(fit_qb)$dispersion
# 3.029673

###

pi_rho_est(total)
#     pi_hat    rho_hat 
# 0.86734694 0.03787178 


#####

# Quasibinomial assumption
qb_pi_total <- quasi_bin_pi(histdat=total, 
                            newsize=rep(10,10))


# PI for binomial proportion
qb_pi_total$lower_prop <- qb_pi_total$lower/10
qb_pi_total$upper_prop <- qb_pi_total$upper/10
qb_pi_total

#    total hist_prob quant_calib  pred_se lower upper lower_prop upper_prop
# 1     10 0.8673469    3.385527 2.712742     0    10          0          1
# 2     10 0.8673469    3.385527 2.712742     0    10          0          1
# 3     10 0.8673469    3.385527 2.712742     0    10          0          1
# 4     10 0.8673469    3.385527 2.712742     0    10          0          1
# 5     10 0.8673469    3.385527 2.712742     0    10          0          1
# 6     10 0.8673469    3.385527 2.712742     0    10          0          1
# 7     10 0.8673469    3.385527 2.712742     0    10          0          1
# 8     10 0.8673469    3.385527 2.712742     0    10          0          1
# 9     10 0.8673469    3.385527 2.712742     0    10          0          1
# 10    10 0.8673469    3.385527 2.712742     0    10          0          1

#####

# Betabinomial assumption
bb_pi_total <- beta_bin_pi(histdat=total, 
                           newsize=rep(10,10))

bb_pi_total$lower_prop <- bb_pi_total$lower/10
bb_pi_total$upper_prop <- bb_pi_total$upper/10
bb_pi_total$Source_Type <- "HCD:p_total"
bb_pi_total

#    total hist_prob quant_calib  pred_se    lower upper lower_prop upper_prop Source_Type
# 1     10 0.8673469    3.678203 1.254676 4.058517    10  0.4058517          1 HCD:p_total
# 2     10 0.8673469    3.678203 1.254676 4.058517    10  0.4058517          1 HCD:p_total
# 3     10 0.8673469    3.678203 1.254676 4.058517    10  0.4058517          1 HCD:p_total
# 4     10 0.8673469    3.678203 1.254676 4.058517    10  0.4058517          1 HCD:p_total
# 5     10 0.8673469    3.678203 1.254676 4.058517    10  0.4058517          1 HCD:p_total
# 6     10 0.8673469    3.678203 1.254676 4.058517    10  0.4058517          1 HCD:p_total
# 7     10 0.8673469    3.678203 1.254676 4.058517    10  0.4058517          1 HCD:p_total
# 8     10 0.8673469    3.678203 1.254676 4.058517    10  0.4058517          1 HCD:p_total
# 9     10 0.8673469    3.678203 1.254676 4.058517    10  0.4058517          1 HCD:p_total
# 10    10 0.8673469    3.678203 1.254676 4.058517    10  0.4058517          1 HCD:p_total


#-------------------------------------------------------------------------------

# Mammary gland
FFCA <-  mutate(hcd,
                No_FFCA=Total-FFCA) %>% 
  select(FFCA,
         No_FFCA)

FFCA <- FFCA[1:9,]

#### 

# Quasi binomial glm
fit_qb_ffca <- glm(cbind(FFCA, No_FFCA) ~ 1,
              family=quasibinomial(),
              data=FFCA)

summary(fit_qb_ffca)$dispersion
# 1.685541

###

pi_rho_est(FFCA)
#     pi_hat    rho_hat 
# 0.64081633 0.01280304 



# Quasibinomial assumption
qb_pi_FFCA <- quasi_bin_pi(histdat=FFCA, 
                           newsize=rep(10,10))


qb_pi_FFCA$lower_prop  <- qb_pi_FFCA$lower/10
qb_pi_FFCA$upper_prop  <- qb_pi_FFCA$upper/10
qb_pi_FFCA$Source_Type <- "HCD:p_mammary"
qb_pi_FFCA

#    total hist_prob quant_calib  pred_se     lower upper lower_prop upper_prop   Source_Type
# 1     10 0.6408163    2.097754 2.861872 0.4046594    10 0.04046594          1 HCD:p_mammary
# 2     10 0.6408163    2.097754 2.861872 0.4046594    10 0.04046594          1 HCD:p_mammary
# 3     10 0.6408163    2.097754 2.861872 0.4046594    10 0.04046594          1 HCD:p_mammary
# 4     10 0.6408163    2.097754 2.861872 0.4046594    10 0.04046594          1 HCD:p_mammary
# 5     10 0.6408163    2.097754 2.861872 0.4046594    10 0.04046594          1 HCD:p_mammary
# 6     10 0.6408163    2.097754 2.861872 0.4046594    10 0.04046594          1 HCD:p_mammary
# 7     10 0.6408163    2.097754 2.861872 0.4046594    10 0.04046594          1 HCD:p_mammary
# 8     10 0.6408163    2.097754 2.861872 0.4046594    10 0.04046594          1 HCD:p_mammary
# 9     10 0.6408163    2.097754 2.861872 0.4046594    10 0.04046594          1 HCD:p_mammary
# 10    10 0.6408163    2.097754 2.861872 0.4046594    10 0.04046594          1 HCD:p_mammary

####



# Betabinomial assumption
bb_pi_FFCA <- beta_bin_pi(histdat=FFCA, 
                          newsize=rep(10,10))

bb_pi_FFCA$lower_prop  <- bb_pi_FFCA$lower/10
bb_pi_FFCA$upper_prop  <- bb_pi_FFCA$upper/10
bb_pi_FFCA$Source_Type <- "HCD:p_mammary"
bb_pi_FFCA

#    total hist_prob quant_calib  pred_se    lower upper lower_prop upper_prop   Source_Type
# 1     10 0.6408163    2.731885 1.618429 1.986801    10  0.1986801          1 HCD:p_mammary
# 2     10 0.6408163    2.731885 1.618429 1.986801    10  0.1986801          1 HCD:p_mammary
# 3     10 0.6408163    2.731885 1.618429 1.986801    10  0.1986801          1 HCD:p_mammary
# 4     10 0.6408163    2.731885 1.618429 1.986801    10  0.1986801          1 HCD:p_mammary
# 5     10 0.6408163    2.731885 1.618429 1.986801    10  0.1986801          1 HCD:p_mammary
# 6     10 0.6408163    2.731885 1.618429 1.986801    10  0.1986801          1 HCD:p_mammary
# 7     10 0.6408163    2.731885 1.618429 1.986801    10  0.1986801          1 HCD:p_mammary
# 8     10 0.6408163    2.731885 1.618429 1.986801    10  0.1986801          1 HCD:p_mammary
# 9     10 0.6408163    2.731885 1.618429 1.986801    10  0.1986801          1 HCD:p_mammary
# 10    10 0.6408163    2.731885 1.618429 1.986801    10  0.1986801          1 HCD:p_mammary

#-------------------------------------------------------------------------------



gg_pi_10 <-  ggplot(hcd, aes(x=Source_Type , y=Proportion)) + 
  
  theme_bw()+
  
  geom_linerange(aes(x=1.1,
                     ymin=qb_pi_total$lower_prop[1],
                     ymax=qb_pi_total$upper_prop[1]),
                 lty="dashed")+
  geom_linerange(aes(x=0.9,
                     ymin=bb_pi_total$lower_prop[1],
                     ymax=bb_pi_total$upper_prop[1]))+

  geom_linerange(aes(x=2.1,
                     ymin=qb_pi_FFCA$lower_prop[1],
                     ymax=qb_pi_FFCA$upper_prop[1]),
                 lty="dashed")+
  geom_linerange(aes(x=1.9,
                     ymin=bb_pi_FFCA$lower_prop[1],
                     ymax=bb_pi_FFCA$upper_prop[1]))+
  
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
  scale_color_manual(values=c("#7CAE00", "#E69F00", "#56B4E9"),
                     name  ="Cluster size")+
  
  scale_x_discrete(name ="", 
                   limits=c("HCD:p_total",
                            "HCD:p_mammary", 
                            "Seralini:p_mammary"),
                   labels=c("HCD total",
                            "HCD mammary",
                            "Seralini"))+
  geom_text(aes(x=3.27, y=0.302, 
                label="Control", ),
            size=3.5)+
  geom_text(aes(x=3.27, y=0.502, 
                label="Treat. min"), 
            size=3.5)+
  geom_text(aes(x=3.27, y=0.802, 
                label="Treat. max"), 
            size=3.5)+
  
  ggtitle("A: PI for future cluster size 10")+
  
  theme(axis.text.x = element_text(face="bold", 
                                   size=12),
        axis.text.y = element_text(face="bold",
                                   size=12),
        axis.title.x = element_text(size=12, 
                                    face="bold"),
        axis.title.y = element_text(size=12, 
                                    face="bold"),
        legend.title = element_text(size=10, 
                                    face="bold"),
        plot.title = element_text(size=12, 
                                  face="bold"))




# Beta binomial PI: Solid line
# Quasi binomial PI: Dashed line
gg_pi_10
ggsave("HCD_Seralini_prediction_intervals_m3.png", width=26, height=18, units="cm")



#-------------------------------------------------------------------------------
#------------------------ PI for n_fut=50 --------------------------------------
#-------------------------------------------------------------------------------

# total numbers of tumors

# Quasibinomial assumption
qb_pi_total_50 <- quasi_bin_pi(histdat=total, 
                               newsize=rep(50,10))


qb_pi_total_50$lower_pi <- qb_pi_total_50$lower/50
qb_pi_total_50$upper_pi <- qb_pi_total_50$upper/50
qb_pi_total_50

#   total hist_prob quant_calib  pred_se    lower upper  lower_pi upper_pi
# 1     50 0.8673469    1.746543 10.68914 24.69831    50 0.4939661        1
# 2     50 0.8673469    1.746543 10.68914 24.69831    50 0.4939661        1
# 3     50 0.8673469    1.746543 10.68914 24.69831    50 0.4939661        1
# 4     50 0.8673469    1.746543 10.68914 24.69831    50 0.4939661        1
# 5     50 0.8673469    1.746543 10.68914 24.69831    50 0.4939661        1
# 6     50 0.8673469    1.746543 10.68914 24.69831    50 0.4939661        1
# 7     50 0.8673469    1.746543 10.68914 24.69831    50 0.4939661        1
# 8     50 0.8673469    1.746543 10.68914 24.69831    50 0.4939661        1
# 9     50 0.8673469    1.746543 10.68914 24.69831    50 0.4939661        1
# 10    50 0.8673469    1.746543 10.68914 24.69831    50 0.4939661        1

###

# Betabinomial assumption
bb_pi_total_50 <- beta_bin_pi(histdat=total, 
                              newsize=rep(50,10))

bb_pi_total_50$lower_pi <- bb_pi_total_50$lower/50
bb_pi_total_50$upper_pi <- bb_pi_total_50$upper/50
bb_pi_total_50$Source_Type <- "HCD:p_total"
bb_pi_total_50

#    total hist_prob quant_calib pred_se    lower upper lower_pi upper_pi Source_Type
# 1     50 0.8673469    4.263555 4.25497 25.22605    50 0.504521        1 HCD:p_total
# 2     50 0.8673469    4.263555 4.25497 25.22605    50 0.504521        1 HCD:p_total
# 3     50 0.8673469    4.263555 4.25497 25.22605    50 0.504521        1 HCD:p_total
# 4     50 0.8673469    4.263555 4.25497 25.22605    50 0.504521        1 HCD:p_total
# 5     50 0.8673469    4.263555 4.25497 25.22605    50 0.504521        1 HCD:p_total
# 6     50 0.8673469    4.263555 4.25497 25.22605    50 0.504521        1 HCD:p_total
# 7     50 0.8673469    4.263555 4.25497 25.22605    50 0.504521        1 HCD:p_total
# 8     50 0.8673469    4.263555 4.25497 25.22605    50 0.504521        1 HCD:p_total
# 9     50 0.8673469    4.263555 4.25497 25.22605    50 0.504521        1 HCD:p_total
# 10    50 0.8673469    4.263555 4.25497 25.22605    50 0.504521        1 HCD:p_total

#-------------------------------------------------------------------------------
# mammary gland

# Quasibinomial assumption
qb_pi_FFCA_50 <- quasi_bin_pi(histdat=FFCA, 
                              newsize=rep(50,10))


qb_pi_FFCA_50$lower_pi <- qb_pi_FFCA_50$lower/50
qb_pi_FFCA_50$upper_pi <- qb_pi_FFCA_50$upper/50
qb_pi_FFCA_50$Source_Type <- "HCD:p_mammary"
qb_pi_FFCA_50


#    total hist_prob quant_calib  pred_se    lower   upper  lower_pi upper_pi   Source_Type
# 1     50 0.6408163    1.297773 11.27676 17.40613 46.6755 0.3481227  0.93351 HCD:p_mammary
# 2     50 0.6408163    1.297773 11.27676 17.40613 46.6755 0.3481227  0.93351 HCD:p_mammary
# 3     50 0.6408163    1.297773 11.27676 17.40613 46.6755 0.3481227  0.93351 HCD:p_mammary
# 4     50 0.6408163    1.297773 11.27676 17.40613 46.6755 0.3481227  0.93351 HCD:p_mammary
# 5     50 0.6408163    1.297773 11.27676 17.40613 46.6755 0.3481227  0.93351 HCD:p_mammary
# 6     50 0.6408163    1.297773 11.27676 17.40613 46.6755 0.3481227  0.93351 HCD:p_mammary
# 7     50 0.6408163    1.297773 11.27676 17.40613 46.6755 0.3481227  0.93351 HCD:p_mammary
# 8     50 0.6408163    1.297773 11.27676 17.40613 46.6755 0.3481227  0.93351 HCD:p_mammary
# 9     50 0.6408163    1.297773 11.27676 17.40613 46.6755 0.3481227  0.93351 HCD:p_mammary
# 10    50 0.6408163    1.297773 11.27676 17.40613 46.6755 0.3481227  0.93351 HCD:p_mammary


#####


# Betabinomial assumption
bb_pi_FFCA_50 <- beta_bin_pi(histdat=FFCA, 
                             newsize=rep(50,10))

bb_pi_FFCA_50$lower_pi <- bb_pi_FFCA_50$lower/50
bb_pi_FFCA_50$upper_pi <- bb_pi_FFCA_50$upper/50
bb_pi_FFCA_50$Source_Type <- "HCD:p_mammary"
bb_pi_FFCA_50

#    total hist_prob quant_calib  pred_se    lower    upper  lower_pi  upper_pi   Source_Type
# 1     50 0.6408163    3.170898 4.543068 17.63521 46.44642 0.3527042 0.9289285 HCD:p_mammary
# 2     50 0.6408163    3.170898 4.543068 17.63521 46.44642 0.3527042 0.9289285 HCD:p_mammary
# 3     50 0.6408163    3.170898 4.543068 17.63521 46.44642 0.3527042 0.9289285 HCD:p_mammary
# 4     50 0.6408163    3.170898 4.543068 17.63521 46.44642 0.3527042 0.9289285 HCD:p_mammary
# 5     50 0.6408163    3.170898 4.543068 17.63521 46.44642 0.3527042 0.9289285 HCD:p_mammary
# 6     50 0.6408163    3.170898 4.543068 17.63521 46.44642 0.3527042 0.9289285 HCD:p_mammary
# 7     50 0.6408163    3.170898 4.543068 17.63521 46.44642 0.3527042 0.9289285 HCD:p_mammary
# 8     50 0.6408163    3.170898 4.543068 17.63521 46.44642 0.3527042 0.9289285 HCD:p_mammary
# 9     50 0.6408163    3.170898 4.543068 17.63521 46.44642 0.3527042 0.9289285 HCD:p_mammary
# 10    50 0.6408163    3.170898 4.543068 17.63521 46.44642 0.3527042 0.9289285 HCD:p_mammary


#-------------------------------------------------------------------------------

gg_pi_50 <- ggplot(hcd, aes(x=Source_Type , y=Proportion)) + 
  
  theme_bw()+
  
  geom_linerange(aes(x=1.1,
                     ymin=qb_pi_total_50$lower_pi[1],
                     ymax=qb_pi_total_50$upper_pi[1]),
                 lty="dashed")+
  geom_linerange(aes(x=0.9,
                     ymin=bb_pi_total_50$lower_pi[1],
                     ymax=bb_pi_total_50$upper_pi[1]))+
  
  geom_linerange(aes(x=2.1,
                     ymin=qb_pi_FFCA_50$lower_pi[1],
                     ymax=qb_pi_FFCA_50$upper_pi[1]),
                 lty="dashed")+
  geom_linerange(aes(x=1.9,
                     ymin=bb_pi_FFCA_50$lower_pi[1],
                     ymax=bb_pi_FFCA_50$upper_pi[1]))+
  
  geom_point(data=filter(hcd,
                         Source_Type=="HCD:p_total" | Source_Type=="HCD:p_mammary"),
             aes(color=factor(Total)), size=3,
             position=position_jitter(height=0, width=0.1))+
  
  ylab("")+
  xlab("")+
  
  ylim(c(0,1))+
  scale_color_manual(values=c("#E69F00", "#56B4E9"),
                     name  ="Cluster size")+
  scale_x_discrete(name ="", 
                   limits=c("HCD:p_total",
                            "HCD:p_mammary"),
                   labels=c("HCD total",
                            "HCD mammary"))+
  
  ggtitle("B: PI for future cluster size 50")+
  
  theme(axis.text.x = element_text(face="bold", 
                                   size=12),
        axis.text.y = element_text(face="bold",
                                   size=12),
        axis.title.x = element_text(size=12, 
                                    face="bold"),
        axis.title.y = element_text(size=12, 
                                    face="bold"),
        legend.title = element_text(size=10, 
                                    face="bold"),
        plot.title = element_text(size=12, 
                                  face="bold"))



gg_pi_50 

#-------------------------------------------------------------------------------

library(ggpubr)

ggarrange(gg_pi_10, gg_pi_50,
          align = "h",
          widths = c(3, 1.8),
          common.legend = TRUE,
          legend="right")

ggsave("HCD_Seralini_prediction_intervals.png", width=26, height=18, units="cm")






