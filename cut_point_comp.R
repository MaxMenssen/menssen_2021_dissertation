
# Sample observations from the binomial proportion

set.seed(234)
y100 <- rbinom(n=100, 50, 0.5)
y100

set.seed(1345)
y20 <- rbinom(n=20, 50, 0.5)
y20

set.seed(3456)
y5 <- rbinom(n=5, 50, 0.5)
y5

set.seed(68)
fut <- rbinom(n=1, 50, 0.5)
fut

obs <- c(y100, y20, y5, fut)

xaxis <- c(rep("y100", 100),
               rep("y20", 20),
               rep("y5", 5),
           "fut")

dat <- data.frame(xaxis, obs)               
str(dat)


#-------------------------------------------------------------------------------


library(tidyverse)

# informal comparison

A <- ggplot(dat, aes(x=xaxis, y=obs))+
  theme_bw()+
  geom_jitter(data=dat[-126,], width=0.1, height=0, alpha=0.2)+
  geom_point(data=dat[126,], color="#7CAE00", size=3)+
  scale_x_discrete(name ="", 
                   limits=c("y100",
                            "y20", 
                            "y5",
                            "fut"),
                   labels=c("h=100",
                            "h=20",
                            "h=5",
                            "future obs."))+
  ylab("")+
  ggtitle("A: Informal comparison")


#-------------------------------------------------------------------------------

# Range


dat_range <- dat %>% 
  group_by(xaxis) %>% 
  summarise(min=min(obs), max=max(obs))

B <- ggplot(dat, aes(x=xaxis, y=obs))+
  theme_bw()+
  geom_jitter(data=dat[-126,], width=0.1, height=0, alpha=0.2)+
  geom_point(data=dat[126,], color="#7CAE00", size=3)+
  geom_linerange(aes(x="y100", ymin=14, ymax=34), color="red", size=1.1)+
  geom_linerange(aes(x="y20", ymin=19, ymax=30), color="red", size=1.1)+
  geom_linerange(aes(x="y5", ymin=24, ymax=29), color="red", size=1.1)+
  scale_x_discrete(name ="", 
                   limits=c("y100",
                            "y20", 
                            "y5",
                            "fut"),
                   labels=c("h=100",
                            "h=20",
                            "h=5",
                            "future obs."))+
  ylab("")+
  ggtitle("B: Range")

#-------------------------------------------------------------------------------
# Quantiles

dat_quant <- dat %>% 
  group_by(xaxis) %>%
  summarize(lower=quantile(obs, 0.25),
            upper=quantile(obs, 0.75))

C <- ggplot(dat, aes(x=xaxis, y=obs))+
  theme_bw()+
  geom_boxplot(data=dat[-126,], alpha=0.2)+
  geom_jitter(data=dat[-126,], width=0.1, height=0, alpha=0.3)+
  geom_point(data=dat[126,], color="#7CAE00", size=3)+
  geom_hline(aes(yintercept=qbinom(0.75,50, 0.5)), color="blue", lty="dashed")+
  geom_hline(aes(yintercept=qbinom(0.25,50, 0.5)), color="blue", lty="dashed")+
  
  geom_linerange(aes(x="y100", ymin=23, ymax=27), color="red", size=1.1)+
  geom_linerange(aes(x="y20", ymin=22.8, ymax=25.2), color="red", size=1.1)+
  geom_linerange(aes(x="y5", ymin=25, ymax=28), color="red", size=1.1)+
  
  scale_x_discrete(name ="", 
                   limits=c("y100",
                            "y20", 
                            "y5",
                            "fut"),
                   labels=c("h=100",
                            "h=20",
                            "h=5",
                            "future obs."))+
  ylab("")+
  ggtitle("C: Quantiles")


#-------------------------------------------------------------------------------


# CI

dat_ci <- dat %>% 
  group_by(xaxis) %>%
  summarize(mean = mean(obs),
            sd=sqrt(var(obs)),
            n=n(),
            lower=mean - 1.96 * sd/sqrt(n),
            upper=mean + 1.96 * sd/sqrt(n))
  
dat_ci


D <- ggplot(dat, aes(x=xaxis, y=obs))+
  theme_bw()+
  geom_jitter(data=dat[-126,], width=0.1, height=0, alpha=0.2)+
  geom_point(data=dat[126,], color="#7CAE00", size=3)+
  geom_hline(aes(yintercept=25), color="blue", lty="dashed")+
  
  geom_linerange(aes(x="y100", ymin=24.6, ymax=25.8), color="red", size=1.1)+
  geom_linerange(aes(x="y20", ymin=23, ymax=25.5), color="red", size=1.1)+
  geom_linerange(aes(x="y5", ymin=24.9, ymax=28.7), color="red", size=1.1)+
  
  scale_x_discrete(name ="", 
                   limits=c("y100",
                            "y20", 
                            "y5",
                            "fut"),
                   labels=c("h=100",
                            "h=20",
                            "h=5",
                            "future obs."))+
  ylab("")+
  ggtitle("D: Confidence intervals")

#-------------------------------------------------------------------------------

# mean +- 2 sd

dat_2sd <- dat %>% 
  group_by(xaxis) %>%
  summarize(mean = mean(obs),
            sd=sqrt(var(obs)),
            lower=mean - 2*sd,
            upper=mean + 2*sd)



E <- ggplot(dat, aes(x=xaxis, y=obs))+
  theme_bw()+
  geom_jitter(data=dat[-126,], width=0.1, height=0, alpha=0.2)+
  geom_point(data=dat[126,], color="#7CAE00", size=3)+
  # geom_hline(aes(yintercept=25), color="red", lty="dashed")+
  
  geom_linerange(aes(x="y100", ymin=18.9, ymax=31.5), color="red", size=1.1)+
  geom_linerange(aes(x="y20", ymin=18.3, ymax=30.2), color="red", size=1.1)+
  geom_linerange(aes(x="y5", ymin=22.5, ymax=31.1), color="red", size=1.1)+
  
  scale_x_discrete(name ="", 
                   limits=c("y100",
                            "y20", 
                            "y5",
                            "fut"),
                   labels=c("h=100",
                            "h=20",
                            "h=5",
                            "future obs."))+
  ylab("")+
  ggtitle("E: Mean plus minus 2 SD")


#-------------------------------------------------------------------------------

library(predint)
?beta_bin_pi

daty100 <- data.frame(x=y100, y=50-y100)
pi_y100 <- beta_bin_pi(histdat=daty100, newsize=50)
pi_y100
# total hist_prob quant_calib  pred_se    lower    upper
#    50    0.5038    1.922148 3.553935 18.35881 32.02119


daty20 <- data.frame(x=y20, y=50-y20)
pi_y20 <- beta_bin_pi(histdat=daty20, newsize=50)
pi_y20
# total hist_prob quant_calib  pred_se    lower    upper
#    50     0.485    1.863613 3.622101 17.49981 31.00019

daty5 <- data.frame(x=y5, y=50-y5)
pi_y5 <- beta_bin_pi(histdat=daty5, newsize=50)
pi_y5
# total hist_prob quant_calib  pred_se    lower    upper
#    50     0.536    1.844102 3.863878 19.67462 33.92538


FG <- ggplot(dat, aes(x=xaxis, y=obs))+
  theme_bw()+
  geom_jitter(data=dat[-126,], width=0.1, height=0, alpha=0.2)+
  geom_point(data=dat[126,], color="#7CAE00", size=3)+
  # geom_hline(aes(yintercept=25), color="red", lty="dashed")+
  
  geom_linerange(aes(x="y100", ymin=18.35881, ymax=32.02119), color="red", size=1.1)+
  geom_linerange(aes(x="y20", ymin=17.49981, ymax=31.00019), color="red", size=1.1)+
  geom_linerange(aes(x="y5", ymin=19.67462, ymax=33.92538), color="red", size=1.1)+
  
  scale_x_discrete(name ="", 
                   limits=c("y100",
                            "y20", 
                            "y5",
                            "fut"),
                   labels=c("h=100",
                            "h=20",
                            "h=5",
                            "future obs."))+
  ylab("")+
  ggtitle("F: Prediction intervals")



#-------------------------------------------------------------------------------

library(ggpubr)

ggarrange(A, B, C, D, E, FG,
          ncol = 2,
          nrow = 3)


ggsave("cut_point_comp.png", width=18, height=26, units="cm")

