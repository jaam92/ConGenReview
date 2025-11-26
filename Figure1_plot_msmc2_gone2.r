# Ethiopian Wolf MSMC2 Analysis and Visualization
# Analysis of demographic history using MSMC2 results
# Comparing different mutation rates (mu) and rho-to-mu ratios (r)
# Testing both 4 haplotypes (n=2 individuals) and 12 haplotypes (n=6 individuals)

library(tidyverse)
library(scales)
setwd("/scratch1/marjanak/wolf_msmc2/msmc2_chr1_rtest")

# Analysis parameters
gen <- 3  # Generation time in years
mu <- 5e-9  # Mutation rate 1: 5 × 10^-9 per site per generation
mu2 <- 1e-8  # Mutation rate 2: 1 × 10^-8 per site per generation  
mu3 <- 5e-8  # Mutation rate 3: 5 × 10^-8 per site per generation

# Load and process MSMC2 results for 6 individuals (12 haplotypes)
# Mutation rate 1 (mu = 5e-9)
msmc_0.1<-read_tsv("EW.chr1.r0.1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu),
         LeftYears = gen*(left_time_boundary/mu),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu))

msmc_1<-read_tsv("EW.chr1.r1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu),
         LeftYears = gen*(left_time_boundary/mu),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu))

msmc_10<-read_tsv("EW.chr1.r10.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu),
         LeftYears = gen*(left_time_boundary/mu),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu))

msmc <- bind_rows("r=0.1"=msmc_0.1,"r=1"=msmc_1,"r=10"=msmc_10,.id="r")

# Load and process MSMC2 results for 2 individuals (4 haplotypes)
# Mutation rate 1 (mu = 5e-9)
msmc_0.1_n2<-read_tsv("EW.n2.chr1.r0.1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu),
         LeftYears = gen*(left_time_boundary/mu),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu))

msmc_1_n2<-read_tsv("EW.n2.chr1.r1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu),
         LeftYears = gen*(left_time_boundary/mu),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu))

msmc_10_n2<-read_tsv("EW.n2.chr1.r10.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu),
         LeftYears = gen*(left_time_boundary/mu),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu))

msmc_n2 <- bind_rows("r=0.1"=msmc_0.1_n2,"r=1"=msmc_1_n2,"r=10"=msmc_10_n2,.id="r")

# Process results with mutation rate 2 (mu = 1e-8)
# 6 individuals (12 haplotypes)
msmc_0.1_2<-read_tsv("EW.chr1.r0.1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu2),
         LeftYears = gen*(left_time_boundary/mu2),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu2))

msmc_1_2<-read_tsv("EW.chr1.r1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu2),
         LeftYears = gen*(left_time_boundary/mu2),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu2))

msmc_10_2<-read_tsv("EW.chr1.r10.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu2),
         LeftYears = gen*(left_time_boundary/mu2),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu2))

msmc2 <- bind_rows("r=0.1"=msmc_0.1_2,"r=1"=msmc_1_2,"r=10"=msmc_10_2,.id="r")

# 2 individuals (4 haplotypes)
msmc_0.1_2_n2<-read_tsv("EW.n2.chr1.r0.1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu2),
         LeftYears = gen*(left_time_boundary/mu2),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu2))

msmc_1_2_n2<-read_tsv("EW.n2.chr1.r1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu2),
         LeftYears = gen*(left_time_boundary/mu2),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu2))

msmc_10_2_n2<-read_tsv("EW.n2.chr1.r10.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu2),
         LeftYears = gen*(left_time_boundary/mu2),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu2))

msmc2_n2 <- bind_rows("r=0.1"=msmc_0.1_2_n2,"r=1"=msmc_1_2_n2,"r=10"=msmc_10_2_n2,.id="r")

# Process results with mutation rate 3 (mu = 5e-8)
# 6 individuals (12 haplotypes)
msmc_0.1_3<-read_tsv("EW.chr1.r0.1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu3),
         LeftYears = gen*(left_time_boundary/mu3),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu3))

msmc_1_3<-read_tsv("EW.chr1.r1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu3),
         LeftYears = gen*(left_time_boundary/mu3),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu3))

msmc_10_3<-read_tsv("EW.chr1.r10.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu3),
         LeftYears = gen*(left_time_boundary/mu3),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu3))

msmc3 <- bind_rows("r=0.1"=msmc_0.1_3,"r=1"=msmc_1_3,"r=10"=msmc_10_3,.id="r")

# 2 individuals (4 haplotypes)
msmc_0.1_3_n2<-read_tsv("EW.n2.chr1.r0.1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu3),
         LeftYears = gen*(left_time_boundary/mu3),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu3))

msmc_1_3_n2<-read_tsv("EW.n2.chr1.r1.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu3),
         LeftYears = gen*(left_time_boundary/mu3),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu3))

msmc_10_3_n2<-read_tsv("EW.n2.chr1.r10.msmc2.final.txt")%>%
  mutate(Ne = (1/lambda)/(2*mu3),
         LeftYears = gen*(left_time_boundary/mu3),
         generationsLeft = LeftYears/gen,
         RightYears = gen*(right_time_boundary/mu3))

msmc3_n2 <- bind_rows("r=0.1"=msmc_0.1_3_n2,"r=1"=msmc_1_3_n2,"r=10"=msmc_10_3_n2,.id="r")

# Combine all mutation rate analyses
msmc_mut <- bind_rows("mu=5e-9"=msmc2, "mu=1e-8"=msmc,"mu=5e-8"=msmc3, .id="mu") 
msmc_mut_n2 <- bind_rows("mu=5e-9"=msmc2_n2, "mu=1e-8"=msmc_n2,"mu=5e-8"=msmc3_n2, .id="mu") 

# Combine all analyses (mutation rates, recombination rates, and sample sizes)
msmc_mut_n <- bind_rows("n=4"= msmc_mut_n2, "n=12"= msmc_mut, .id="size") %>% 
  mutate(mu = factor(mu, levels=c("mu=5e-9","mu=1e-8","mu=5e-8")))%>% 
  mutate(r = factor(r, levels=c("r=0.1","r=1","r=10")))

# Create demographic history plot
# Comparing effective population size over time across different parameters
pm<-ggplot()+
  geom_step(data=subset(msmc_mut_n,LeftYears>2000),aes(x=LeftYears, y=Ne,color=size),linewidth=0.7)+
  scale_y_log10(labels=comma)+
  scale_x_log10(labels=comma)+
  facet_grid(rows=vars(r),cols=vars(mu))+
  theme_bw()+
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=13),
        strip.background = element_blank(),
        plot.margin = margin(0.1,1,0.1,0.1, "cm"),
        panel.spacing.x = unit(0.5, "cm"),
        strip.text = element_text(size=13),
        legend.position="top")+  
  scale_color_manual(values=c("royalblue4","steelblue1"))+
  scale_alpha_manual(values=c(0.3,1,0.5))+
  labs(x = "Log Years", y =  "Log Effective Population Size")+
  guides(color=guide_legend(title="Haplotypes"))


# Analysis of demographic history using Gone2 results
# Comparing different recombination rates (R)
# Testing both 6 haplotypes (n=3 individuals) and 12 haplotypes (n=6 individuals)

# Load Gone2 results for 6 individuals (12 haplotypes)

r1 <- read_tsv("/scratch1/marjanak/wolf_msmc2/EW_chr1_r1_GONE2_Ne")
r1.5 <- read_tsv("/scratch1/marjanak/wolf_msmc2/EW_chr1_r1.5_GONE2_Ne")
r0.75 <- read_tsv("/scratch1/marjanak/wolf_msmc2/EW_chr1_r0.75_GONE2_Ne")

goner <- bind_rows("R=1"=r1,"R=1.5"=r1.5,"R=0.75"=r0.75,.id="Recomb") %>% 
  mutate(Years=Generation*gen)

# Load and process Gone2 results for 3 individuals (6 haplotypes)

r1_n3 <- read_tsv("/scratch1/marjanak/wolf_msmc2/EW_chr1_i3_r1_GONE2_Ne")
r1.5_n3 <- read_tsv("/scratch1/marjanak/wolf_msmc2/EW_chr1_i3_r1.5_GONE2_Ne")
r0.75_n3 <- read_tsv("/scratch1/marjanak/wolf_msmc2/EW_chr1_i3_r0.75_GONE2_Ne")

goner3 <- bind_rows("R=1"=r1_n3,"R=1.5"=r1.5_n3,"R=0.75"=r0.75_n3,.id="Recomb") %>% 
  mutate(Years=Generation*gen)

# Combine results

gone_r_n <- bind_rows("n=12"=goner,"n=6"=goner3,.id="size")


# Create demographic history plot
# Comparing effective population size over time across different parameters

pg<-ggplot()+
  geom_step(data=subset(gone_r_n,Years<100),aes(x=Years, y=Ne_diploids,color=size),linewidth=0.7)+
  scale_y_log10(labels=comma)+
  scale_x_log10(labels=comma)+
  facet_grid(cols=vars(Recomb))+
  theme_bw()+
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=13),
        strip.background = element_blank(),
        plot.margin = margin(0.1,1,0.1,0.1, "cm"),
        panel.spacing.x = unit(0.5, "cm"),
        strip.text = element_text(size=13),
        legend.position="top")+  
  scale_color_manual(values=c("royalblue4","steelblue1"))+
  labs(x = "Log Years", y =  "Log Effective Population Size")+
  guides(color=guide_legend(title="Haplotypes"))


# Plot Figure 1

ggarrange(pm,pg,nrow=2,heights = c(1,0.5), labels=c("A) Historic Demography","B) Recent Demography"),label.x = -0.1)



