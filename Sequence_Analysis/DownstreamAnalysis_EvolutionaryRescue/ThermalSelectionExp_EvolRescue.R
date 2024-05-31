################################################################################
##### Thermal Selection Experiment Evolutionary Rescue Simulation  #############  
##################### Written by Lisa Couper  ##################################
##################### Last updated: current  ###################################

# The code below estimates adaptive potential 
# using the evolutionary rescue framework without plasticity
# It estimates the max rate of environmental change the population could tolerate
# incorporating heritabilty, rmax, generation time, strength of selection, and Vp

library(akima)
library(ggplot2)

#### Rescue model approach: Hold r_max constant, simluate over other paramters ####

# gamma (y): 0.578 across all rounds; 0.463, 0.606, 0.590 for round 1-3 individually
# LD-pruned h2 = 0.370, Vp: 0.278 -> h2Vp =  0.10286
# Unpruned: h2 = 0.228, Vp = 0.258 -> h2Vp = 0.058824
# average h2vp = 0.080842

gamma = seq(0.3931, 0.6931, by = 0.025)
h2 = seq(0.1, 0.33, by = 0.025)
Vp = seq(0.1, 0.33, by = 0.025)
h2Vp = h2*Vp
cmb = expand.grid(gamma,h2,Vp)
colnames(cmb) = c("gamma", "h2", "Vp")
gridh2vP <- cmb$h2 * cmb$Vp

obvs <- c(0.578, mean(c(0.370,0.228)), mean(c(0.278, 0.258)))
rescue1(obvs)
rescue2(obvs)
rescue3(obvs)

###### R_max = 0.15 ######

rescue1 = function(row) {
  rmax = 0.15
  gamma = as.numeric(row[1])
  h2Vp = as.numeric(row[2]) *as.numeric(row[3])
  GenT = 1
  N_c = sqrt((2 * rmax * gamma)/GenT) * h2Vp
  return(N_c)
}

Nc_rmax1 <- apply(cmb,1,rescue1)

df_rmax1 <- cbind.data.frame(gridh2vP, cmb$gamma, Nc_rmax1)
colnames(df_rmax1) <- c("h2Vp", "gamma", "envr_change")
my.df.interp_rmax1 <- interp(x = df_rmax1$h2Vp, y = df_rmax1$gamma, z = df_rmax1$envr_change, 
                             nx = 200, ny = 200, duplicate = "strip") # removes duplicate values
df_new_rmax1 <- as.data.frame(interp2xyz(my.df.interp_rmax1))
colnames(df_new_rmax1) <- c("h2Vp", "gamma", "envr_change")


# Plot 1: R_max = 0.15
col_breaks <- seq(0, 0.08, by = 0.005)
coloptions <- brewer.pal(11, "RdBu")[1:7]
mycolors <- rev(colorRampPalette(coloptions)(17))

p1 <- ggplot(df_new_rmax1, aes(x = h2Vp, y = gamma, fill = envr_change)) + 
  geom_tile() + theme_minimal() + 
  scale_fill_gradientn(limits = c(0, 0.08), colours = mycolors, breaks = col_breaks) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + 
 scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("heritability x phenotypic variance") + 
  ylab("strength of selection") + 
  geom_contour(aes(z = envr_change), breaks = 0.027, col = 'darkblue', lwd = 1) +
  geom_contour(aes(z = envr_change), breaks = 0.024, col = '#743089', lwd = 1) +
  geom_contour(aes(z = envr_change), breaks = 0.029, col = 'darkgreen', lwd = 1) +
  geom_point(aes(x = 0.080842, y = 0.578), pch = 16, col = "black", size = 4) + 
  geom_segment(aes(x=0.058824,xend=0.10286,y=0.578,yend=0.578)) + # vary across h2 estimates
  geom_segment(aes(x=0.080842,xend=0.080842,y=0.463,yend=0.606)) # vary across selection strength estimates

###### R_max = 0.25 #######

rescue2 = function(row) {
  rmax = 0.25
  gamma = as.numeric(row[1])
  h2Vp = as.numeric(row[2]) *as.numeric(row[3])
  GenT = 1
  N_c = sqrt((2 * rmax * gamma)/GenT) * h2Vp
  return(N_c)
}

Nc_rmax2 <- apply(cmb,1,rescue2)

df_rmax2 <- cbind.data.frame(gridh2vP, cmb$gamma, Nc_rmax2)
colnames(df_rmax2) <- c("h2Vp", "gamma", "envr_change")
my.df.interp_rmax2 <- interp(x = df_rmax2$h2Vp, y = df_rmax2$gamma, z = df_rmax2$envr_change, 
                             nx = 200, ny = 200, duplicate = "strip") # removes duplicate values
df_new_rmax2 <- as.data.frame(interp2xyz(my.df.interp_rmax2))
colnames(df_new_rmax2) <- c("h2Vp", "gamma", "envr_change")


# Plot 2: R_max = 0.25
col_breaks <- seq(0, 0.08, by = 0.005)
coloptions <- brewer.pal(11, "RdBu")[1:7]
mycolors <- rev(colorRampPalette(coloptions)(17))

p2 <- ggplot(df_new_rmax2, aes(x = h2Vp, y = gamma, fill = envr_change)) + 
  geom_tile() + theme_minimal() + 
  scale_fill_gradientn(limits = c(0, 0.08), colours = mycolors, breaks = col_breaks) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("heritability x phenotypic variance") + 
  ylab("strength of selection") + 
  geom_contour(aes(z = envr_change), breaks = 0.027, col = 'darkblue', lwd = 1) +
  geom_contour(aes(z = envr_change), breaks = 0.024, col = '#743089', lwd = 1) +
  geom_contour(aes(z = envr_change), breaks = 0.029, col = 'darkgreen', lwd = 1) +
  geom_point(aes(x = 0.080842, y = 0.578), pch = 16, col = "black", size = 4) + 
  geom_segment(aes(x=0.058824,xend=0.10286,y=0.578,yend=0.578)) + # vary across h2 estimates
  geom_segment(aes(x=0.080842,xend=0.080842,y=0.463,yend=0.606)) # vary across selection strength estimates


###### R_max = 0.35 #######

rescue3 = function(row) {
  rmax = 0.35
  gamma = as.numeric(row[1])
  h2Vp = as.numeric(row[2]) *as.numeric(row[3])
  GenT = 1
  N_c = sqrt((2 * rmax * gamma)/GenT) * h2Vp
  return(N_c)
}

Nc_rmax3 <- apply(cmb,1,rescue3)

df_rmax3 <- cbind.data.frame(gridh2vP, cmb$gamma, Nc_rmax3)
colnames(df_rmax3) <- c("h2Vp", "gamma", "envr_change")
my.df.interp_rmax3 <- interp(x = df_rmax3$h2Vp, y = df_rmax3$gamma, z = df_rmax3$envr_change, 
                             nx = 200, ny = 200, duplicate = "strip") # removes duplicate values
df_new_rmax3 <- as.data.frame(interp2xyz(my.df.interp_rmax3))
colnames(df_new_rmax3) <- c("h2Vp", "gamma", "envr_change")

# Plot 3: R_max = 0.35
col_breaks <- seq(0, 0.08, by = 0.005)
coloptions <- brewer.pal(11, "RdBu")[1:7]
mycolors <- rev(colorRampPalette(coloptions)(17))

p3 <- ggplot(df_new_rmax3, aes(x = h2Vp, y = gamma, fill = envr_change)) + 
  geom_tile() + theme_minimal() + 
  scale_fill_gradientn(limits = c(0, 0.08), colours = mycolors, breaks = col_breaks) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("heritability x phenotypic variance") + 
  ylab("strength of selection") + 
  geom_contour(aes(z = envr_change), breaks = 0.027, col = 'darkblue', lwd = 1) +
  geom_contour(aes(z = envr_change), breaks = 0.024, col = '#743089', lwd = 1) +
  geom_contour(aes(z = envr_change), breaks = 0.029, col = 'darkgreen', lwd = 1) +
  geom_point(aes(x = 0.080842, y = 0.578), pch = 16, col = "black", size = 4) + 
  geom_segment(aes(x=0.058824,xend=0.10286,y=0.578,yend=0.578)) + # vary across h2 estimates
  geom_segment(aes(x=0.080842,xend=0.080842,y=0.463,yend=0.606)) # vary across selection strength estimates

