##########################################################################
# TITLE: Using a coupled integral projection model to investigate 
# interspecific competition during an invasion: An application to silver 
# carp and gizzard shad
# AUTHORS: J. Peirce, G. Sandland, D. Schumann, H. Thompson,  R. Erickson
# CODED BY: J. Peirce
##########################################################################

# Remember to set the working directory:
# Session -> Set Working Directory -> To Source Location
library(tidyverse)
library(scales)

#####################
## Fish Parameters ##
#####################

# Gizzard shad (S) 
m_par_S <- tibble(
  # von-Bertalanaffy growth and standard deviation
  ## growth from Michaletz (2017) paper - Table 1
  grow_rate =  0.26 , # von-Bertalanaffy K
  grow_max  =  394.3, # von-Berralanaffy Linf 
  grow_sd   =   25,  # growth sd
  # adult survival
  surv_min  =  0.002, # Bodola (1955)
  # Then in mm: 1 - 8.804128*grow_rate^.73*grow_max^(-.33)
  surv_max  =  1 - 8.804128*grow_rate^(.73)*grow_max^(-.33), 
  # Estimated from LTRM - main channel
  surv_infl = 80.0136, # infection point 
  surv_slope = -139.9312, # slope 
  # effect of carp density on growth and reproduction
  g = 0, # original model - not going to change
  # reproduction bla
  prob_spawn = 0.9, # probability of spawning (in a year)
  # egg production - 3 parameter logistic
  ## From Jons and Miranda
  egg_max =  737.512105,
  egg_slope = -7.181562,
  egg_infl = 313.742703,
  ## From Bodola (1955):
  egg_viable = 0.002, # probability of egg -> age-0 fish
  # age-0 survival 
  ## Estimated from Michaletz (2009)
  surv0 = 0.2686, # surv_alpha - survival % without density effect
  f_S = 0.0030, # surv_beta - affect of age-0 gizzard density on survival
  # AND interspecific competition coeff:
  f_C = 0, # for original model
  # recruits
  ## New recruit from Michaletz (2017)
  recruit_mean = 105,
  recruit_sd = 25,
  # length-weight: log(W) = lw_alpha*log(z)+lw_beta 
  lw_alpha = 2.952,
  lw_beta = -4.904,
  # simulation setup
  L = 0, # smallest simulated length
  U = 600, # largest simulated length
  N = 400, # number of length bins
  delta_z = (U-L)/N
)

# Silver carp (C) 
m_par_C<- tibble(
  # von-Bertalanaffy growth and standard deviation
  grow_rate =   0.173, # Hayer et al. (2014)
  grow_max  =  1224, # Hayer et al. (2014)
  grow_sd   =   40,  # Williamson and Garvey (2005)
  # adult survival
  surv_min  =  0.0, 
  # Then in mm: 1 - 8.872*grow_rate^.73*grow_max^(-.33)
  surv_max = 1 - 8.804128*grow_rate^(.73)*grow_max^(-.33),
  surv_infl = 450, # inflection - Shireman et al. (1978) 
  surv_slope = -0.015,  
  # effect of carp density on growth and reproduction
  g = 5*10^(-9), # original model - not going to change
  # reproduction
  prob_spawn = 0.9, # probability of spawning (in a year)
  # egg production - 3 parameter logistic
  egg_max =  6636993,
  egg_slope = 0.002565,
  egg_infl = 11965.07,  
  egg_viable = 0.00005, # probability of egg -> age-0 fish
  # age-0 survival 
  ## no dependence on density
  surv0 = 1, # surv_alpha - survival % already build into egg_viable_C
  f_C = 0, # for original model - not going to change
  # AND interspecific competition coeff:
  f_S = 0, # 
  # recruits
  ## New recruit from Williamson and Garvey 2008 for t=1 (2017)
  recruit_mean = 319,
  recruit_sd = 40, # same as grow_sd
  # length-weight: log(W) = lw_alpha*log(z)+lw_beta
  lw_alpha = 3.122,
  lw_beta = -5.294,
  # simulation setup
  L = 0, # smallest simulated length
  U = 1200, # largest simulated length
  N = 400, # number of length bins
  delta_z = (U-L)/N
)

# compute zmesh (representative values from each length bin)
zmesh_S <- with(m_par_S, L + ((1:N) - 1 / 2) * (U - L) / N)
zmesh_C <- with(m_par_C, L + ((1:N) - 1 / 2) * (U - L) / N)

###################################
## Functions Common to Both Fish ##
###################################

## length-weight function
length_weight <- function(z, m_par){
  return(10^(m_par$lw_beta) * z^(m_par$lw_alpha))
}

## adult survival - 4 parameter logistic
s_z <- function(z, m_par) {
  m_par$surv_min + (m_par$surv_max - m_par$surv_min) /
    (1 + exp(m_par$surv_slope * ((z) - (m_par$surv_infl))))
}

## growth function 
# given you are size z now returns the pdf of size z1 next time.
G_z1z <- function(z1, z, m_par) {
  mu <- m_par$grow_max * (1 - exp(- m_par$grow_rate)) +
    exp(-m_par$grow_rate) * z           # mean size next year
  sig <- m_par$grow_sd                       # sd about mean
  new_length <- dnorm(z1, mean = mu, sd = sig)    # pdf that you are size z1
  # given you were size z
  return(new_length)
}

C_1z1 <- function(z1, m_par) {
  mu <- m_par$recruit_mean
  sig <- m_par$recruit_sd
  recruit <- dnorm(z1, mean = mu, sd = sig)
  recruit <- recruit / (sum(recruit * m_par$delta_z))
  return(recruit)
}

##############################
## Egg Production Functions ##
##############################

# gizzard shad
egg_S <- function(z,m_par){
  # Eggs produced (note: data are in thousands)
  return(1000 * m_par$egg_max / (1 + exp(m_par$egg_slope 
                                         * (log(z)-log(m_par$egg_infl)))))
}

# silver carp
egg_C <- function(z,m_par){
  weight = length_weight(z,m_par)
  return(m_par$egg_max/(1+exp(-1*m_par$egg_slope*(weight-m_par$egg_infl))))
}

###################
## Density Terms ##
###################

biomass <- function(z, n, m_par){
  weight_dist <- length_weight(z, m_par) * n
  return(sum(weight_dist * m_par$delta_z))
}

# age-0 density
d_0 <- function(z, n, egg, m_par){
  # distribution of VIABLE age0 from current population n
  age0dist <-  m_par$prob_spawn * egg(z, m_par) *
    m_par$egg_viable * n
  return(10^(-3) * (sum(age0dist * m_par$delta_z)))
}

########################
## No-density Kernels ##
########################

# growth and survival kernel - NO DENSITY
P_z1z <- function(z1, z, m_par) { # note: must have n = n_C
  G_matrix <- matrix(0, m_par$N, m_par$N)
  for (x in 1:m_par$N) {
    G_matrix[, x] <- G_z1z(z, rep(z[x], times = m_par$N), m_par)
    G_matrix[, x] <- G_matrix[, x] / (sum(G_matrix[, x]) * m_par$delta_z)
  }
  return(G_matrix %*% diag(s_z(z, m_par)))
}

# fecundity kernel - NO DENSITY
F_z1z <- function(z1, z, egg, m_par) { # note: must have n = n_S
  age1_per <- m_par$prob_spawn * egg(z, m_par) *
    m_par$egg_viable * m_par$surv0
  #returns fecundity kernel (as a matrix). Recruits= F.dot(n*delta_z)
  return(outer(C_1z1(z1, m_par), age1_per))
}

###########
## Model ##
###########

tf <- 60 # number of years
# initial gizzard shad length distibution
n_S0 <- rep(0, length(zmesh_S))
n_S0_total <- 3437.333
n_S0 <- dnorm(zmesh_S, mean = 0.5*m_par_S$grow_max, sd = 30)
#normal like LTRM 1994
n_S0 <- (n_S0 / sum(n_S0)) * n_S0_total / m_par_S$delta_z
# Note: sum(n[,1])*delta_z = n0_total

# initial silver carp length distribution
n_C0 <- rep(0, length(zmesh_C))
n_C0_total <- 816.381
n_C0 <- dnorm(zmesh_C, mean = 0.5*m_par_C$grow_max, sd = 30)
n_C0 <- (n_C0 / sum(n_C0)) * n_C0_total / m_par_C$delta_z

length_dist <- function(zS, zC, nS0, nC0){
  n_S <- matrix(0,length(zS), tf)
  n_S[,1] <- nS0
  n_C <- matrix(0,length(zC), tf)
  n_C[,1] <- nC0
  
for (i in 1:(tf-1)){
    # compute biomass and age-0 densities
    bio_C <- biomass(zmesh_C, n_C[,i], m_par_C)
    d0_S <- d_0(zmesh_S, n_S[,i], egg_S, m_par_S)
    d0_C <- d_0(zmesh_C, n_C[,i], egg_C, m_par_C)
    
    k_iter_S <- (P_z1z(zmesh_S, zmesh_S, m_par_S) * # growth and survival
                   exp(-1*m_par_S$g*bio_C) + # intraspecific comp on P
                   F_z1z(zmesh_S, zmesh_S, egg_S, m_par_S) * # fecundity
                   exp(-1*(m_par_S$f_S*d0_S + m_par_S$f_C*d0_C))) # comp on F   
    n_S[, i + 1] <- k_iter_S %*% n_S[, i] * m_par_S$delta_z
    ####
    k_iter_C <- (P_z1z(zmesh_C, zmesh_C, m_par_C) * # growth and survival
                   exp(-1*m_par_C$g*bio_C) + # intraspecific comp on P
                   F_z1z(zmesh_C, zmesh_C, egg_C, m_par_C) * # fecundity
                   exp(-1*(m_par_C$f_S*d0_S + m_par_C$f_C*d0_C))) # comp on F
    n_C[, i + 1] <- k_iter_C %*% n_C[, i] * m_par_C$delta_z
  }
  return(tibble(nS = n_S, nC = n_C))
}

model <- length_dist(zmesh_S, zmesh_C, n_S0, n_C0)

fish_totals <- tibble(
  t = 1:tf,
  shad = colSums(model$nS * m_par_S$delta_z),
  carp = colSums(model$nC * m_par_C$delta_z),
  C2S = rep(m_par_S$f_C, times = tf),
  S2C = rep(m_par_C$f_S, times = tf)
)

#########################################
## Graphs of Total Population vs. Time ##
#########################################

fish_totals_graph <- gather(fish_totals, "fish", "total", -t, -C2S, -S2C)
ggplot(data = fish_totals_graph) +
  geom_line(aes(x=t, y =total, color = fish)) +
  labs(x = "years",
       y = "total") +
  labs(x = "time (in years)",
       y = "total density",
       color = "Fish") +
  scale_y_continuous(label=comma) +
  scale_color_manual(values=c("firebrick4", "dodgerblue4"))+
  theme_bw() +
  theme(text = element_text(size=12),
        aspect.ratio = .7)

##################################
## And now with m_par$f vectors ##
##################################
C2S_vec <- 10^c(-4, -2, -1)
S2C_vec <- 0

fish_totals <- NULL
for (i in C2S_vec){
  for (j in S2C_vec){
    m_par_S$f_C <- i
    m_par_C$f_S <- j
    model <- length_dist(zmesh_S, zmesh_C, n_S0, n_C0)
    fish_totals <- bind_rows(fish_totals,
                             tibble( t = 1:tf,
                                     shad = colSums(model$nS * m_par_S$delta_z),
                                     carp = colSums(model$nC * m_par_C$delta_z),
                                     C2S = rep(signif(m_par_S$f_C, digits=3), times = tf),
                                     S2C = rep(signif(m_par_C$f_S, digits=3), times = tf)))
  }
}

fish_totals_graph <- gather(fish_totals, "fish", "total", -t, -C2S, -S2C)
fish_totals_graph$C2S <- factor(fish_totals_graph$C2S, levels = c(10^(-4), 10^(-2), 10^(-1)),
                  labels = c("low", "med", "high"))
ggplot(data = fish_totals_graph) +
  geom_line(aes(x=t, y =total, color = fish)) +
  labs(x = "years",
       y = "total") +
  labs(x = "time (in years)",
       y = "total density",
       color = "Fish") +
  facet_grid(~C2S)+
  scale_y_continuous(label=comma) +
  scale_color_manual(values=c("firebrick4", "dodgerblue4"))+
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(text = element_text(size=12),
        aspect.ratio = .7)
ggsave("carp_on_shad.png")

##################################
## And now with m_par$f vectors ##
##################################
C2S_vec <- 0
S2C_vec <- c(10^(-5), 10^(-4), 0.002)

fish_totals <- NULL
for (i in C2S_vec){
  for (j in S2C_vec){
    m_par_S$f_C <- i
    m_par_C$f_S <- j
    model <- length_dist(zmesh_S, zmesh_C, n_S0, n_C0)
    fish_totals <- bind_rows(fish_totals,
                             tibble( t = 1:tf,
                                     shad = colSums(model$nS * m_par_S$delta_z),
                                     carp = colSums(model$nC * m_par_C$delta_z),
                                     C2S = rep(signif(m_par_S$f_C, digits=3), times = tf),
                                     S2C = rep(signif(m_par_C$f_S, digits=3), times = tf)))
  }
}

fish_totals_graph <- gather(fish_totals, "fish", "total", -t, -C2S, -S2C)
fish_totals_graph$S2C <- factor(fish_totals_graph$S2C, levels = S2C_vec,
                                labels = c("low", "med", "high"))
ggplot(data = fish_totals_graph) +
  geom_line(aes(x=t, y =total, color = fish)) +
  labs(x = "years",
       y = "total") +
  labs(x = "time (in years)",
       y = "total density",
       color = "Fish") +
  facet_grid(~S2C) +
  scale_y_continuous(label=comma) +
  scale_color_manual(values=c("firebrick4", "dodgerblue4"))+
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(text = element_text(size=12),
        aspect.ratio = .7)
ggsave("shad_on_carp.png")

##########
## Both ##
##########
C2S_vec <- 10^c(-4, -2, -1)
S2C_vec <- c(10^(-5), 10^(-4),  0.002)

fish_totals <- NULL
for (i in C2S_vec){
  for (j in S2C_vec){
    m_par_S$f_C <- i
    m_par_C$f_S <- j
    model <- length_dist(zmesh_S, zmesh_C, n_S0, n_C0)
    fish_totals <- bind_rows(fish_totals,
                             tibble( t = 1:tf,
                                     shad = colSums(model$nS * m_par_S$delta_z),
                                     carp = colSums(model$nC * m_par_C$delta_z),
                                     C2S = rep(signif(m_par_S$f_C, digits=3), times = tf),
                                     S2C = rep(signif(m_par_C$f_S, digits=3), times = tf)))
  }
}

fish_totals_graph <- gather(fish_totals, "fish", "total", -t, -C2S, -S2C)
fish_totals_graph$C2S <- factor(fish_totals_graph$C2S, levels = C2S_vec,
                                 labels = c("carp on shad: low", 
                                            "carp on shad: med", 
                                            "carp on shad: high"))
fish_totals_graph$S2C <- factor(fish_totals_graph$S2C, levels = S2C_vec,
                                 labels = c("shad on carp: low", 
                                            "shad on carp: med", 
                                            "shad on carp: high"))
ggplot(data = fish_totals_graph) +
  geom_line(aes(x=t, y =total, color = fish)) +
  labs(x = "time (in years)",
       y = "total density",
       color = "Fish") +
  facet_grid(rows = vars(S2C), cols = vars(C2S))+
  scale_y_continuous(label=comma) +
  scale_color_manual(values=c("firebrick4", "dodgerblue4"))+
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(text = element_text(size=12),
        aspect.ratio = .8)
ggsave("carp_and_shad.png")


###############
## Threshold ##
###############
C2S_vec <- seq(from=10^(-4), to=10^(-1), length.out = 100)
S2C_vec <- seq(from=10^(-5), to=3*10^(-3), length.out = 100)

min_pop <-10000 # extinction threshold value
years2ave <- 10  

fish_state <- NULL
for (i in C2S_vec){ 
  for (j in S2C_vec){
    m_par_S$f_C <- i
    m_par_C$f_S <- j
    model <- length_dist(zmesh_S, zmesh_C, n_S0, n_C0)
    nS <- colSums(model$nS * m_par_S$delta_z)
    nC <- colSums(model$nC * m_par_C$delta_z)
    aveS <- mean(nS[(tf-years2ave):tf])
    aveC <- mean(nC[(tf-years2ave):tf])
  state = if(aveS>min_pop & aveC > min_pop){
     "both"
   } else if (aveS>min_pop & aveC <= min_pop){
     "shad only"
   } else if (aveS<=min_pop & aveC > min_pop){
     "carp only"
   } else {
     "none"
       }
    fish_state <- bind_rows(fish_state, tibble(
          C2S = signif(i, digits=3),
          S2C = signif(j, digits=3),
          aveS = signif(aveS, digits=3),
          aveC = aveC,
          state= state
          ))
  }
}

colors <- c("both" = "#56B4E9", 
            "carp only" = "#D55E00", 
            "shad only" =  "#009E73",
            "none" = "white")
del_x <- 1.5*(C2S_vec[2]-C2S_vec[1])
del_y <- 1.5*(S2C_vec[2]-S2C_vec[1])
ggplot(data = fish_state, 
       aes(x=C2S, y=S2C, fill= state)) + 
  geom_tile(width = del_x, height = del_y) +
  scale_fill_manual(values = colors) +
  labs(x = expression(~f[21]~silver~carp~effect~on~gizzard~shad))+
  labs(y = expression(~f[12]~gizzard~effect~on~silver~carp))+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme(text = element_text(size=18),
        aspect.ratio = .7)
ggsave("bifurcation.png")
