library(shiny)
library(haven)
library(mice)
library(ggplot2)
library(gganimate)

airq <- airquality
airq$Month <- format(ISOdate(2004,1:12,1),"%B")[airq$Month]

ggplot(airq, aes(Day, Temp, group = Month)) + 
  geom_line() + 
  geom_segment(aes(xend = 31, yend = Temp), linetype = 2, colour = 'grey') + 
  geom_point(size = 2) + 
  geom_text(aes(x = 31.1, label = Month), hjust = 0) + 
  transition_reveal(Day) + 
  coord_cartesian(clip = 'off') + 
  labs(title = 'Temperature in New York', y = 'Temperature (°F)') + 
  theme_minimal() + 
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))


niter <- 20
nimp <- 5

dat <- read_sav(file="Backpain50 MI missing.sav")
dat$Group <- factor(ifelse(is.na(dat$Tampascale), 1,0))
levels(dat$Group) <- c("Observed", "Imputed")
imp0 <- mice(data=dat, m=1, method="mean", maxit = 0)
pred <- imp0$predictorMatrix
pred[, c("ID", "Group", "Disability", "Radiation")] <- 0
imp <- mice(data=dat, predictorMatrix = pred, m=nimp, maxit = niter, seed=1232)

#tamp_imp1 <- imp$chainMean[3,,1]

imp_mean <- apply(imp$chainMean, 1L, c)[, "Tampascale"]
iternr <- rep(1:niter, nimp)
impnr <- rep(1:nimp, each=niter)
dat_imp <- data.frame(impnr, iternr, imp_mean)
dat_imp

#myX <- scale_x_continuous(name="Iteration number", limits=c(0,10))
#myY <- scale_y_continuous(name="Tampa scale", limits=c(35,45))

ggplot(dat_imp, aes(iternr, imp_mean, group = impnr)) +
  geom_line(aes(colour=rep(c("blue", "red", "purple", "orange", "black"), each=20))) +
  geom_segment(aes(xend = niter, yend = imp_mean), linetype = 2, colour = 'grey') + 
  geom_point(size = 2) + 
  geom_text(aes(x = niter+0.1, label = impnr), hjust = 0) + 
  transition_reveal(iternr) + 
  coord_cartesian(clip = 'off') + 
  labs(title = 'MI algorithm', y = 'Mean Imputed') + 
  theme(legend.position="none") #theme_minimal() #+ 
  #theme(plot.margin = margin(5.5, 40, 5.5, 5.5))


  