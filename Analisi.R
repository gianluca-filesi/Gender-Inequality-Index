rm(list = ls()) 
graphics.off()

library(BAS)
library(car)
library(corrplot)
library(dplyr)
library(ellipse)
library(faraway)
library(fastDummies)
library(GGally)
library(leaps)
library(MASS)
library(Matrix)
library(misc3d)
library(plot3D)
library(qpcR)
library(RColorBrewer)
library(readxl)
library(rgl)

verifica_res <- function(g){
  
  par(mfrow = c(2,1))
  plot(g$fitted.values, g$residuals, xlab = "Fitted", ylab = "Residuals", 
       main = "Residuals vs Fitted Values", pch = 16 )
  abline( h = 0, lwd = 2, lty = 2, col = 'red' ) #Sembrano omoscehdastici
  qqnorm(g$residuals, ylab = 'Residuals', pch=16) # C'è sicuramente qualche outlier
  qqline(g$residuals) 
  
  S=shapiro.test(g$residuals)
  print(S$p)
  par(mfrow = c(1,1))
}

DF<-read.csv("Gender_Inequality_Index.csv",row.names = 1)[-3]
colSums(is.na(DF))
colnames(DF) <- c("HDI","GII","Mat_mortality","Adolescent_birth_rate","Seats_parliament","F_edu","M_edu","F_lab","M_lab")
str(DF)



DF$HDI=as.factor(DF$HDI)
#ggpairs(data = df[,-1], title ="Relationships between predictors & response",lower = list(continuous=wrap("points", alpha = 0.5, size=0.1)))

## 80% of the sample size
smp_size <- floor(0.80 * nrow(DF))

## set the seed to make your partition reproducible
set.seed(400)
train_ind <- sample(seq_len(nrow(DF)), size = smp_size)

df <- DF[train_ind, ]
test <- DF[-train_ind, ]

#correlation

X = df[,-1]


## Modello 1
g1 = lm(GII ~ ., data = df)
summary(g1) 
verifica_res(g1) 

# Nonlinearity 
d <-lm(GII ~ .-Seats_parliament, data = df)$res
m <-lm(Seats_parliament ~ .-GII, data = df)$res
plot(m,d,xlab="Seats_parliament",ylab="GII", main="Partial Regression")
abline(0,g1$coef['Seats_parliament'])


## Modello 2 - vif selection
vif(g1)
prplot(g1,6)

g2_1 = update(g1, .~. -M_edu)
summary(g2_1) # 0.8963
verifica_res(g2_1) # 0.002030099
vif(g2_1)

g2_2 = update(g1, .~. -F_edu)
summary(g2_2) # 0.8937
verifica_res(g2_2) # 0.003802853

g2 = update(g1, .~. -F_edu-M_edu)
summary(g2) # 0.8806
verifica_res(g2) # 0.004346862

g2=g2_2
vif(g2)
## Modello 3 - step and R2adj 
x = model.matrix(g2) [ , -1 ]
y = df$GII

adjr = leaps( x, y, method = "adjr2" )
bestmodel_adjr2_ind = which.max( adjr$adjr2 )
adjr$which[ bestmodel_adjr2_ind, ] 
maxadjr( adjr, 5 ) # Non ne vale la pena

step(g2_1, direction = "both" , trace = T) #È il modello g2_2


# Modello 3 - leverages and Cook's distance
{lev = hatvalues(g2)
  r = g2$rank  
  p = r-1
  n = dim(df)[1] 
  plot(g2$fitted.values, lev, ylab = "Leverages", main = "Plot of Leverages", 
       pch = 16, col = 'black' )
  abline( h = 2 * r/n, lty = 2, col = 'red' )
  watchout_points_lev = lev[ which( lev > 2 * r/n  ) ]
  watchout_ids_lev = seq_along( lev )[ which( lev > 2 * r/n ) ]
  points( g2$fitted.values[ watchout_ids_lev ], watchout_points_lev, col = 'red', pch = 16 )
  lev [ lev >  2 * r/n ]
  sum(lev[lev >  2 * r/n ] )
} #Leverages
g3_1 = lm(GII ~ .-F_edu, data = df, subset = (lev < 2*r/n) )
summary(g3_1) # 0.8993
verifica_res(g3_1) #0.001180344

# Standardized 
{
  gs = summary(g2)
  res_std = g2$res/gs$sigma
  watchout_ids_rstd = which( abs( res_std ) > 2 )
  watchout_rstd = res_std[ watchout_ids_rstd ]
  watchout_rstd
  
  # Plot Residui standardizzati 
  
  plot( g2$fitted.values, res_std, ylab = "Standardized Residuals", main = "Standardized Residuals" )
  abline( h = c(-2,2), lty = 2, col = 'orange' )
  points( g2$fitted.values[watchout_ids_rstd], 
          res_std[watchout_ids_rstd], col = 'red', pch = 16 )
  points( g2$fitted.values[watchout_ids_lev], 
          res_std[watchout_ids_lev], col = 'orange', pch = 16 )
  legend('topleft', col = c('red','orange'), 
         c('Standardized Residuals', 'Leverages'), pch = rep( 16, 2 ), bty = 'n' ) } #Standardization

# Studentized
{stud = rstandard(g2)
  watchout_ids_stud = which( abs( stud ) > 2 )
  watchout_stud = stud[ watchout_ids_stud ]
  plot( g2$fitted.values, stud, ylab = "Studentized Residuals", main = "Studentized Residuals", pch = 16 )
  points( g2$fitted.values[watchout_ids_stud], 
          stud[watchout_ids_stud], col = 'pink', pch = 16 )
  points( g2$fitted.values[watchout_ids_lev], 
          stud[watchout_ids_lev], col = 'orange', pch = 16 )
  abline( h = c(-2,2), lty = 2, col = 'orange' )
  legend('topleft', col = c('pink','orange'), 
         c('Studentized Residual', 'Leverages'), pch = rep( 16, 3 ), bty = 'n' )} #Studentization

## Modello 3 - Cook's distance
# Cook's distance
{Cdist = cooks.distance(g2)
  
  watchout_ids_Cdist = which( Cdist > 4/(n-p-1) ) 
  watchout_Cdist = Cdist[ watchout_ids_Cdist ]
  
  
  {par( mfrow = c( 1, 3 ) )
    plot( g2$fitted.values, Cdist, pch = 16, xlab = 'Fitted values', 
          ylab = 'Cooks Distance', main = 'Cooks Distance' )
    points( g2$fitted.values[ watchout_ids_Cdist ], Cdist[ watchout_ids_Cdist ], 
            col = 'green', pch = 16 )
    plot( g2$fitted.values, stud, pch = 16, xlab = 'Fitted values', 
          ylab = 'Studentized Residuals', main = 'Studentized Residuals' )
    points( g2$fitted.values[ watchout_ids_stud ], stud[ watchout_ids_stud ], 
            col = 'pink', pch = 16 )
    plot( g2$fitted.values, lev, pch = 16, xlab = 'Fitted values', 
          ylab = 'Leverages', main = 'Leverages' )
    points( g2$fitted.values[ watchout_ids_lev ], lev[ watchout_ids_lev ],
            col = 'orange', pch = 16 )
    
    par( mfrow = c( 1, 1 ) )} #plot
  
  id_to_keep = !( 1:n %in% watchout_ids_Cdist )} #Cook's distance
g3 = lm(GII ~ .-F_edu, data = df[ id_to_keep, ])
summary(g3) # 0.9324
verifica_res(g3) #0.1669743

## Modello 4 - Cook's distance and levergaes
g4 = lm(GII ~ .-F_edu, data = df[ id_to_keep, ],subset = (lev < 2*p/n) )
summary(g4) # 0.9181
verifica_res(g4) #0.1442702

plot(g1, which = 5)

# Attualmente il migliore è g3

# Modello 5 - box cox
b1 = boxcox(GII ~ .-F_edu, data = df[ id_to_keep, ])
best_lambda_ind1 = which.max( b1$y )
best_lambda1 = b1$x[ best_lambda_ind1] #1.39394

g5 = lm((GII^best_lambda1-1)/best_lambda1 ~ .-F_edu, data = df[ id_to_keep, ]) 
summary(g5) #0.9369
verifica_res(g5) # 0.7450989

# Modello 6 - box cox con valore decente 1.4
best_lambda2 = 1.5
g6 = lm((GII^best_lambda2-1)/best_lambda2 ~ .-F_edu, data = df[ id_to_keep, ]) 
summary(g6) #0.9371
verifica_res(g6) # 0.7545415

(summary(g6)$adj-summary(g3)$adj)/summary(g3)$adj*100 # Migliora di solo l'1.9%

# Sum up
summary(g1)$r.squared
summary(g1)$adj
verifica_res(g1) 
summary(g2)$r.squared
summary(g2)$adj
verifica_res(g2) 
summary(g3_1)$r.squared
summary(g3_1)$adj
verifica_res(g3_1) 
summary(g3)$r.squared
summary(g3)$adj
verifica_res(g3) 
summary(g4)$r.squared
summary(g4)$adj
verifica_res(g4) 
summary(g5)$r.squared
summary(g5)$adj
verifica_res(g5) 
summary(g6)$r.squared
summary(g6)$adj
verifica_res(g6) 

# Uso g3

#Prediction 

y.conf = predict(g3,test,interval = "confidence", se = T )
y.pred = predict(g3,test,interval = "prediction", se = T )

azzurro = brewer.pal(3,'Paired')[1] 
blu = brewer.pal(3,'Paired')[2] 
bluscuro = brewer.pal(9,'PuBu')[9] 

matplot( 1:dim(test)[1], y.conf$fit, lty = c( 1, 2, 2 ), col = c(blu,azzurro,azzurro), type = "l",
         xlab = "ing", ylab = "GII", main = 'Intervallo di Confidenza' )
points( 1:dim(test)[1], test$GII, col = bluscuro, pch = 1)
matplot( (1:dim(test)[1]), y.pred$fit, lty = c( 1, 2, 2 ), col = c(blu,azzurro,azzurro), type = "l",
         xlab = "ing", ylab = "GII", main = 'Intervallo di previsione' )
points( 1:dim(test)[1], test$GII, col = bluscuro, pch = 1)


## ===================Grafici per presentazione=====================

colnames(DF) <- c("HDI","GII","Mortalità al parto","Adolescenti in gravidanza","Seggi parlamentari","Educazione Femminile","Educazione Maschile","Donne lavoratrici","Uomini lavoratori")
colnames(df) <- c("HDI","GII","Mortalità al parto","Adolescenti in gravidanza","Seggi parlamentari","Educazione Femminile","Educazione Maschile","Donne lavoratrici","Uomini lavoratori")
colnames(test) <- c("HDI","GII","Mortalità al parto","Adolescenti in gravidanza","Seggi parlamentari","Educazione Femminile","Educazione Maschile","Donne lavoratrici","Uomini lavoratori")

# Colori
lilla = brewer.pal(9,'Paired')[9] 
viola = brewer.pal(10,'Paired')[10] 
azzurro = brewer.pal(3,'Paired')[1] 
blu = brewer.pal(3,'Paired')[2]
cref = '#675D76'

# Pie chart by HDI
png(filename = "NationsByHDI.png", width = 1200, height = 1200, pointsize = 30,units = "px", res = NA)
pie(table(DF$HDI), col = c(lilla,cref,lilla,cref), main = 'Nazioni per Indice di Sviluppo umano' )
dev.off()

# GGpairs iniziale
png(filename = "GGpairs1.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
ggpairs(data = df[,-1], title ="Relazione tra predditori e risposta",
        lower = list(continuous=wrap("points", alpha = 0.5, size=0.1)))
dev.off()

# Correlazione tra covariate
png(filename = "CorrelationNum.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
corrplot(cor(X), method = 'number')
dev.off()
png(filename = "CorrelationCol.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
corrplot(cor(X), method='color')
dev.off()
png(filename = "Heatmap.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
heatmap( cor( X ), Rowv = NA, Colv = NA, symm = TRUE, keep.dendro = F)
dev.off()
png(filename = "Correlation.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
image( as.matrix( cor( X ) ), main = 'Correlazione di X' )
dev.off()


#GGpairs dopo modello
png(filename = "GGpairs2.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
ggpairs(X)
dev.off()

# Nonlinearity 
png(filename = "Partial regression.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
plot(m,d,xlab="Seggi parlamentari",ylab="GII", main="Regressione parziale")
abline(0,g1$coef['Seats_parliament'])
dev.off()

#Prplot
png(filename = "Prplot.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
prplot(g1,6)
dev.off()

#Leverages
png(filename = "Leverages.png", width = 2052, height = 1155, pointsize = 50,units = "px", res = NA)
plot(g2$fitted.values, lev, xlab = "Valori fittati", ylab = "Leverages", main = "Grafico dei Leverages", 
     pch = 1, col = 'black' )
abline( h = 2 * r/n, lty = 2, col = cref)
points( g2$fitted.values[ watchout_ids_lev ], watchout_points_lev, col = lilla, pch = 16 )
dev.off()

#Residui standardizzati
png(filename = "Standard_Lev.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
plot( g2$fitted.values, res_std, xlab = "Valori fittati", ylab = "Residui standardizzati", main = "Residui strandardizzati" )
abline( h = c(-2,2), lty = 2, col = cref )
points( g2$fitted.values[watchout_ids_rstd], 
        res_std[watchout_ids_rstd], col = lilla, pch = 16 )
points( g2$fitted.values[watchout_ids_lev], 
        res_std[watchout_ids_lev], col = cref, pch = 16 )
legend('topleft', col = c(lilla,cref), 
       c('Residui standardizzati', 'Leverages'), pch = rep( 16, 2 ), bty = 'n' ) 
dev.off()

#Residui studentizzati
png(filename = "Stud_Lev.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
plot( g2$fitted.values, stud, xlab = "Valori fittati", ylab = "Residui Studentizzati", main = "Residui Studentizzati", pch = 1 )
points( g2$fitted.values[watchout_ids_stud], 
        stud[watchout_ids_stud], col = lilla, pch = 16 )
points( g2$fitted.values[watchout_ids_lev], 
        stud[watchout_ids_lev], col = cref, pch = 16 )
abline( h = c(-2,2), lty = 2, col = cref)
legend('topleft', col = c(lilla,cref), 
       c('Residui Studentizzati', 'Leverages'), pch = rep( 16, 3 ), bty = 'n' )
dev.off()

# Cook's distance 
png(filename = "Cook.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
plot( g2$fitted.values, Cdist, pch = 1, xlab = "Valori fittati",
      ylab = 'Distanza di Cook', main = 'Distanza di Cook' )
points( g2$fitted.values[ watchout_ids_Cdist ], Cdist[ watchout_ids_Cdist ], 
        col = lilla, pch = 16 )
abline( h = 4/(n-p-1), lty = 2, col = cref)
dev.off()

# Comparison 
png(filename = "CSS.png", width = 2052, height = 1155, pointsize = 50,units = "px", res = NA)
par( mfrow = c( 1, 3 ) )
plot( g2$fitted.values, Cdist, pch = 1, xlab = "Valori fittati",
      ylab = 'Distanza di Cook', main = 'Cooks Distance' )
points( g2$fitted.values[ watchout_ids_Cdist ], Cdist[ watchout_ids_Cdist ], 
        col = lilla, pch = 16 )
abline( h = 4/(n-p-1), lty = 2, col = cref)
plot( g2$fitted.values, stud, pch = 1, xlab = "Valori fittati",
      ylab = 'Residui Studentizzati', main = 'Residui Studentizzati' )
points( g2$fitted.values[ watchout_ids_stud ], stud[ watchout_ids_stud ], 
        col = lilla, pch = 16 )
abline( h = c(-2,2), lty = 2, col = cref )
plot( g2$fitted.values, res_std, xlab = "Valori fittati", ylab = "Residui Standardizzati", main = "Residui Strandardizzati" )
abline( h = c(-2,2), lty = 2, col = cref )
points( g2$fitted.values[watchout_ids_rstd], 
        res_std[watchout_ids_rstd], col = lilla, pch = 16 )
par( mfrow = c( 1, 1 ) )
dev.off()

#Boxcox 1
colnames(df) <- c("HDI","GII","Mat_mortality","Adolescent_birth_rate","Seats_parliament","F_edu","M_edu","F_lab","M_lab")
png(filename = "BoxCox.png", width = 1200, height = 1200, pointsize = 40,units = "px", res = NA)
b1 = boxcox(GII ~ .-F_edu, data = df[ id_to_keep, ])
axis(1, at = best_lambda1, labels = 1.23, col = 'black', pos = -653)
dev.off()
colnames(df) <- c("HDI","GII","Mortalità al parto","Adolescenti in gravidanza","Seggi parlamentari","Educazione Femminile","Educazione Maschile","Donne lavoratrici","Uomini lavoratori")


#Prediction and confidence
png(filename = "Prediction.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
matplot( (1:dim(test)[1]), y.pred$fit, lty = c( 1, 2, 2 ), col = c(cref,cref,cref), type = "l", xlab = "Ranking",
        ylab = "GII", main = 'Intervallo di previsione' )
points( 1:dim(test)[1], test$GII, col = c(cref,cref), pch = 19)
dev.off()
png(filename = "Confidence.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
matplot( 1:dim(test)[1], y.conf$fit, lty = c( 1, 2, 2 ), col = c(cref,cref,cref), type = "l", xlab = "Ranking",
        ylab = "GII", main = 'Intervallo di Confidenza' )
points( 1:dim(test)[1], test$GII, col = cref, pch = 19)
dev.off()


summary(g1)$adj #0.9931913
verifica_res(g1) #9.626529e-12
summary(g2)$adj #0.8805752
verifica_res(g2) #0.004346862
summary(g3_1)$adj #0.8841983
verifica_res(g3_1) #0.001180344
summary(g3)$adj #0.9198484
verifica_res(g3) #0.1669743
summary(g4)$adj #0.0.9181057
verifica_res(g4) #0.1442702
summary(g5)$adj #0.9369103
verifica_res(g5) #0.7450989
summary(g6)$adj #0.9370818
verifica_res(g6) #0.7545415
summary(g3)

png(filename = "QQnorm_1.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
par(mfrow = c(2,1))
plot(g1$fitted.values, g1$residuals, xlab = "Fitted", ylab = "Residuals", 
     main = "Residuals vs Fitted Values", pch = 16 )
abline( h = 0, lwd = 2, lty = 2, col = 'red' ) #Sembrano omoscehdastici
qqnorm(g1$residuals, ylab = 'Residuals', pch=16) # C'è sicuramente qualche outlier
qqline(g1$residuals) 

par(mfrow = c(1,1))
dev.off()

png(filename = "QQnorm_2.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
par(mfrow = c(2,1))
plot(g2$fitted.values, g2$residuals, xlab = "Fitted", ylab = "Residuals", 
     main = "Residuals vs Fitted Values", pch = 16 )
abline( h = 0, lwd = 2, lty = 2, col = 'red' ) #Sembrano omoscehdastici
qqnorm(g2$residuals, ylab = 'Residuals', pch=16) # C'è sicuramente qualche outlier
qqline(g2$residuals) 
par(mfrow = c(1,1))
dev.off()

png(filename = "QQnorm_3.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
par(mfrow = c(2,1))
plot(g3_1$fitted.values, g3_1$residuals, xlab = "Fitted", ylab = "Residuals", 
     main = "Residuals vs Fitted Values", pch = 16 )
abline( h = 0, lwd = 2, lty = 2, col = 'red' ) #Sembrano omoscehdastici
qqnorm(g3_1$residuals, ylab = 'Residuals', pch=16) # C'è sicuramente qualche outlier
qqline(g3_1$residuals) 
par(mfrow = c(1,1))
dev.off()

png(filename = "QQnorm_4.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
par(mfrow = c(2,1))
plot(g3$fitted.values, g3$residuals, xlab = "Fitted", ylab = "Residuals", 
     main = "Residuals vs Fitted Values", pch = 16 )
abline( h = 0, lwd = 2, lty = 2, col = 'red' ) #Sembrano omoscehdastici
qqnorm(g3$residuals, ylab = 'Residuals', pch=16) # C'è sicuramente qualche outlier
qqline(g3$residuals) 
par(mfrow = c(1,1))
dev.off()

png(filename = "QQnorm_5.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
par(mfrow = c(2,1))
plot(g4$fitted.values, g4$residuals, xlab = "Fitted", ylab = "Residuals", 
     main = "Residuals vs Fitted Values", pch = 16 )
abline( h = 0, lwd = 2, lty = 2, col = 'red' ) #Sembrano omoscehdastici
qqnorm(g4$residuals, ylab = 'Residuals', pch=16) # C'è sicuramente qualche outlier
qqline(g4$residuals) 
par(mfrow = c(1,1))
dev.off()

png(filename = "QQnorm_6.png", width = 1200, height = 1200, pointsize = 50,units = "px", res = NA)
par(mfrow = c(2,1))
plot(g6$fitted.values, g6$residuals, xlab = "Fitted", ylab = "Residuals", 
     main = "Residuals vs Fitted Values", pch = 16 )
abline( h = 0, lwd = 2, lty = 2, col = 'red' ) #Sembrano omoscehdastici
qqnorm(g6$residuals, ylab = 'Residuals', pch=16) # C'è sicuramente qualche outlier
qqline(g6$residuals) 
par(mfrow = c(1,1))
dev.off()

(summary(g6)$adj-summary(g3)$adj)/summary(g3)$adj*100 # Migliora di solo l'1.2%
boxplot( g3$res, main = "Boxplot dei residui", pch = 16, col = cref )
hist( g3$res, 10, probability = TRUE, col = cref, main = 'Residui'  )

