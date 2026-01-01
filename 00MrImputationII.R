###################
thetahat <- function(y, x, iter = 5, convplot = TRUE){
df2 <- data.frame(y, x)
n1 <- nrow(df2)
y2 <- df2$y
x1 <- df2$x
df4 <- na.omit(df2)

#OLS: Starting Values
model0 <- lm(y ~ x -1, data = df4)
ehat1 <- residuals(model0)
ehat2 <- ehat1^2

#Estimating Theta: First Stage
x12 <- df4$x^2
model3b <- lm(log(ehat2) ~ log(x12))
that0c <- summary(model3b)$coef[2, 1]
that1c <- summary(model3b)$coef[2, 1]

#Estimating Theta: Iterations
b6hatj <- NULL
that1d <- NULL

for(j in 1:iter){
b1 <- sum(df4$x^(1-2*that1c)*df4$y)/sum(df4$x^(2*(1-that1c)))
ehat1b <- (df4$y - b1*df4$x)
ehat2b <- ehat1b^2
model3b <- lm(log(ehat2b) ~ log(x12))
that1c <- summary(model3b)$coef[2, 1]
that1d[j] <- summary(model3b)$coef[2, 1]
}
that1e <- c(that0c, that1d)
iteration1 <- c(1:(iter+1))

if(convplot == TRUE){
plot(iteration1, that1e, type = "l", xlab = "iteration", ylab = "theta", main = "Convergence")
}
estimated.theta <- data.frame(round(that1c, 4))
rownames(estimated.theta) <- "output"
colnames(estimated.theta) <- "theta.hat"
estimated.theta
}

###################
BPTtest <- function(y, x, theta = 0){
df2 <- data.frame(y, x)
y2 <- df2$y
x1 <- df2$x
df3 <- na.omit(df2)
n1 <- nrow(df3)
b1 <- sum(df3$x^(1-2*theta)*df3$y)/sum(df3$x^(2*(1-theta)))
u1 <- (df3$y - b1*df3$x)/(df3$x^theta)
u2 <- u1^2
logu2 <- log(u1^2)
model1 <- lm(logu2 ~ df3$x)
BPT <- n1*summary(model1)$r.squared
df <- 1
p.value <- pchisq(BPT, 1, lower.tail = FALSE)
output <- round(data.frame(BPT, df, p.value), 4)
rownames(output) <- "output"
output
}

###################
mrimpute2 <- function(y, x, M = 5, iter = 5){

df2 <- data.frame(y, x)
n1 <- nrow(df2)
y2 <- df2$y
x1 <- df2$x

dfimp <- matrix(NA, n1, M + 2)
df3 <- na.omit(df2)

#Bootstrap sample framework
sampleframe <- matrix(sample(nrow(df3), nrow(df2)*M, replace = TRUE),
                      nrow = nrow(df2), ncol = M) 

#Multiple Imputation
for(m in 1:M){

#Bootstrap
df4 <- df3[sampleframe[,m],]

#OLS: Starting Values
model0 <- lm(y ~ x -1, data=df4)
ehat1 <- residuals(model0)
ehat2 <- ehat1^2

#Estimating Theta: First Stage
x12 <- df4$x^2
model3b <- lm(log(ehat2) ~ log(x12))
that0b <- summary(model3b)$coef[2, 1]
that1c <- summary(model3b)$coef[2, 1]

#Estimating Theta: Iterations
b6hatj <- NULL
that1d <- NULL

for(j in 1:iter){
b1 <- sum(df4$x^(1-2*that1c)*df4$y)/sum(df4$x^(2*(1-that1c)))
ehat1b <- (df4$y - b1*df4$x)
ehat2b <- ehat1b^2
model3b <- lm(log(ehat2b) ~ log(x12))
that1c <- summary(model3b)$coef[2, 1]
that1d[j] <- summary(model3b)$coef[2, 1]
w1 <- 1/(df4$x^(2*that1c))
model7 <- lm(y ~ x -1, weights = w1, data = df4)
b6hatj[j] <- summary(model7)$coef[1, 1]
that0c <- summary(model7)$sigma
}

#Imputation
that1e <- c(that0b, that1d)
iteration1 <- c(1:(iter+1))
b1hat <- b6hatj[iter]

e1hat <- rnorm(n1, 0, that0c*x1^that1c)
y1hat <- b1hat*x1 + e1hat
y4 <- y2
y4[is.na(y2)==TRUE] <- y1hat[is.na(y2)==TRUE]
dfimp[, m + 2] <- y4
par(mfrow=c(2, 2))
plot(iteration1, that1e, type = "l", xlab = "iteration", ylab = "theta", main = "Convergence Check for Theta")
barplot(m/M*100, ylim = c(0, 100), main = "Progress Bar")
}
dfimp[, 1] <- y2
dfimp[, 2] <- x1
dfimp
}

###################
mranalyze2 <- function(data, alpha = 0.05){

data1 <- data
M <- ncol(data1) - 2
n1 <- nrow(data1)

ymean <- NULL
ymeanvar <- NULL
corxy <- NULL

for(m in 1:M){
y4 <- data1[ , m + 2]
x1 <- data1[ , 2]
ymean[m] <- mean(y4)
ymeanvar[m] <- var(y4)/n1
corxy[m] <- cor(y4, x1)
}
ybar <- mean(ymean)
wbar <- mean(ymeanvar)
bbar <- (1/(M-1)*sum((ymean - mean(ymean))^2))
se <- sqrt(wbar + (1+1/M)*bbar)

zm <- 1/2*log((1 + corxy)/(1 - corxy))
zmbar <- mean(zm)
rbar <- (exp(2*zmbar)-1)/(exp(2*zmbar)+1)
corr <- rbar

L1 <- (bbar + bbar/M)/(wbar + (1+1/M)*bbar)
k1 <- 1
v1 <- (M-1)/(L1^2)
v2 <- (n1 - k1 + 1)/(n1 - k1 +3)*(n1 - k1)*(1 - L1)
df <- (v1*v2)/(v1 + v2)

c <- qt(1-alpha/2, df)
CI.UL <- ybar + c*se
CI.LL <- ybar - c*se

output <- round(data.frame(ybar, se, CI.LL, CI.UL, corr, df), 4)
rownames(output) <- "output"
output
}

###################
mrdiag <- function(data){

data1 <- data
M <- ncol(data1) - 2
n1 <- nrow(data1)

y1 <- na.omit(data1[ , 1])

minx <- min(min(y1), min(data1[ , 3:(M+2)]))
maxx <- max(max(y1), max(data1[ , 3:(M+2)]))

maxy1 <- max(density(y1)$y)
maxy2 <- NULL
for(m in 1:M){
maxy2[m] <- max(density(data1[ , m + 2])$y)
}
maxy <- max(maxy1, maxy2)

plot(density(y1), ylim = c(0, maxy), xlim = c(minx, maxx), type = "n", xlab = "Values of Target Variable", main = "Diagnostic Plot")

for(m in 1:M){
y4 <- data1[ , m + 2]
par(new = TRUE)
plot(density(y4), ylim = c(0, maxy), xlim = c(minx, maxx), xlab = "", main = "", col = 2)
}
par(new = TRUE)
plot(density(y1), ylim = c(0, maxy), xlim = c(minx, maxx), xlab = "", main = "")

}
