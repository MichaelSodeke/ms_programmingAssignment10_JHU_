# [1_regScript.r]

source("my_theme.r")
windowsFonts()
loadfonts(device = "win")


# [case 1]: full model
# “Is there an association between college major category and income?”
require(collegeIncome)
data(college)
	#--[1]. model: one-way (unbalanced) model with q = 3 covariates
df <- college
df <- college %>% dplyr::select(median, major_category, perc_college_jobs, perc_non_college_jobs)
df <- na.omit(df)
df <- arrange(df, major_category)
df[,2] <- as.factor(df[,2])
lapply(df, head)
df %>% group_by(major_category, .add=TRUE) %>% group_nest()

Z <- as.matrix( ( model.matrix( median ~ major_category, data=df) ) ) # full-rank
X <- as.matrix( cbind(df$perc_college_jobs, df$perc_non_college_jobs) ) # full-rank
colnames(X) <- c("perc_college_jobs", "perc_non_college_jobs")
Y <- df$median
	#--[2]. estimating parameters
P <- Z %*% ginv(t(Z) %*% Z) %*% t(Z)
I <- diag(dim(P)[1])
B_ <- solve(t(X) %*% (I - P) %*% X) %*% (t(X) %*% (I - P) %*% Y); colnames(B_) <- "Estimates"
A_ <- ginv(t(Z) %*% Z) %*% (t(Z) %*% Y) - ginv(t(Z) %*% Z) %*% (t(Z) %*% (X %*% B_))
rownames(A_) <- colnames(Z)
rownames(B_) <- colnames(X)
U <- cbind(Z, X)
GAMMA <- rbind(A_, B_)
Y_ <- U %*% GAMMA
	#--[3]. degrees of freedom
n <- nrow(U); k <- qr(U)$rank
	#--[4]. sum of squares
ss.res <- (t(Y) %*% Y) - (t(GAMMA) %*% t(U) %*% Y)
ss.reg <- t(Y_ - mean(Y)) %*% (Y_ - mean(Y))
ss.tot <- ss.res + ss.reg
R.2 <- ss.res/ss.reg
	#--[5]. mean squares
ms.res <- ss.res/(n-k)
SE.. <- matrix( sqrt( diag( ms.res[1] * solve(t(U) %*% U) ) ) );  colnames(SE..) <- "Std. Error"
	#--[6]. hypothesis test | t-test (linear relation)
t.B_ <- GAMMA / SE..;  colnames(t.B_) <- "t value"
t.alpha <- qt(1-(0.05/2), df=n-k)
p_val.B_ <- 2*pt(-abs(t.B_), df=n-k);  colnames(p_val.B_) <- "Pr(>|t|)"
p_val.alpha <- 2*pt(-abs(t.alpha), df=n-k)
	#--[7].  influencial observations and leverage | https://data.princeton.edu/wws509/r/c2s9
H <- X %*% solve(t(X) %*% X) %*% t(X)
create.hii <- function(mat=NULL)
{
	vec <- matrix(rep(0, nrow(mat)))
	for (i in 1:nrow(mat))
	{
		for (j in 1:ncol(mat))
		{
			if ( i == j)
			{
				vec[i] <- mat[i,j] 
			}
		}
	}
	vec
}
hii <- create.hii(mat=H) # leverage
e_ <- Y - Y_; # residual
r_ <- e_/sqrt(ms.res[1] * (1 - hii)) # studentized residual
D_ <- ((r_^2)/(k+1)) * ((hii)/(1-hii)) # Cook's distance
inf.measures <- cbind(head(Y), head(Y_), head(e_), head(hii), head(r_),  head(D_))
colnames(inf.measures) <- c("yi", "yi_", "ei_", "hii", "ri_", "Di_")
infl <- c( "hii_LEVERAGE"=(2*(k + 1))/n )
	#--[8]. VARIFY RESULTS
r.fit1 <- summary( lm(median ~ major_category + perc_college_jobs + perc_non_college_jobs, data=df) ) # test for [linear relation] with t-test
summary( aov(median ~ major_category * perc_college_jobs, data = df) ) # test for [effect] with f-test
my.fit1 <- list(Y.values.vs.fitted=head(cbind(Y, Y_)), predictors=head(cbind(Z, X)), coefficients=cbind(GAMMA, SE.., t.B_, p_val.B_))
options(scipen = 999)
my.fit2 <- list(infl, inf.measures)
options(scipen = 0)
	#--[9]. MODEL SELECTION
fit1 <- lm(median ~ major_category, data=df)
fit2 <- lm(median ~ major_category + perc_college_jobs, data=df)
fit3 <- lm(median ~ major_category + perc_college_jobs + perc_non_college_jobs, data=df)
mv1 <- anova(fit1, fit2, fit3)

fit4 <- aov(median ~ major_category + perc_college_jobs, data = df)
fit5 <- aov(median ~ major_category * perc_college_jobs, data = df)
anova(fit4, fit5)

VIF1 <- vif( lm(median ~ major_category + perc_college_jobs + perc_non_college_jobs, data=df) ) # variance inflation factor
VIF2 <- vif( lm(median ~ major_category + perc_college_jobs, data=df) )
VIF <- list(VIF1, VIF2) 
	#--[10]. PLOT
options(scipen = 999)
grid <- data.frame(x=X[,1], y=X[,2], z=Y)
png(file="./report_and_figures/figures/plot1.png", width=45, height=25, units="cm", res=500) # regression [plot]
p.reg <- cloud(z ~ x + y,
	  data=grid,
	  pch=20,
	  alpha=0.4,
	  col="black",
	  scales=list(arrows=FALSE, cex=0.5, tck=1, lty=0, family="Segoe UI"),
	  aspect=c(1,1),
	  screen=list(z=-140, x=-65),
	  main=list(label="ANCOVA (one-way) with q = 2 covariates", family="Segoe UI", adj=-1),
	  xlab=list(label="Degree Jobs", cex=0.60, family="Segoe UI", rot=10),
	  ylab=list(label="Non-Degree Jobs", cex=0.60, family="Segoe UI", rot=-15),
	  zlab=list(label="Median income", cex=0.60, family="Segoe UI", rot=90) )
print(p.reg)
dev.off()

grid <- data.frame(x=X[,1], y=X[,2], z=e_); colnames(grid) <- c("x", "y", "z")
png(file="./report_and_figures/figures/plot2.png", width=45, height=25, units="cm", res=500) # residuals vs. fitted [plot]
p.reg <- cloud(z ~ x + y,
	  data=grid,
	  pch=20,
	  alpha=0.4,
	  col="black",
	  scales=list(arrows=FALSE, cex=0.5, tck=1, lty=0, family="Segoe UI"),
	  aspect=c(1,1),
	  screen=list(z=-140, x=-65),
	  main=list(label="Residuals vs. Fitted", family="Segoe UI", adj=-1),
	  xlab=list(label="Degree Jobs", cex=0.60, family="Segoe UI", rot=10),
	  ylab=list(label="Non-Degree Jobs", cex=0.60, family="Segoe UI", rot=-15),
	  zlab=list(label="Residuals", cex=0.60, family="Segoe UI", rot=90) )
print(p.reg)
dev.off()




# [case 2]: adjusted model
# “Is there an association between college major category and income?”
# one-way ANCOVA
require(collegeIncome)
data(college)
	#--[1]. model: one-way (unbalanced) model with q = 3 covariates
df <- college
df <- college %>% dplyr::select(median, major_category, perc_college_jobs)
df <- na.omit(df)
df <- arrange(df, major_category)
df[,2] <- as.factor(df[,2])
lapply(df, head)
df %>% group_by(major_category, .add=TRUE) %>% group_nest()
Z <- as.matrix( ( model.matrix( median ~ major_category, data=df) ) ) # full-rank
X <- as.matrix( cbind(df$perc_college_jobs) ) # full-rank
colnames(X) <- c("perc_college_jobs")
Y <- df$median
	#--[2]. estimating parameters
P <- Z %*% ginv(t(Z) %*% Z) %*% t(Z)
I <- diag(dim(P)[1])
B_ <- solve(t(X) %*% (I - P) %*% X) %*% (t(X) %*% (I - P) %*% Y); colnames(B_) <- "Estimates"
A_ <- ginv(t(Z) %*% Z) %*% (t(Z) %*% Y) - ginv(t(Z) %*% Z) %*% (t(Z) %*% (X %*% B_))
rownames(A_) <- colnames(Z)
rownames(B_) <- colnames(X)
U <- cbind(Z, X)
GAMMA <- rbind(A_, B_)
Y_ <- U %*% GAMMA
	#--[3]. degrees of freedom
n <- nrow(U); k <- qr(U)$rank
	#--[4]. sum of squares
ss.res <- (t(Y) %*% Y) - (t(GAMMA) %*% t(U) %*% Y)
ss.reg <- t(Y_ - mean(Y)) %*% (Y_ - mean(Y))
ss.tot <- ss.res + ss.reg
R.2 <- ss.res/ss.reg
	#--[5]. mean squares
ms.res <- ss.res/(n-k)
SE.. <- matrix( sqrt( diag( ms.res[1] * solve(t(U) %*% U) ) ) );  colnames(SE..) <- "Std. Error"
	#--[6]. hypothesis test | t-test (linear relation)
t.B_ <- GAMMA / SE..;  colnames(t.B_) <- "t value"
t.alpha <- qt(1-(0.05/2), df=n-k)
p_val.B_ <- 2*pt(-abs(t.B_), df=n-k);  colnames(p_val.B_) <- "Pr(>|t|)"
p_val.alpha <- 2*pt(-abs(t.alpha), df=n-k)
	#--[7].  influencial observations and leverage | https://data.princeton.edu/wws509/r/c2s9
H <- X %*% solve(t(X) %*% X) %*% t(X)
create.hii <- function(mat=NULL)
{
	vec <- matrix(rep(0, nrow(mat)))
	for (i in 1:nrow(mat))
	{
		for (j in 1:ncol(mat))
		{
			if ( i == j)
			{
				vec[i] <- mat[i,j] 
			}
		}
	}
	vec
}
hii <- create.hii(mat=H) # leverage
e_ <- Y - Y_; # residual
r_ <- e_/sqrt(ms.res[1] * (1 - hii)) # studentized residual
D_ <- ((r_^2)/(k+1)) * ((hii)/(1-hii)) # Cook's distance
inf.measures <- cbind(head(Y), head(Y_), head(e_), head(hii), head(r_),  head(D_)); colnames(inf.measures) <- c("yi", "yi_", "ei_", "hii", "ri_", "Di_")
infl <- c( "hii_LEVERAGE"=(2*(k + 1))/n )
	#--[8]. VARIFY RESULTS
r.fit2 <- summary( lm(median ~ major_category + perc_college_jobs, data=df) )
my.fit3 <- list(Y.values.vs.fitted=head(cbind(Y, Y_)), predictors=head(cbind(Z, X)), coefficients=cbind(GAMMA, SE.., t.B_, p_val.B_))
options(scipen = 999)
my.fit4 <- list(infl, inf.measures)
options(scipen = 0)
	#--[9]. PLOT
p.hist <- ggplot(data=df, aes(x=major_category)) # histogram [plot]
p.hist <- p.hist + geom_bar()

options(scipen = 999)
p.reg <- ggplot(data=data.frame(Y_, Y, X), aes(x=X, y=Y)) # regression [plot]
p.reg <- p.reg + geom_segment(aes(xend=X, yend=Y_), color="black", alpha=0.2)
p.reg <- p.reg + geom_point(alpha=0.4)
p.reg <- p.reg + geom_abline(aes(intercept=GAMMA[1,], slope=GAMMA[16,]), color='grey30') # reguard regressor X
p.reg <- p.reg + geom_abline(aes(intercept=GAMMA[1,] + GAMMA[2,] + GAMMA[3,] + GAMMA[4,] +
					   GAMMA[5,] + GAMMA[6,] + GAMMA[7,] + GAMMA[8,] +
					   GAMMA[9,] + GAMMA[10,] + GAMMA[11,] + GAMMA[12,] +
					   GAMMA[13,] + GAMMA[14,] + GAMMA[15,],
					   slope=GAMMA[16,]), color='grey30') # reguard regressor X
p.reg <- p.reg + labs(title = "ANCOVA (one-way) with q = 2 covariates",
		x = "Percentage of Jobs requiring a Degree",
		y = "Median income")
p.reg <- p.reg + my_theme()
ggsave("./report_and_figures/figures/plot4.png", p.reg, dpi=500, width = 45, height = 25, units="cm", device="png")

p.res <- ggplot(data=data.frame(X, e_), aes(x=X, y=e_)) # residuals vs. fitted [plot]
p.res <- p.res + geom_segment(aes(xend=X, yend=0), color="grey30", alpha=0.2)
p.res <- p.res + geom_point(alpha=0.4, color="black")
p.res <- p.res + geom_hline(yintercept=0, color='black')
p.res <- p.res + labs(title = "Residuals vs. Fitted",
		x = "Fitted",
		y = "Residuals")
p.res <- p.res + my_theme()

p.sd <- ggplot(data=data.frame(X, r_), aes(x=X, y=r_)) # scale-location [plot]
p.sd <- p.sd + geom_segment(aes(xend=X, yend=0), color="grey30", alpha=0.2)
p.sd <- p.sd + geom_point(alpha=0.4, color="black")
p.sd <- p.sd + geom_hline(yintercept=0, color='black')
p.sd <- p.sd + labs(title = "Scale-Location",
		x = "Fitted",
		y = "Standardized Residuals")
p.sd <- p.sd + my_theme()

p.lv <- ggplot(data=data.frame(r_, hii), aes(x=hii, y=r_)) # standardized residuals vs. leverage [plot]
p.lv <- p.lv + geom_segment(aes(xend=hii, yend=0), color="grey30", alpha=0.2)
p.lv <- p.lv + geom_point(alpha=0.4, color="black")
p.lv <- p.lv + geom_hline(yintercept=0, color='black')
p.lv <- p.lv + labs(title = "Standardized Residuals vs. Leverage",
		x = "Leverage",
		y = "Standardized Residuals")
p.lv <- p.lv + my_theme()

plot <- ggarrange(p.res, p.sd, p.lv, nrow=1, ncol=3)
ggsave("./report_and_figures/figures/plot3.png", plot, dpi=500, width = 45, height = 25, units="cm", device="png")