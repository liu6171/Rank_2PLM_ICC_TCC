#Supplemental Online Material
#R Code to Create ICC and TCC plots Based on the MUPP-2PL and Triplet-2PL Models


#Set working directory, for example, “C:\ICC TCC”
setwd("C:\\ICC TCC")

#The following code has been tested to work with the mirt package version 1.38.1

#load libraries

library(mirt)
library(mvtnorm)
library(psych)
library(lattice)
library(Deriv)
library(latticeExtra)

#In the code below, the functions I wrote are explicitly specified; all the other functions are from the mirt package or R. 
#You can always use the "dump" function to see the source code of a function. For example,

dump("mirt", file="", control= NULL)

#simulate MUPP-2PL data using simdata() function in mirt.

#There are 5 traits
Ntrait.2 <- 5

#Each trait measured by 8 statements
Nstate.per.trait.2 <- 8

#the Nitem X 2 matrix is 2 X possible combinations of 2 out of 5 traits, containing Trait ID measured by each statement in each item
Dmatrix.2 <- rbind(t(combn(Ntrait.2,2)), t(combn(Ntrait.2,2))) 
colnames(Dmatrix.2) <- c("dim1", "dim2")
#create the Q matrix
Q.2 <- matrix(0, nr=nrow(Dmatrix.2), ncol=Ntrait.2, dimnames = list(NULL, paste0('F', 1:Ntrait.2)))
for (i in 1:nrow(Q.2)){
	Q.2[i,Dmatrix.2[i,]] <-1
}

#all intertrait correlations are set to 0.5
Trait.cor.2 <- matrix(0.5, nc=Ntrait.2, nr=Ntrait.2)
diag(Trait.cor.2) <-1

#Both the discrimination and intercept parameters are drawn from a standard normal 
#based on Equation 4 in the paper
a.2 <- matrix(round(rnorm(nrow(Dmatrix.2)*Ntrait.2),2), nc=Ntrait.2)*Q.2
d.2 <- matrix(round(rnorm((2*nrow(Dmatrix.2))),2), nc=2)

#Or draw the discrimination parameters for a lognormal with mean =1 and std=0.25 at the log scale
#change the following discrimination parameters to negative: both statements in Items 1 and 2 and the first statement in Item 3 
#a.2 <- matrix(round(lognorm(nrow(Dmatrix.2)*Ntrait.2, 1, 0.5^2),2), nc=Ntrait.2)*Q.2
#a.2[1:2,] <- -a.2[1:2,]
#a.2[3,1] <- -a.2[3,1]

#Transform to M2PL parameters in mirt
a.M2PL <- a.2
for (i in 1:nrow(Dmatrix.2)) a.M2PL[i,Dmatrix.2[i,2]] <- -a.M2PL[i,Dmatrix.2[i,2]]
d.M2PL <- d.2[,1]-d.2[,2]


#simulate 1000 test takers and 20 items (40 statements); 1 = selecting the first statement,
#0 = selecting the second statement
sim.data.2 <- simdata(a.M2PL, d.M2PL, N=1000, itemtype='dich', sigma=Trait.cor.2)

#create model object for mirt function
COV <- matrix(T, nrow=Ntrait.2, nc=Ntrait.2)
diag(COV) <- F
Model.MUPP.2PL.5 <- mirt.model(Q.2, COV=COV)

#For a model with a direct estimation of item parameters and intertrait correlations 
#set initial values to -1 for negative discrimination, 1 for positive discrimination, and 0 for intercept
tmp <- mirt(sim.data.2, Model.MUPP.2PL.5, '2PL', pars = 'values')
tmp[substr(tmp$name,1,1)=="a","value"] <- as.vector(t(sign(a.M2PL)))
tmp[substr(tmp$name,1,1)=="d" & tmp$est,"value"] <- 0
sv.MUPP.2PL.5.est <- tmp

#For a model with fixed item parameters and intertrait correlations
#set initial values to the true values or previous estimates from a Likert scale
tmp <- mirt(sim.data.2, Model.MUPP.2PL.5, '2PL', pars = 'values')
tmp[substr(tmp$name,1,1)=="a","value"] <- as.vector(t(a.M2PL))
tmp[substr(tmp$name,1,1)=="d" & tmp$est,"value"] <- d.M2PL

#Fix intertrait correlations
tmp[grepl("COV", tmp$name, fixed=T), ]$value <- Trait.cor.2[lower.tri(Trait.cor.2,T)]
tmp[grepl("COV", tmp$name, fixed=T), ]$est <- FALSE

sv.MUPP.2PL.5.fix <- tmp

#run mirt to estimate item parameters using "QMCEM" method; 'MHRM' method produces better estimates 
#but takes a much longer time
MUPP.2PL.5.est <- mirt(sim.data.2, Model.MUPP.2PL.5, '2PL', pars = sv.MUPP.2PL.5.est, method='QMCEM', 
SE=T, large=T) 

#run mirt by fixing item parameters and using "QMCEM" method; the 'MHRM' method cannot be 
#used in the fixed run. 

MUPP.2PL.5.fix <- mirt(sim.data.2, Model.MUPP.2PL.5, '2PL', pars = sv.MUPP.2PL.5.fix, TOL=NaN, method='QMCEM', 
SE=F, large=T) 

##################################################################
#Function to extract item parameter estimates of MUPP-2PL 
#
#Description:
#Given a mirt object of M2PL estimation, extract item parameter estimates and transform to those
#based on Equation 4 in the paper. 
#
#Usage:
#extract.par.MUPP.2PL(mirt.obj, dim, d1)
#
#Arguments:
#mirt.obj -- an object returned by the mirt function estimating M2PL
#dim -- a Nitem-row matrix including trait ID measured by each statement in each item, with column 
#		names as "dim1", "dim2".
#d1 -- a Nitem vector including the intercept of the first statement in each item 
#
#Return:
#a Nitem X 6 matrix including item parameter estimates and trait ID for each statement in each item
#with column names a1, b1, a2, b2, dim1, dim2. 
#
#Note:
#In the output file, intercept is labeled by "b".   
#In Equation 4 in the paper, the intercept of the first statement in each item
#is fixed during estimation
##########################################################################

extract.par.MUPP.2PL <-
function (mirt.obj, dim, d1) 
{
    d1 <- round(d1, 2)
    tmp <- coef(mirt.obj, simplify = T)
    tmp1 <- round(extract.par.est2(tmp[[1]], dim[, c("dim1", 
        "dim2")]), 2)
    tmp1[, "a2"] <- -tmp1[, "a2"]
    tmp1 <- cbind(tmp1, d1)
    tmp1[, 3] <- d1 - tmp1[, 3]
    tmp1 <- cbind(tmp1, dim[, c("dim1", "dim2")])
    tmp1 <- tmp1[, c(1, 4, 5, 2, 3, 6)]
    colnames(tmp1) <- c("a1", "b1", "dim1", "a2", "b2", "dim2")
    tmp1
}

extract.par.est2 <-
function (par, dim) 
{
    tmp1 <- NULL
    for (i in 1:nrow(par)) {
        tmp1 <- rbind(tmp1, par[i, dim[i, ]])
    }
    tmp1 <- cbind(tmp1, par[, "d"])
    tmp1
}


MUPP.2PL.5.par.ets <- extract.par.MUPP.2PL(MUPP.2PL.5.est, Dmatrix.2, d.2[,1])
MUPP.2PL.5.par.fix <- extract.par.MUPP.2PL(MUPP.2PL.5.fix, Dmatrix.2, d.2[,1])

#estimate MAP trait scores
MUPP.2PL.5.MAP.ets <- fscores(MUPP.2PL.5.est, method="MAP", QMC=T, full.scores.SE=T) 
MUPP.2PL.5.MAP.fix <- fscores(MUPP.2PL.5.fix, method="MAP", QMC=T, full.scores.SE=T) 

#item fit: Orlando and Thissen’s (2000) chi-squared statistic, S_X2
MUPP.2PL.5.est.X2 <- itemfit(MUPP.2PL.5.est, fit_stats ="S_X2", Theta=MUPP.2PL.5.MAP.est[,1:Ntrait.2], QMC=TRUE)
MUPP.2PL.5.fix.X2 <- itemfit(MUPP.2PL.5.fix, fit_stats ="S_X2", Theta=MUPP.2PL.5.MAP.fix[,1:Ntrait.2], QMC=TRUE)


#########################################################################
#Function to create MUPP-2PL ICC data that are used for drawing ICC plots 
#
#Description:
#Create observed and expected ICC data for one item.
#
#Usage:
#ICC.MUPP.2PL(item, model.est, dimcor, par.est, theta.up, theta.low, 
#    theta.interval = NULL, nquad = NULL, data, theta.score)
#
#Arguments:
#item -- numeric, item sequence number (e.g., 1) of the item whose ICC is calculated
#model.est -- an object returned by the mirt function estimating M2PL
#dimcor -- the covariance/correlation matrix of trait scores
#par.est -- a Nitem-row matrix containing the following colnames: 
#           discrimination and intercept parameters for Statements 1 and 2 in each item: "a1"   "b1"   "a2"   "b2"  
#		trait number(1:Ntrait) measured by each statement in each item: "dim1" "dim2"
#theta.up -- the upper bound of the theta range
#theta.low -- the lower bound of the theta range 
#theta.interval -- equal spaced interval in the theta range
#nquad -- number of theta points in the theta range; if both theta.interval and nquad are provided, theta.interval overwrites nquad
#data -- a matrix (Nexaminee X Nitem)containing item score (1=select the first statement,0=select the second statement)
#theta.score -- a matrix (Nexaminee X Ntrait) containing trait score(s) for each examinee
#
#Return:
#A list with 6 elements. 
#
#1. matrix  :  item true score at a theta vector   
#  ..$ F1       : trait point measured by the first statement
#  ..$ F2       : trait point measured by the second statement
#  ..$ weight   : weight (probability) of the theta vector (F1, F2)
#  ..$ Item&item : expected item scores at (F1, F2)
#  ..$ Item&item.w : Item&item*weight
#
#2. matrix  : item true score for each examinee 
#  ..$ Item&item  : observed item score
#  ..$ Item&item.t  : expected item score 
#  ..$ F1      : trait score estimate by the first statement
#  ..$ F2      : trait score estimate by the second statement
#
#3. data.frame  : expected (marginal) item true score at a theta point 
#  ..$ Theta  : a theta point
#  ..$ F1.Item&item  : expected item true score at a point of the trait measured by the first statement
#  ..$ F2.Item&item  : expected item true score at a point of the trait measured by the second statement
#
#4. data.frame  :  item true score at a theta vector with negative statement(s) changed to positive one(s)    
#  ..$ F1       : trait point measured by the first statement
#  ..$ F2       : trait point measured by the second statement
#  ..$ weight   : weight (probability) of the theta vector (F1, F2)
#  ..$ Item&item : expected item scores at (F1, F2)
#  ..$ Item&item.w : Item&item*weight
#
#5. data.frame  : expected (marginal) item true score at a theta point with negative statement(s) changed to positive one(s)
#  ..$ Theta  : a theta point
#  ..$ F1.Item&item  : expected item true score at a point of the trait measured by the first statement
#  ..$ F2.Item&item  : expected item true score at a point of the trait measured by the second statement
#
#6. matrix  : adjusted observed item score for each examinee with negative statement(s) changed to positive one(s) 
#  ..$ Item&item  : adjusted observed item score
#  ..$ F1      : trait score estimate by the first statement
#  ..$ F2      : trait score estimate by the second statement
#
#
##########################################################################


ICC.MUPP.2PL <-
function (item, model.est, dimcor, par.est, theta.up, theta.low, 
    theta.interval = NULL, nquad = NULL, data, theta.score) 
{
    block.size <- 2
    nf <- nrow(dimcor)
    if (is.null(theta.interval)) {
        tmp <- round(seq(theta.low, theta.up, length.out = nquad), 
            2)
    }
    else {
        tmp <- round(seq(theta.low, theta.up, by = theta.interval), 
            2)
    }
    if (block.size == 3) {
        tmp1 <- as.matrix(expand.grid(tmp, tmp, tmp))
    }
    else if (block.size == 2) {
        tmp1 <- as.matrix(expand.grid(tmp, tmp))
    }
    tmp1.1 <- as.matrix(expand.grid(tmp, tmp))
    tmp1.dim <- par.est[item, paste0("dim", 1:block.size)]
    colnames(tmp1) <- paste0("F", 1:length(tmp1.dim))
    tmp1.d <- dmvnorm(tmp1, sigma = dimcor[tmp1.dim, tmp1.dim])
    tmp1.d <- tmp1.d/sum(tmp1.d)
    tmp2.d <- dnorm(tmp)
    tmp2.d <- cbind(tmp, tmp2.d/sum(tmp2.d))
    colnames(tmp2.d) <- c("Group", "Weight.s")
    theta <- matrix(0, nc = nf, nr = nrow(tmp1), byrow = T)
    theta[, par.est[item, paste0("dim", 1:block.size)]] <- tmp1
    i.exp <- expected.item(extract.item(model.est, item), Theta = theta.score[, 
        1:nf], min = 0)
    icc.exp <- expected.item(extract.item(model.est, item), Theta = theta, 
        min = 0)
    icc.exp.t <- cbind(tmp1, tmp1.d, icc.exp, icc.exp * tmp1.d)
    colnames(icc.exp.t) <- c(paste0("F", 1:length(tmp1.dim)), 
        "weight", paste0("Item", item), paste0("Item", item, 
            ".w"))
    i.exp.t <- cbind(data, i.exp, theta.score[, tmp1.dim])
    colnames(i.exp.t) <- c(paste0("Item", item), paste0("Item", 
        item, ".t"), paste0("F", 1:length(tmp1.dim)))
    for (i in 1:block.size) {
        tmp3 <- describeBy(icc.exp.t[, paste0("Item", item, ".w")], 
            group = icc.exp.t[, paste0("F", i)], mat = T, digits = 15)[, 
            c("group1", "vars", "n", "mean")]
        tmp4 <- data.frame(tmp3[, 1:2], mean = tmp3[, 4] * tmp3[, 
            3])
        tmp5 <- cbind(tmp, matrix(tmp4[, 3], nr = length(tmp)))
        colnames(tmp5) <- c("Group", paste0("F", i, paste0(".Item", 
            item)))
        if (i == 1) 
            icc.E <- tmp5
        else icc.E <- merge(icc.E, tmp5, by = "Group", sort = F)
    }
    icc.E1 <- merge(tmp2.d, icc.E, sort = F)
    icc.E2 <- cbind(Theta = icc.E1[, 1], icc.E1[, 3:ncol(icc.E1)]/icc.E1[, 
        2])
    sign.tmp <- sign(par.est[item, paste0("a", 1:block.size)])
    sign.tmp[sign.tmp == 0] <- 1
    sign.n <- which(sign.tmp == -1)
    tmp8 <- icc.exp.t
    icc.E2.pa <- icc.E2
    if (sum(sign.tmp) == -2) {
        tmp8[, paste0("Item", item)] <- 1 - tmp8[, paste0("Item", 
            item)]
        icc.E2.pa[, 2:3] <- 1 - icc.E2[, 2:3]
    }
    else if (sum(sign.tmp) == 0) {
        a <- abs(par.est[item, paste0("a", 1:block.size)])
        b <- par.est[item, paste0("b", 1:block.size)]
        b[sign.n] <- -b[sign.n]
        tmp8[, paste0("Item", item)] <- P.MUPP.2PL.s(c(a, b), 
            tmp8[, "F1"], tmp8[, "F2"])[, 2]
        tmp8[, paste0("Item", item, ".w")] <- tmp8[, paste0("Item", 
            item)] * tmp8[, "weight"]
        for (i in 1:block.size) {
            tmp3 <- aggregate(tmp8[, paste0("Item", item, ".w")], 
                by = list(tmp8[, paste0("F", i)]), FUN = sum)
            colnames(tmp3) <- c("Group", paste0("F", i, ".Item", 
                item))
            if (i == 1) {
                tmp5 <- tmp3
            }
            else {
                tmp5 <- merge(tmp5, tmp3, by = "Group", sort = F)
            }
        }
        icc.E1.pa <- merge(tmp2.d, tmp5, sort = F)
        icc.E2.pa <- cbind(Theta = icc.E1.pa[, 1], icc.E1.pa[, 
            3:ncol(icc.E1.pa)]/icc.E1.pa[, 2])
    }
    tmp <- i.exp.t[, c(paste0("Item", item), paste0("F", 1:length(tmp1.dim)))]
    if (sum(sign.tmp) == -2) {
        tmp[, 1] <- 1 - tmp[, 1]
    }
    else if (sum(sign.tmp) == 0) {
        a <- abs(par.est[item, paste0("a", 1:block.size)])
        b <- par.est[item, paste0("b", 1:block.size)]
        b[sign.n] <- -b[sign.n]
        tmp[, 1] <- P.MUPP.2PL.s(c(a, b), tmp[, "F1"], tmp[, 
            "F2"])[, 2]
    }
    i.exp.pa <- tmp
    return(list(icc.exp.t, i.exp.t, icc.E2, tmp8, icc.E2.pa, 
        i.exp.pa))
}

#Item response function (Equation 4 in the paper) for the MUPP-2PL
P.MUPP.2PL.s <- function(par,theta1, theta2){
	pp <- 1/(1+exp(-(par[1]*theta1-par[2]*theta2+par[3]-par[4])))	
	cbind(1-pp, pp)
}
#Create ICC data for Item 1

MUPP.2PL.5.fix.ICC.1 <- ICC.MUPP.2PL(item=1, 
	model.est=MUPP.2PL.5.fix, 
	dimcor=Trait.cor.2, 
	par.est=MUPP.2PL.5.par.fix, 
	theta.up=3, 
	theta.low=-3, 
	theta.interval=NULL, 
	nquad=31,  
	data= sim.data.2[, 1], 
	theta.score=MUPP.2PL.5.MAP.fix)

#Create ICC data for all items in a test

MUPP.2PL.5.fix.ICC <- vector('list', 20)
for ( i in 1:20){
	MUPP.2PL.5.fix.ICC[[i]] <- ICC.MUPP.2PL(item=i, 
	model.est=MUPP.2PL.5.fix, 
	dimcor=Trait.cor.2, 
	par.est=MUPP.2PL.5.par.fix, 
	theta.up=3, 
	theta.low=-3, 
	theta.interval=NULL, 
	nquad=31,  
	data= sim.data.2[, i], 
	theta.score=MUPP.2PL.5.MAP.fix)
}

#########################################################################
#Function to draw MUPP-2PL ICC plots using the ICC data output from the ICC.MUPP.2PL function 
#
#Description:
#Draw one or multiple MUPP-2PL ICC plots using the ICC data output from the ICC.MUPP.2PL function
#
#Usage:
#ICC.plot.MUPP.2PL(ICC.data, plot.item, par.est, PageNo = T, fit) 
#
#Arguments:
#ICC.data -- a list with Nitem elements being the ICC data output of the ICC.MUPP.2PL function for all items in a test 
#plot.item -- a numeric vector containing item number(s) for which the ICC(s) is drawn 
#par.est -- a Nitem-row matrix containing the following colnames: 
#           discrimination and intercept parameters for Statements 1 and 2 in each item: "a1"   "b1"   "a2"   "b2"  
#		trait number(1:Ntrait) measured by each statement in each item: "dim1" "dim2"
#PageNo -- logic; T = print page number at the bottom center, F = no page number
#fit -- the output of S_X2 item fit statistics from the itemfit function in the mirt package
#
#Return:
#Output the ICC plot(s) for the input item(s)
##########################################################################

ICC.plot.MUPP.2PL <-
function (ICC.data, plot.item, par.est, PageNo = T, fit) 
{
    block.size <- 2
    for (i in plot.item) {
        i.dim <- par.est[i, paste0("dim", 1:block.size)]
        par.item <- matrix(par.est[i, paste0(c("a", "b"), rep(1:block.size, 
            each = 2))], ncol = 2, byrow = T)
        state.text <- c("A", "B", "C")[1:block.size]
        group.text <- paste0("F", i.dim, "-Statement ", state.text, 
            " (a= ", par.item[, 1], ", b= ", par.item[, 2], ")")
        group.text <- factor(group.text, levels = group.text)
        i.dat <- ICC.data[[i]][[2]][, paste0("Item", i)]
        i.dat <- i.dat[!is.na(i.dat)]
        if (block.size == 3) 
            key.text <- paste0("Item", i, " N=", length(i.dat), 
                " P=", round(mean(i.dat), 2), " Zh=", round(fit[i, 
                  2], 2), " P(Zh)=", round(pnorm(fit[i, 2]), 
                  2))
        else key.text <- paste0("Item", i, " N=", length(i.dat), 
            " P=", round(mean(i.dat), 2), " SX2=", round(fit[i, 
                2], 2), " P(SX2)=", round(fit[i, 5], 2))
        i.exp.tmp <- ICC.data[[i]][[2]][!is.na(i.dat), c(paste0("Item", 
            i), paste0("F", 1:block.size))]
        i.exp.tmp1 <- reshape(data.frame(i.exp.tmp), direction = "long", 
            varying = paste0("F", 1:block.size), v.names = "Theta", 
            times = paste0("F", 1:block.size), timevar = "group")
        i.exp.tmp1[i.exp.tmp1$group == "F2", 1] <- 1 - i.exp.tmp1[i.exp.tmp1$group == 
            "F2", 1]
        i.exp.tmp1$group <- rep(group.text, each = length(i.dat))
        f.tmp <- as.formula(paste(paste0("Item", i), "~Theta|group"))
        col.tmp <- sum(par.item[, 1] >= 0) + 1
        bg.col <- c(2, 3, 1)[col.tmp]
        plot1 <- xyplot(f.tmp, data = i.exp.tmp1, panel = function(x, 
            y, ...) {
            panel.grid(h = -1, v = -1, col.line = grey(0.9))
            panel.smoother(x, y, ..., family = "g", n = 100, 
                span = 0.9, lty = 2)
        }, ylab = list(label = "Expected Item Trait Score", col = 1, 
            cex = 1), strip = strip.custom(strip.names = c(F, 
            TRUE), par.strip.text = list(cex = 0.7, col = 1), 
            bg = trellis.par.get("strip.background")$col[[bg.col]]), 
            key = list(text = list(c("Expected", "Observed(with 95% CI)"), 
                cex = 0.7, col = c(2, 4)), lines = list(lty = 1:2, 
                col = c(2, 4), cex = 0.7), space = "top", border = F, 
                between = 1, adj = 0, column = 2), layout = c(1, 
                2), main = list(label = key.text, cex = 1, fontface = 1), 
            as.table = T, sub = if (PageNo) 
                list(label = i, cex = 0.7)
            else "")
        icc.E2.tmp <- reshape(ICC.data[[i]][[3]][, c("Theta", 
            paste0("F", 1:block.size, ".Item", i))], direction = "long", 
            varying = paste0("F", 1:block.size, ".Item", i), 
            v.names = paste0("Item", i), times = paste0("F", 
                1:block.size), timevar = "group")
        icc.E2.tmp[icc.E2.tmp$group == "F2", 3] <- 1 - icc.E2.tmp[icc.E2.tmp$group == 
            "F2", 3]
        plot2 <- xyplot(f.tmp, icc.E2.tmp, type = "l", col = 2, 
            layout = c(1, 2), as.table = T)
        print(plot1 + as.layer(plot2))
    }
}

#draw the first item's ICC plot to screen
ICC.plot.MUPP.2PL(ICC.data=MUPP.2PL.5.fix.ICC, plot.item=1, par.est=MUPP.2PL.5.par.fix,PageNo=F, fit=MUPP.2PL.5.fix.X2)

#draw the ICC plots of all items in test and save to a PDF file

pdf("MUPP_2PL_5_fix_ICC_plot.pdf")

ICC.plot.MUPP.2PL(ICC.data=MUPP.2PL.5.fix.ICC,  plot.item=1:20, par.est=MUPP.2PL.5.par.fix,PageNo=F, fit=MUPP.2PL.5.fix.X2)

dev.off()

#########################################################################
#Function to create TCC data for the MUPP-2PL
#
#Description:
#Create TCC data for the MUPP-2PL using the ICC data output from the ICC.MUPP.2PL function
#
#Usage:
#TCC.MUPP.2PL(ICC.data, par.est, theta.score)
#
#Arguments:
#ICC.data -- a list with Nitem elements being the ICC data output of the ICC.MUPP.2PL function for all items in a test 
#par.est -- a Nitem-row matrix containing the following colnames: 
#           discrimination and intercept parameters for Statements 1 and 2 in each item: "a1"   "b1"   "a2"   "b2"  
#		trait number(1:Ntrait) measured by each statement in each item: "dim1" "dim2"
#theta.score -- a matrix (Nexaminee X Ntrait) containing trait score(s) for each examinee
#
#Return: 
#A list with 4 elements. 
#
#1. matrix  :  expected trait test true score at a theta value with negative statement(s) changed to positive one(s)   
#  ..$ Theta    : trait score point defined by theta.up, theta.low, theta.interval, or nquad in the ICC.MUPP.2PL function  
#  ..$ F1.subscore : Trait 1's test scores
#  ..$ F2.subscore : Trait 2's test scores
#  ... 
#  ..$ F&Ntrait.subscore : last Trait's test scores
#  
#2. matrix  :  observed trait test scores for each examinee without changing the direction of negative statement(s)   
#  ..$ F1    : Trait 1's score estimate 
#  ..$ F2    : Trait 2's score estimate 
#  ... 
#  ..$ F&Ntrait : last Trait's score estimate
#  ..$ F1.subscore : Trait 1's observed test scores
#  ..$ F2.subscore : Trait 2's observed test scores
#  ... 
#  ..$ F&Ntrait.subscore : last Trait's observed test scores
#3. matrix  :  expected trait test true score at a theta value without changing the direction of negative statement(s)   
#  ..$ Theta    : trait score point defined by theta.up, theta.low, theta.interval, or nquad in the ICC.MUPP.2PL function  
#  ..$ F1.subscore : Trait 1's test scores
#  ..$ F2.subscore : Trait 2's test scores
#  ... 
#  ..$ F&Ntrait.subscore : last Trait's test scores
#  
#4. matrix  :  observed trait test scores for each examinee with negative statement(s) changed to positive one(s)   
#  ..$ F1    : Trait 1's score estimate 
#  ..$ F2    : Trait 2's score estimate 
#  ... 
#  ..$ F&Ntrait : last Trait's score estimate
#  ..$ F1.subscore : Trait 1's observed test scores
#  ..$ F2.subscore : Trait 2's observed test scores
#  ... 
#  ..$ F&Ntrait.subscore : last Trait's observed test scores
##########################################################################

TCC.MUPP.2PL <-
function (ICC.data, par.est, theta.score) 
{
    nf <- max(par.est[, c("dim1", "dim2")])
    block.size <- 2
    subscore.t.adj <- matrix(0, nr = nrow(ICC.data[[1]][[5]]), 
        nc = nf)
    subscore.t <- matrix(0, nr = nrow(ICC.data[[1]][[3]]), nc = nf)
    subscore.obs <- subscore.obs.adj <- matrix(0, nr = nrow(ICC.data[[1]][[2]]), 
        nc = nf)
    for (i in 1:length(ICC.data)) {
        tmp1.dim <- par.est[i, paste0("dim", 1:block.size)]
        tmp2 <- ICC.data[[i]][[5]]
        tmp3 <- ICC.data[[i]][[2]]
        tmp4 <- ICC.data[[i]][[3]]
        tmp5 <- ICC.data[[i]][[6]]
        subscore.t.adj[, tmp1.dim] <- subscore.t.adj[, tmp1.dim] + 
            cbind(tmp2[, 2], 1 - tmp2[, 3])
        subscore.obs[, tmp1.dim] <- subscore.obs[, tmp1.dim] + 
            cbind(tmp3[, 1], 1 - tmp3[, 1])
        subscore.t[, tmp1.dim] <- subscore.t[, tmp1.dim] + cbind(tmp4[, 
            2], 1 - tmp4[, 3])
        subscore.obs.adj[, tmp1.dim] <- subscore.obs.adj[, tmp1.dim] + 
            cbind(tmp5[, 1], 1 - tmp5[, 1])
    }
    subscore.t.adj <- cbind(ICC.data[[1]][[5]][, 1], subscore.t.adj)
    subscore.obs <- cbind(theta.score[, 1:nf], subscore.obs)
    subscore.obs.adj <- cbind(theta.score[, 1:nf], subscore.obs.adj)
    subscore.t <- cbind(ICC.data[[1]][[3]][, 1], subscore.t)
    colnames(subscore.t) <- colnames(subscore.t.adj) <- c("Theta", 
        paste0("F", 1:nf, ".subscore"))
    colnames(subscore.obs) <- colnames(subscore.obs.adj) <- c(paste0("F", 
        1:nf), paste0("F", 1:nf, ".subscore"))
    return(list(subscore.t.adj, subscore.obs, subscore.t, subscore.obs.adj))
}

#create TCC data
MUPP.2PL.5.fix.TCC <- TCC.MUPP.2PL(ICC.data=MUPP.2PL.5.fix.ICC, par.est=MUPP.2PL.5.par.fix, theta.score=MUPP.2PL.5.MAP.fix)

#########################################################################
#Function to draw TCC plots for both MUPP-2PL and Triplet-2PL 
#
#Description:
#Draw TCC plots using the TCC data output from the TCC.MUPP.2PL or TCC.Triplet.2PL function
#
#Usage:
#TCC.plot(TCC.data, plot.Factors, par.est, trait, PageNo = T)
#
#Arguments:
#TCC.data -- the TCC data output from the TCC.MUPP.2PL or TCC.Triplet.2PL function. 
#plot.Factors -- a numeric vector containing the trait number(s) for which the TCC(s) is drawn  
#par.est -- a Nitem-row matrix containing the following colnames: 
#		for the MUPP-2PL:
#           discrimination and intercept parameters for Statements 1 and 2 in each item: "a1"   "b1"   "a2"   "b2"  
#		trait number(1:Ntrait) measured by each statement in each item: "dim1" "dim2"
#		for the Triplet-2PL:
#           discrimination and intercept parameters for Statements 1-3 in each item: "a1"   "b1"   "a2"   "b2"  "a3"   "b3" 
#		trait number(1:Ntrait) measured by each statement in each item: "dim1" "dim2"  "dim3"
#
#trait -- a string vector containing the labels for all traits measured by a test
#PageNo -- logic; T = print page number at the bottom center, F = no page number
#adj.a.sign -- logic; T = changing negative statement(s) to positive one(s), F = no change 
#
#Return: 
#Output the TCC plot(s) for the selected trait(s). If adj.a.sign = T and any item measuring the plotted trait contains 
#statements with mixed keyed-directions, only the expected TCC is plotted. Otherwise, both the expected TCC and the empirical TCC with 95% CIs
#are plotted.   
##########################################################################

TCC.plot <-
function (TCC.data, plot.Factors, par.est, trait, PageNo = T, 
    adj.a.sign = F) 
{
    for (i in plot.Factors) {
        if (!adj.a.sign) {
            data <- TCC.data[[3]][, c(1, i + 1)]
            data.obs <- TCC.data[[2]][, c(paste0("F", i), paste0("F", 
                i, ".subscore"))]
            data <- data.frame(data, group = paste0("F", i, "-", 
                trait[i]))
            data.obs <- data.frame(data.obs, group = paste0("F", 
                i, "-", trait[i]))
            colnames(data)[2] <- "F"
            colnames(data.obs) <- c("Theta", "F", "group")
            plot1 <- xyplot(F ~ Theta | group, data.obs, panel = function(x, 
                y, ...) {
                panel.grid(h = -1, v = -1, col.line = grey(0.9))
                panel.smoother(x, y, ..., family = "g", n = 100, 
                  span = 0.9, lty = 2)
            }, ylab = list(label = "Expected Test Trait Score", 
                col = 1, cex = 1), strip = strip.custom(strip.names = c(F, 
                TRUE), par.strip.text = list(cex = 0.7)), key = list(text = list(c("Expected", 
                "Observed(with 95% CI)"), cex = 0.7, col = c(2, 
                4)), lines = list(lty = 1:2, col = c(2, 4), cex = 0.7), 
                space = "top", border = F, between = 1, adj = 0, 
                column = 2), layout = c(1, 1), as.table = T, 
                sub = if (PageNo) 
                  list(label = i, cex = 0.7)
                else "")
            plot2 <- xyplot(F ~ Theta | group, data, type = "l", 
                col = 2, grid = T, strip = strip.custom(strip.names = c(F, 
                  T), sep = " ", par.strip.text = list(cex = 0.7)), 
                xlab = list("Theta", cex = 1), ylab = list(label = "Expected Test Trait Score", 
                  cex = 1), scale = list(cex = 0.7, alternating = 1), 
                panel = function(x, y, ...) {
                  panel.xyplot(x, y, ...)
                }, layout = c(1, 1), as.table = T)
            print(plot1 + as.layer(plot2))
        }
        else {
            if ("dim3" %in% colnames(par.est)) 
                block.size <- 3
            else block.size <- 2
            if (block.size == 2) {
                par.est.tmp <- subset(par.est, par.est[, c("dim1")] == 
                  i | par.est[, c("dim2")] == i)
            }
            else {
                par.est.tmp <- subset(par.est, par.est[, c("dim1")] == 
                  i | par.est[, c("dim2")] == i | par.est[, c("dim3")] == 
                  i)
            }
            sign.tmp <- sign(par.est.tmp[, grep("a", colnames(par.est))])
            sign.tmp[sign.tmp == 0] <- 1
            sign.dir <- apply(sign.tmp, 1, sum)
            if ((any(sign.dir == 0) && block.size == 2) || (any(sign.dir %in% 
                c(-1, 1)) && block.size == 3)) {
                data <- TCC.data[[1]][, c(1, i + 1)]
                data <- data.frame(data, group = paste0("F", 
                  i, "-", trait[i]))
                colnames(data)[2] <- "F"
                plot2 <- xyplot(F ~ Theta | group, data, type = "l", 
                  col = 2, grid = T, strip = strip.custom(strip.names = c(F, 
                    T), sep = " ", par.strip.text = list(cex = 0.7)), 
                  xlab = list("Theta", cex = 1), ylab = list(label = "Expected Adjusted Test Trait Score", 
                    cex = 1), scale = list(cex = 0.7, alternating = 1), 
                  panel = function(x, y, ...) {
                    panel.xyplot(x, y, ...)
                  }, layout = c(1, 1), as.table = T)
                print(plot2)
            }
            else {
                if (all(sign.dir == 2) || all(sign.dir == 3)) 
                  ylab <- "Expected Test Trait Score"
                else ylab <- "Expected Adjusted Test Trait Score"
                data <- TCC.data[[1]][, c(1, i + 1)]
                data.obs <- TCC.data[[4]][, c(paste0("F", i), 
                  paste0("F", i, ".subscore"))]
                data <- data.frame(data, group = paste0("F", 
                  i, "-", trait[i]))
                data.obs <- data.frame(data.obs, group = paste0("F", 
                  i, "-", trait[i]))
                colnames(data)[2] <- "F"
                colnames(data.obs) <- c("Theta", "F", "group")
                plot1 <- xyplot(F ~ Theta | group, data.obs, 
                  panel = function(x, y, ...) {
                    panel.grid(h = -1, v = -1, col.line = grey(0.9))
                    panel.smoother(x, y, ..., family = "g", n = 100, 
                      span = 0.9, lty = 2)
                  }, ylab = list(label = ylab, col = 1, cex = 1), 
                  strip = strip.custom(strip.names = c(F, TRUE), 
                    par.strip.text = list(cex = 0.7)), key = list(text = list(c("Expected", 
                    "Observed(with 95% CI)"), cex = 0.7, col = c(2, 
                    4)), lines = list(lty = 1:2, col = c(2, 4), 
                    cex = 0.7), space = "top", border = F, between = 1, 
                    adj = 0, column = 2), layout = c(1, 1), as.table = T, 
                  sub = if (PageNo) 
                    list(label = i, cex = 0.7)
                  else "")
                plot2 <- xyplot(F ~ Theta | group, data, type = "l", 
                  col = 2, grid = T, strip = strip.custom(strip.names = c(F, 
                    T), sep = " ", par.strip.text = list(cex = 0.7)), 
                  xlab = list("Theta", cex = 1), ylab = list(label = ylab, 
                    cex = 1), scale = list(cex = 0.7, alternating = 1), 
                  panel = function(x, y, ...) {
                    panel.xyplot(x, y, ...)
                  }, layout = c(1, 1), as.table = T)
                print(plot1 + as.layer(plot2))
            }
        }
    }
}

#Trait lables
trait.label <- paste0("Trait ", 1:5)

#draw Trait 3's TCC plot to screen
TCC.plot(TCC.data=MUPP.2PL.5.fix.TCC, plot.Factors=3, par.est=MUPP.2PL.5.par.fix, trait=trait.label, PageNo = F)
TCC.plot(TCC.data=MUPP.2PL.5.fix.TCC, plot.Factors=3, par.est=MUPP.2PL.5.par.fix, trait=trait.label, PageNo = F, adj.a.sign=T)

#draw the TCC plots of all traits in a test and save to a PDF file

pdf("MUPP_2PL_5_fix_TCC_plot.pdf")

TCC.plot(TCC.data=MUPP.2PL.5.fix.TCC, plot.Factors=1:5, par.est=MUPP.2PL.5.par.fix, trait=trait.label, PageNo = T)

dev.off()

#draw the TCC plots of all traits in a test with negative statement(s) changed to positive one(s),
#and save to a PDF file

pdf("MUPP_2PL_5_fix_TCC_plot_adj.pdf")

TCC.plot(TCC.data=MUPP.2PL.5.fix.TCC, plot.Factors=1:5, par.est=MUPP.2PL.5.par.fix, trait=trait.label, PageNo = T, adj.a.sign=T)

dev.off()

#Simulate data based on Triplet-2PLM (Equation 9)
 
##################################################################
#Function to simulate data based on Triplet-2PLM 
#
#Description:
#Provide specifications to simulate Triplet-2PLM data. 
#
#Usage:
#simdata.Triplet.2PL(a=NULL, d=NULL, a.mean=NULL, a.std=NULL, d.mean=NULL, d.std=NULL, nf, nstate.per.f, theta.cor,
#dim, HPHN=T, Mix=T, mix.index=NULL, return.par=T, nsample=1000)
#
#Arguments:
#a -- vector of discrimination (a) parameters of the statements by the order in a test, provided only if "a" is fixed
#b -- vector of intercept (d) parameters of the statements by the order in a test, provided only if "d" is fixed
#a.mean, a.std -- discrimination (a) parameters' mean and std at normal (not log) scale, provided only if "a" is randomly draw by a lognormal distribution
#d.mean, d.std -- intercept (d) parameters' mean and std, provided only if "d" is randomly draw by a normal distribution
#nf --  number of factors/traits
#nstate.per.f --  number of statements measuring each factor/trait
#theta.cor -- theta/trait correlation matrix; theta/trait scores are assumed a standard multivariate normal distribution
#dim -- matrix of dimension/trait ID (1:nf) for each statement (column) in each item/block (row)
#HPHN -- T if the first half of all items/blocks are all positive statements and the second half are negative
#Mix -- T if also simulate a dataset with items/blocks containing mixed positive and negative statements
#mix.index -- a two-column index matrix to indicate which statements' "a" parameters should change signs: 
#	first column is the item/block number, and the second column is the statement number within an item
#return.par -- T if return item parameters, F otherwise

#Return:
#if Mix = F
#	if return.par = T, then return a list with three components: 
#		1. True item parameters: a (Nitem X 9) matrix with columns: a1 (discrimination for statement 1) d1 (slope for statement 1) dim1 (Trait ID measured by statement 1) a2 d2 dim2 a3 d3 dim3
#		2. True trait scores: a (Nsample X nf) matrix
#		3. simulated data: a (Nsample X Nitem) matrix with the following response coding (an item with Statements A, B, C in order) 
#			Statement Ranking		Response code 
#			A>B>C 				1
#			A>C>B 				2
#			B>A>C 				3
#			B>C>A 				4
#			C>A>B 				5
#			C>B>A 				6
#
#	else return a two-component list with #2 and #3
#else if return.par = T, then return a list with five components: #1-3 for nonmixed data, and #1 and #3 for mixed data.  
#	else return a three-component list with #2 and #3 for nonmixed data, and #3 for mixed data 
#		  
#############################################################   

simdata.Triplet.2PL <-
function (a = NULL, d = NULL, a.mean = NULL, a.std = NULL, d.mean = NULL, 
    d.std = NULL, nf, nstate.per.f, theta.cor, dim, HPHN = T, 
    Mix = F, mix.index = NULL, return.par = T, nsample = 1000) 
{
    fail <- T
    while (fail) {
        fail <- F
        if (is.null(a)) 
            a <- round(lognorm(nf * nstate.per.f, a.mean, a.std^2), 
                2)
        if (is.null(d)) 
            d <- round(rnorm(nf * nstate.per.f, d.mean, d.std), 
                2)
        a1 <- matrix(a, nc = 3)
        d1 <- matrix(d, nc = 3)
        if (HPHN) {
            a1[-(1:(nrow(a1)/2)), ] <- a1[-(1:(nrow(a1)/2)), 
                ] * (-1)
        }
        par <- cbind(a1, d1)[, c(1, 4, 2, 5, 3, 6)]
        theta <- matrix(round(rmvnorm(nsample, rep(0, nf), theta.cor), 
            2), ncol = nf)
        simdata <- matrix(NA, nr = nsample, nc = nrow(par))
        for (i1 in 1:nrow(par)) {
            prob <- P.Triplet.2PL.s(par[i1, ], theta[, dim[i1, 
                1]], theta[, dim[i1, 2]], theta[, dim[i1, 3]])
            for (i2 in 1:nsample) {
                simdata[i2, i1] <- sample(6, 1, replace = F, 
                  prob = prob[i2, ])
            }
            tmp <- table(simdata[, i1])
            if (any(tmp < 5) | length(tmp) != 6) {
                fail <- T
                break
            }
        }
    }
    colnames(simdata) <- paste0("Item", 1:nrow(par))
    par <- cbind(par, dim)[, c(1, 2, 7, 3, 4, 8, 5, 6, 9)]
    colnames(par) <- paste0(c("a", "d", "dim"), rep(1:3, each = 3))
    if (Mix) {
        a.mix <- a1
        for (i1 in (1:nrow(mix.index))) {
            a.mix[mix.index[i1, 1], mix.index[i1, 2]] <- a.mix[mix.index[i1, 
                1], mix.index[i1, 2]] * (-1)
        }
        par.mix <- cbind(a.mix, d1)[, c(1, 4, 2, 5, 3, 6)]
        simdata.mix <- matrix(NA, nr = nsample, nc = nrow(par.mix))
        for (i1 in 1:nrow(par.mix)) {
            prob.mix <- P.Triplet.2PL.s(par.mix[i1, ], theta[, 
                dim[i1, 1]], theta[, dim[i1, 2]], theta[, dim[i1, 
                3]])
            repeat {
                for (i2 in 1:nsample) {
                  simdata.mix[i2, i1] <- sample(6, 1, replace = F, 
                    prob = prob.mix[i2, ])
                }
                tmp <- table(simdata.mix[, i1])
                if (all(tmp >= 5) & length(tmp) == 6) {
                  break
                }
            }
        }
        colnames(simdata.mix) <- paste0("Item", 1:nrow(par))
        par.mix <- cbind(par.mix, dim)[, c(1, 2, 7, 3, 4, 8, 
            5, 6, 9)]
        colnames(par.mix) <- paste0(c("a", "d", "dim"), rep(1:3, 
            each = 3))
        if (return.par) 
            return(list(par, theta, simdata, par.mix, simdata.mix))
        else return(list(theta, simdata, simdata.mix))
    }
    else {
        if (return.par) 
            return(list(par, theta, simdata))
        else return(list(theta, simdata))
    }
}

lognorm <-
function (N, M, V) 
{
    Sdlog <- sqrt(log(V/M^2 + 1))
    Meanlog <- log(M) - Sdlog^2/2
    rlnorm(N, Meanlog, Sdlog)
}

#P.Triplet.2PL.s is the Triplet_2PL item response function (Equation 9 in the paper).
P.Triplet.2PL.s <-
function (par, theta1, theta2, theta3) 
{
    p1 <- 1/((1 + exp(par[3] * theta2 + par[4] - par[1] * theta1 - 
        par[2]) + exp(par[5] * theta3 + par[6] - par[1] * theta1 - 
        par[2])) * (1 + exp(par[5] * theta3 + par[6] - par[3] * 
        theta2 - par[4])))
    p2 <- 1/((1 + exp(par[5] * theta3 + par[6] - par[1] * theta1 - 
        par[2]) + exp(par[3] * theta2 + par[4] - par[1] * theta1 - 
        par[2])) * (1 + exp(-par[5] * theta3 - par[6] + par[3] * 
        theta2 + par[4])))
    p3 <- 1/((1 + exp(-par[3] * theta2 - par[4] + par[1] * theta1 + 
        par[2]) + exp(par[5] * theta3 + par[6] - par[3] * theta2 - 
        par[4])) * (1 + exp(par[5] * theta3 + par[6] - par[1] * 
        theta1 - par[2])))
    p4 <- 1/((1 + exp(-par[3] * theta2 - par[4] + par[5] * theta3 + 
        par[6]) + exp(par[1] * theta1 + par[2] - par[3] * theta2 - 
        par[4])) * (1 + exp(-par[5] * theta3 - par[6] + par[1] * 
        theta1 + par[2])))
    p5 <- 1/((1 + exp(-par[5] * theta3 - par[6] + par[1] * theta1 + 
        par[2]) + exp(par[3] * theta2 + par[4] - par[5] * theta3 - 
        par[6])) * (1 + exp(-par[1] * theta1 - par[2] + par[3] * 
        theta2 + par[4])))
    p6 <- 1/((1 + exp(par[3] * theta2 + par[4] - par[5] * theta3 - 
        par[6]) + exp(par[1] * theta1 + par[2] - par[5] * theta3 - 
        par[6])) * (1 + exp(-par[3] * theta2 - par[4] + par[1] * 
        theta1 + par[2])))
    cbind(p1, p2, p3, p4, p5, p6)
}

#There are 5 traits
Ntrait.3 <- 5

#Each trait measured by 6 statements
Nstate.per.trait.3 <- 6

#the Nitem X 3 matrix is the possible combinations of 3 out of 5 traits, containing Trait ID measured by each statement in each item
Dmatrix.3 <- t(combn(Ntrait.3,3)) 

#all intertrait correlations are set to 0.5
Trait.cor.3 <- matrix(0.5, nc=Ntrait.3, nr=Ntrait.3)
diag(Trait.cor.3) <-1

#Discrimination parameters are drawn from a lognormal distribution with mean=1 and std=0.5 (at normal scale)
#intercept parameters are drawn from a standard normal
#simulate 1000 test takers and 10 items (30 statements)
#change the signs of the discrimination parameters of all statements in the last 5 items to negative
#change the signs of the discrimination parameters of the first statement in Item 1, the first and second statements in Item 2 to negative 

sim.data.3 <- simdata.Triplet.2PL(a.mean=1, a.std=0.5, d.mean=0, d.std=1, nf=5, nstate.per.f=6, theta.cor=Trait.cor.3,
dim=Dmatrix.3, HPHN=T, Mix=T, mix.index =matrix(c(1,2,2,1,1,2),nc=2), nsample=1000)

#the current program requires no missing score category in every item.
#check minimum number of responses in a score category >0 
min(apply(sim.data.3[[5]], 2, table))


#create new item type, "Triplet_2PL", in mirt
name <- 'Triplet_2PL'
par0 <- rep(c(1,0), Ntrait.3)
names(par0)<- c(paste0(c("a", "d"), rep(1:Ntrait.3, each=2)))
est <- rep(F, 2*Ntrait.3)

#P.Triplet.2PL.s is the Triplet_2PL item response function (Equation 9 in the paper)for the new mirt item "Triplet_2PL".

P.Triplet.2PL <-
function (par, Theta, ncat, est) 
{
    theta <- Theta[, est, drop = F]
    est1 <- c((est - 1) * 2 + 1, (est - 1) * 2 + 2)
    est1 <- sort(est1)
    par1 <- par[est1]
    P.Triplet.2PL.s(par1, theta[, 1], theta[, 2], theta[, 3])
}

#P.Triplet.2PL.gr is the first derivative function of the Triplet_2PL item response function 
#P.Triplet.2PL.der1 is the function output from the symbolic derivative function, Deriv.
P.Triplet.2PL.gr <-
function (x, Theta) 
{
    P <- P.Triplet.2PL(x@par, Theta, ncat, x@userdata[[1]])
    nfactor <- ncol(Theta)
    ThetaLength <- nrow(Theta)
    r_P <- x@dat/P
    est <- x@userdata[[1]]
    est1 <- c((est - 1) * 2 + 1, (est - 1) * 2 + 2)
    est1 <- sort(est1)
    dp1 <- array(P.Triplet.2PL.der1(x@par[est1], Theta[, est[1]], 
        Theta[, est[2]], Theta[, est[3]]), c(ThetaLength, x@ncat, 
        length(est1)))
    grad <- numeric(length(x@par))
    j <- 1
    for (i in est1) {
        grad[i] <- sum(r_P * dp1[, , j])
        j <- j + 1
    }
    grad
}

P.Triplet.2PL.der1 <-
function (par, theta1, theta2, theta3) 
{
    .e1 <- par[2] + theta1 * par[1]
    .e2 <- par[4] + theta2 * par[3]
    .e3 <- par[6] + theta3 * par[5]
    .e4 <- exp(.e1 - .e2)
    .e5 <- exp(.e1 - .e3)
    .e6 <- exp(.e2 - .e1)
    .e7 <- exp(.e2 - .e3)
    .e8 <- exp(.e3 - .e1)
    .e9 <- exp(.e3 - .e2)
    .e10 <- 1 + .e4
    .e11 <- 1 + .e5
    .e12 <- 1 + .e6
    .e13 <- .e10 + .e9
    .e14 <- .e11 + .e7
    .e15 <- .e12 + .e8
    .e16 <- 1 + .e7
    .e17 <- 1 + .e8
    .e18 <- 1 + .e9
    .e19 <- (.e13 * .e11)^2
    .e20 <- (.e13 * .e17)^2
    .e21 <- (.e10 * .e14)^2
    .e22 <- (.e14 * .e12)^2
    .e23 <- (.e15 * .e16)^2
    .e24 <- (.e15 * .e18)^2
    .e25 <- .e13 * .e5
    .e26 <- .e13 * .e8
    .e27 <- .e14 * .e4
    .e28 <- .e14 * .e6
    .e29 <- .e15 * .e7
    .e30 <- .e15 * .e9
    .e31 <- .e4 + .e9
    .e32 <- .e5 + .e7
    .e33 <- .e6 + .e8
    .e34 <- .e25 + .e11 * .e4
    .e35 <- .e26 + .e17 * .e9
    .e37 <- .e10 * .e5 + .e27
    .e39 <- .e10 * .e7 - .e27
    .e40 <- .e28 + .e12 * .e7
    .e42 <- .e11 * .e9 - .e25
    .e43 <- .e29 + .e16 * .e6
    .e44 <- .e30 + .e18 * .e8
    .e46 <- .e12 * .e5 - .e28
    .e48 <- .e16 * .e8 - .e29
    .e50 <- .e17 * .e4 - .e26
    .e52 <- .e18 * .e6 - .e30
    c(par_1 = cbind(theta1 * .e18 * .e33/.e24, theta1 * .e16 * 
        .e33/.e23, -(theta1 * .e50/.e20), -(theta1 * .e34/.e19), 
        -(theta1 * .e46/.e22), -(theta1 * .e37/.e21)), par_2 = cbind(.e18 * 
        .e33/.e24, .e16 * .e33/.e23, -(.e50/.e20), -(.e34/.e19), 
        -(.e46/.e22), -(.e37/.e21)), par_3 = cbind(-(theta2 * 
        .e52/.e24), -(theta2 * .e43/.e23), theta2 * .e17 * .e31/.e20, 
        theta2 * .e11 * .e31/.e19, -(theta2 * .e40/.e22), -(theta2 * 
            .e39/.e21)), par_4 = cbind(-(.e52/.e24), -(.e43/.e23), 
        .e17 * .e31/.e20, .e11 * .e31/.e19, -(.e40/.e22), -(.e39/.e21)), 
        par_5 = cbind(-(theta3 * .e44/.e24), -(theta3 * .e48/.e23), 
            -(theta3 * .e35/.e20), -(theta3 * .e42/.e19), theta3 * 
                .e12 * .e32/.e22, theta3 * .e10 * .e32/.e21), 
        par_6 = cbind(-(.e44/.e24), -(.e48/.e23), -(.e35/.e20), 
            -(.e42/.e19), .e12 * .e32/.e22, .e10 * .e32/.e21))
}

#P.Triplet.2PL.hss is the second derivative function of the Triplet_2PL item response function, 
#P.Triplet.2PL.der2 is the function output from the symbolic derivative function, Deriv.
P.Triplet.2PL.hss <-
function (x, Theta) 
{
    P <- P.Triplet.2PL(x@par, Theta, ncat, x@userdata[[1]])
    nfactor <- ncol(Theta)
    xLength <- length(x@par)
    ThetaLength <- nrow(Theta)
    est <- x@userdata[[1]]
    est1 <- c((est - 1) * 2 + 1, (est - 1) * 2 + 2)
    est1 <- sort(est1)
    npar <- length(est1)
    dp1 <- array(P.Triplet.2PL.der1(x@par[est1], Theta[, est[1]], 
        Theta[, est[2]], Theta[, est[3]]), c(ThetaLength, x@ncat, 
        npar))
    dp2 <- array(P.Triplet.2PL.der2(x@par[est1], Theta[, est[1]], 
        Theta[, est[2]], Theta[, est[3]]), c(ThetaLength, x@ncat, 
        npar, npar))
    H <- matrix(0, xLength, xLength)
    H1 <- matrix(NA, npar, npar)
    P2 <- P^2
    for (i in 1:npar) {
        for (j in i:npar) {
            H1[i, j] <- sum(x@dat * dp2[, , i, j]/P + x@dat * 
                dp1[, , i] * (-dp1[, , j]/P2))
            H1[j, i] <- H1[i, j]
        }
    }
    H[est1, est1] <- H1
    H
}

P.Triplet.2PL.der2 <-
function (par, theta1, theta2, theta3) 
{
    .e1 <- par[6] + theta3 * par[5]
    .e2 <- par[4] + theta2 * par[3]
    .e3 <- par[2] + theta1 * par[1]
    .e4 <- exp(.e3 - .e1)
    .e5 <- exp(.e2 - .e1)
    .e6 <- exp(.e1 - .e2)
    .e7 <- exp(.e3 - .e2)
    .e8 <- exp(.e2 - .e3)
    .e9 <- exp(.e1 - .e3)
    .e10 <- 1 + .e7
    .e11 <- 1 + .e4
    .e12 <- 1 + .e8
    .e13 <- .e10 + .e6
    .e14 <- .e11 + .e5
    .e15 <- .e12 + .e9
    .e16 <- 1 + .e5
    .e17 <- 1 + .e9
    .e18 <- 1 + .e6
    .e19 <- .e7 + .e6
    .e20 <- .e4 + .e5
    .e21 <- .e8 + .e9
    .e22 <- (.e13 * .e11)^2
    .e23 <- (.e13 * .e17)^2
    .e24 <- (.e10 * .e14)^2
    .e25 <- (.e14 * .e12)^2
    .e26 <- (.e15 * .e16)^2
    .e27 <- (.e15 * .e18)^2
    .e28 <- .e13 * .e4
    .e29 <- .e13 * .e9
    .e30 <- .e14 * .e7
    .e31 <- .e14 * .e8
    .e32 <- .e15 * .e5
    .e33 <- .e15 * .e6
    .e34 <- .e10 * .e4
    .e35 <- .e11 * .e7
    .e36 <- .e12 * .e5
    .e37 <- .e16 * .e8
    .e38 <- .e17 * .e6
    .e39 <- .e18 * .e9
    .e40 <- .e28 + .e35
    .e41 <- .e29 + .e38
    .e42 <- .e34 + .e30
    .e44 <- .e10 * .e5 - .e30
    .e45 <- .e31 + .e36
    .e47 <- .e11 * .e6 - .e28
    .e48 <- .e32 + .e37
    .e49 <- .e33 + .e39
    .e51 <- .e12 * .e4 - .e31
    .e53 <- .e16 * .e9 - .e32
    .e55 <- .e17 * .e7 - .e29
    .e57 <- .e18 * .e8 - .e33
    .e58 <- .e19 * .e4
    .e59 <- .e19 * .e9
    .e60 <- .e20 * .e8
    .e61 <- .e21 * .e5
    .e62 <- .e21 * .e6
    .e63 <- .e7 * .e20
    .e64 <- (1 + 2 * .e4) * .e7
    .e65 <- (1 + 2 * .e5) * .e8
    .e66 <- (1 + 2 * .e6) * .e9
    .e67 <- 1 + 2 * .e7
    .e68 <- 1 + 2 * .e8
    .e69 <- 1 + 2 * .e9
    .e70 <- theta1 * theta2
    .e71 <- theta1 * theta3
    .e72 <- theta2 * theta3
    .e146 <- .e67 * .e4 + 2 * (.e40 * .e47 * .e13 * .e11/.e22)
    .e147 <- .e64 + 2 * (.e42 * .e44 * .e10 * .e14/.e24)
    .e149 <- .e68 * .e5 + 2 * (.e48 * .e53 * .e15 * .e16/.e26)
    .e150 <- .e65 + 2 * (.e45 * .e51 * .e14 * .e12/.e25)
    .e152 <- .e69 * .e6 + 2 * (.e49 * .e57 * .e15 * .e18/.e27)
    .e153 <- .e66 + 2 * (.e41 * .e55 * .e13 * .e17/.e23)
    .e154 <- 2 * (.e40 * .e13 * .e11 * .e19/.e22)
    .e155 <- 2 * (.e41 * .e13 * .e17 * .e19/.e23)
    .e156 <- 2 * (.e42 * .e10 * .e14 * .e20/.e24)
    .e157 <- 2 * (.e44 * .e10 * .e14 * .e20/.e24)
    .e158 <- 2 * (.e45 * .e14 * .e12 * .e20/.e25)
    .e159 <- 2 * (.e47 * .e13 * .e11 * .e19/.e22)
    .e160 <- 2 * (.e48 * .e15 * .e16 * .e21/.e26)
    .e161 <- 2 * (.e49 * .e15 * .e18 * .e21/.e27)
    .e162 <- 2 * (.e51 * .e14 * .e12 * .e20/.e25)
    .e163 <- 2 * (.e53 * .e15 * .e16 * .e21/.e26)
    .e164 <- 2 * (.e55 * .e13 * .e17 * .e19/.e23)
    .e165 <- 2 * (.e57 * .e15 * .e18 * .e21/.e27)
    .e166 <- theta1^2
    .e167 <- theta2^2
    .e168 <- theta3^2
    .e216 <- (.e67 + .e6) * .e4 + .e64 - 2 * (.e40^2 * .e13 * 
        .e11/.e22)
    .e219 <- (.e68 + .e9) * .e5 + .e65 - 2 * (.e48^2 * .e15 * 
        .e16/.e26)
    .e222 <- (.e69 + .e8) * .e6 + .e66 - 2 * (.e49^2 * .e15 * 
        .e18/.e27)
    .e224 <- (1 + 3 * .e4 + .e5) * .e7 + .e10 * (.e4 - 2 * (.e42^2 * 
        .e14/.e24))
    .e226 <- (1 + 3 * .e5 + .e4) * .e8 + .e12 * (.e5 - 2 * (.e45^2 * 
        .e14/.e25))
    .e228 <- (1 + 3 * .e6 + .e7) * .e9 + .e17 * (.e6 - 2 * (.e41^2 * 
        .e13/.e23))
    .e236 <- .e10 * (.e156 - .e4) - .e63
    .e238 <- .e10 * (.e157 - .e5) + .e63
    .e240 <- .e10 * (.e4 - .e156) + .e63
    .e242 <- .e10 * (.e5 - .e157) - .e63
    .e244 <- .e34 + .e6 - 2 * (.e47^2 * .e13 * .e11/.e22)
    .e253 <- .e11 * (.e154 - .e7) - .e58
    .e255 <- .e11 * (.e159 - .e6) + .e58
    .e257 <- .e11 * (.e7 - .e154) + .e58
    .e259 <- .e11 * (.e6 - .e159) - .e58
    .e261 <- .e35 + .e5 - 2 * (.e44^2 * .e10 * .e14/.e24)
    .e269 <- .e12 * (.e158 - .e5) - .e60
    .e271 <- .e12 * (.e162 - .e4) + .e60
    .e273 <- .e12 * (.e4 - .e162) - .e60
    .e275 <- .e12 * (.e5 - .e158) + .e60
    .e277 <- .e36 + .e9 - 2 * (.e53^2 * .e15 * .e16/.e26)
    .e279 <- .e16 * (.e160 - .e8) - .e61
    .e281 <- .e16 * (.e163 - .e9) + .e61
    .e283 <- .e16 * (.e8 - .e160) + .e61
    .e285 <- .e16 * (.e9 - .e163) - .e61
    .e287 <- .e37 + .e4 - 2 * (.e51^2 * .e14 * .e12/.e25)
    .e289 <- .e17 * (.e155 - .e6) - .e59
    .e291 <- .e17 * (.e164 - .e7) + .e59
    .e293 <- .e17 * (.e7 - .e164) - .e59
    .e295 <- .e17 * (.e6 - .e155) + .e59
    .e297 <- .e38 + .e8 - 2 * (.e57^2 * .e15 * .e18/.e27)
    .e299 <- .e18 * (.e161 - .e9) - .e62
    .e301 <- .e18 * (.e165 - .e8) + .e62
    .e303 <- .e18 * (.e8 - .e165) - .e62
    .e305 <- .e18 * (.e9 - .e161) + .e62
    .e307 <- .e39 + .e7 - 2 * (.e55^2 * .e13 * .e17/.e23)
    .e309 <- 2 * (.e13 * .e11^2 * .e19/.e22) - 1
    .e311 <- 2 * (.e13 * .e17^2 * .e19/.e23) - 1
    .e313 <- 2 * (.e10^2 * .e14 * .e20/.e24) - 1
    .e315 <- 2 * (.e14 * .e12^2 * .e20/.e25) - 1
    .e317 <- 2 * (.e15 * .e16^2 * .e21/.e26) - 1
    .e319 <- 2 * (.e15 * .e18^2 * .e21/.e27) - 1
    .e332 <- .e146/.e22
    .e333 <- .e147/.e24
    .e334 <- .e149/.e26
    .e335 <- .e150/.e25
    .e336 <- .e152/.e27
    .e337 <- .e153/.e23
    .e338 <- cbind(-(theta2 * .e297/.e27), -(theta2 * .e219/.e26), 
        theta2 * .e17 * .e311 * .e19/.e23, theta2 * .e11 * .e309 * 
            .e19/.e22, -(theta2 * .e226/.e25), -(theta2 * .e261/.e24))
    .e339 <- cbind(-(theta3 * .e222/.e27), -(theta3 * .e277/.e26), 
        -(theta3 * .e228/.e23), -(theta3 * .e244/.e22), theta3 * 
            .e12 * .e315 * .e20/.e25, theta3 * .e10 * .e313 * 
            .e20/.e24)
    .e340 <- cbind(theta1 * .e18 * .e319 * .e21/.e27, theta1 * 
        .e16 * .e317 * .e21/.e26, -(theta1 * .e307/.e23), -(theta1 * 
        .e216/.e22), -(theta1 * .e287/.e25), -(theta1 * .e224/.e24))
    .e342 <- theta1 * .e146/.e22
    .e344 <- theta1 * .e147/.e24
    .e346 <- theta1 * .e150/.e25
    .e348 <- theta1 * .e153/.e23
    .e350 <- .e70 * .e147/.e24
    .e352 <- .e70 * .e150/.e25
    .e354 <- .e71 * .e146/.e22
    .e356 <- .e71 * .e153/.e23
    .e358 <- theta2 * .e147/.e24
    .e360 <- theta2 * .e149/.e26
    .e362 <- theta2 * .e150/.e25
    .e364 <- theta2 * .e152/.e27
    .e366 <- .e72 * .e149/.e26
    .e368 <- .e72 * .e152/.e27
    .e370 <- theta3 * .e146/.e22
    .e372 <- theta3 * .e149/.e26
    .e374 <- theta3 * .e152/.e27
    .e376 <- theta3 * .e153/.e23
    c(par_1 = c(par_1 = cbind(.e166 * .e18 * .e319 * .e21/.e27, 
        .e166 * .e16 * .e317 * .e21/.e26, -(.e166 * .e307/.e23), 
        -(.e166 * .e216/.e22), -(.e166 * .e287/.e25), -(.e166 * 
            .e224/.e24)), par_2 = .e340, par_3 = cbind(-(.e70 * 
        .e301/.e27), -(.e70 * .e279/.e26), .e70 * .e293/.e23, 
        .e70 * .e257/.e22, .e352, .e350), par_4 = cbind(-(theta1 * 
        .e301/.e27), -(theta1 * .e279/.e26), theta1 * .e293/.e23, 
        theta1 * .e257/.e22, .e346, .e344), par_5 = cbind(-(.e71 * 
        .e299/.e27), -(.e71 * .e281/.e26), .e356, .e354, .e71 * 
        .e273/.e25, .e71 * .e240/.e24), par_6 = cbind(-(theta1 * 
        .e299/.e27), -(theta1 * .e281/.e26), .e348, .e342, theta1 * 
        .e273/.e25, theta1 * .e240/.e24)), par_2 = c(par_1 = .e340, 
        par_2 = cbind(.e18 * .e319 * .e21/.e27, .e16 * .e317 * 
            .e21/.e26, -(.e307/.e23), -(.e216/.e22), -(.e287/.e25), 
            -(.e224/.e24)), par_3 = cbind(-(theta2 * .e301/.e27), 
            -(theta2 * .e279/.e26), theta2 * .e293/.e23, theta2 * 
                .e257/.e22, .e362, .e358), par_4 = cbind(-(.e301/.e27), 
            -(.e279/.e26), .e293/.e23, .e257/.e22, .e335, .e333), 
        par_5 = cbind(-(theta3 * .e299/.e27), -(theta3 * .e281/.e26), 
            .e376, .e370, theta3 * .e273/.e25, theta3 * .e240/.e24), 
        par_6 = cbind(-(.e299/.e27), -(.e281/.e26), .e337, .e332, 
            .e273/.e25, .e240/.e24)), par_3 = c(par_1 = cbind(.e70 * 
        .e303/.e27, .e70 * .e283/.e26, -(.e70 * .e291/.e23), 
        -(.e70 * .e253/.e22), .e352, .e350), par_2 = cbind(theta2 * 
        .e303/.e27, theta2 * .e283/.e26, -(theta2 * .e291/.e23), 
        -(theta2 * .e253/.e22), .e362, .e358), par_3 = cbind(-(.e167 * 
        .e297/.e27), -(.e167 * .e219/.e26), .e167 * .e17 * .e311 * 
        .e19/.e23, .e167 * .e11 * .e309 * .e19/.e22, -(.e167 * 
        .e226/.e25), -(.e167 * .e261/.e24)), par_4 = .e338, par_5 = cbind(.e368, 
        .e366, -(.e72 * .e289/.e23), -(.e72 * .e255/.e22), .e72 * 
            .e275/.e25, .e72 * .e242/.e24), par_6 = cbind(.e364, 
        .e360, -(theta2 * .e289/.e23), -(theta2 * .e255/.e22), 
        theta2 * .e275/.e25, theta2 * .e242/.e24)), par_4 = c(par_1 = cbind(theta1 * 
        .e303/.e27, theta1 * .e283/.e26, -(theta1 * .e291/.e23), 
        -(theta1 * .e253/.e22), .e346, .e344), par_2 = cbind(.e303/.e27, 
        .e283/.e26, -(.e291/.e23), -(.e253/.e22), .e335, .e333), 
        par_3 = .e338, par_4 = cbind(-(.e297/.e27), -(.e219/.e26), 
            .e17 * .e311 * .e19/.e23, .e11 * .e309 * .e19/.e22, 
            -(.e226/.e25), -(.e261/.e24)), par_5 = cbind(.e374, 
            .e372, -(theta3 * .e289/.e23), -(theta3 * .e255/.e22), 
            theta3 * .e275/.e25, theta3 * .e242/.e24), par_6 = cbind(.e336, 
            .e334, -(.e289/.e23), -(.e255/.e22), .e275/.e25, 
            .e242/.e24)), par_5 = c(par_1 = cbind(.e71 * .e305/.e27, 
        .e71 * .e285/.e26, .e356, .e354, -(.e71 * .e271/.e25), 
        -(.e71 * .e236/.e24)), par_2 = cbind(theta3 * .e305/.e27, 
        theta3 * .e285/.e26, .e376, .e370, -(theta3 * .e271/.e25), 
        -(theta3 * .e236/.e24)), par_3 = cbind(.e368, .e366, 
        .e72 * .e295/.e23, .e72 * .e259/.e22, -(.e72 * .e269/.e25), 
        -(.e72 * .e238/.e24)), par_4 = cbind(.e374, .e372, theta3 * 
        .e295/.e23, theta3 * .e259/.e22, -(theta3 * .e269/.e25), 
        -(theta3 * .e238/.e24)), par_5 = cbind(-(.e168 * .e222/.e27), 
        -(.e168 * .e277/.e26), -(.e168 * .e228/.e23), -(.e168 * 
            .e244/.e22), .e168 * .e12 * .e315 * .e20/.e25, .e168 * 
            .e10 * .e313 * .e20/.e24), par_6 = .e339), par_6 = c(par_1 = cbind(theta1 * 
        .e305/.e27, theta1 * .e285/.e26, .e348, .e342, -(theta1 * 
        .e271/.e25), -(theta1 * .e236/.e24)), par_2 = cbind(.e305/.e27, 
        .e285/.e26, .e337, .e332, -(.e271/.e25), -(.e236/.e24)), 
        par_3 = cbind(.e364, .e360, theta2 * .e295/.e23, theta2 * 
            .e259/.e22, -(theta2 * .e269/.e25), -(theta2 * .e238/.e24)), 
        par_4 = cbind(.e336, .e334, .e295/.e23, .e259/.e22, -(.e269/.e25), 
            -(.e238/.e24)), par_5 = .e339, par_6 = cbind(-(.e222/.e27), 
            -(.e277/.e26), -(.e228/.e23), -(.e244/.e22), .e12 * 
                .e315 * .e20/.e25, .e10 * .e313 * .e20/.e24)))
}

#create new item type, "Triplet_2PL"
Item.Triplet.2PL.5 <- createItem(name, par=par0, est=est, P=P.Triplet.2PL, gr=P.Triplet.2PL.gr, hss=P.Triplet.2PL.hss)

#create model object for mirt function

Q <- matrix(0, nr=nrow(Dmatrix.3), ncol=Ntrait.3, dimnames = list(NULL, paste0('F', 1:Ntrait.3)))
for (i in 1:nrow(Q)){
	Q[i,Dmatrix.3[i,]] <-1
}
COV <- matrix(T, nrow=Ntrait.3, nc=Ntrait.3)
diag(COV) <- F
Model.Triplet.2PL.5 <- mirt.model(Q, COV=COV)

#create user data used in mirt function: a Nitem list with each component being one row of Dmatrix.3
user.data.3 <- vector('list', nrow(Dmatrix.3))
for ( i in 1:nrow(Dmatrix.3)){
	user.data.3[[i]] <- Dmatrix.3[i,]
}

#fix/free item parameters for a model with direct item parameter estimation: 
#	1. fix the intercept of the first statement in each item to the true value
#	2. fix intertrait correlation matrix to the true values
#set initial values: -1 for negative discrimination, 1 for positive discrimination, and 0 for intercept
 
tmp <- mirt(sim.data.3[[5]], Model.Triplet.2PL.5, 'Triplet_2PL', customItems=list(Triplet_2PL=Item.Triplet.2PL.5), pars = 'values')
par <- sim.data.3[[4]]
tmp1 <- NULL
for ( i in 1:nrow(par)){		
	tmp1 <- rbind(tmp1, data.frame(item=paste0("Item", i), name=paste0(c("a", "d"), rep(Dmatrix.3[i,], each=2)), value= c(sign(par[i,"a1"]),par[i,"d1"],sign(par[i,"a2"]),0,sign(par[i,"a3"]),0),est=c(T, F, T, T, T, T)))		
}

tmp3 <- merge (tmp, tmp1, by=c("item", "name"), all.x=T, sort=F)
tmp3 <- tmp3[order(tmp3$parnum),]

tmp3$value <- ifelse(is.na(tmp3$value.y), tmp3$value.x, tmp3$value.y)
tmp3$est <- ifelse(is.na(tmp3$est.y), tmp3$est.x, tmp3$est.y)
tmp4 <- tmp3[, colnames(tmp)]

#skip the following two lines if you do not want to fix the correlation matrix to the true values
tmp4[grepl("COV", tmp4$name, fixed=T), ]$value <- Trait.cor.3[lower.tri(Trait.cor.3,T)]
tmp4[grepl("COV", tmp4$name, fixed=T), ]$est <-F

sv.Triplet.2PL.5.est <- tmp4

#fix/free item parameters for a model with fixed item parameters: 
#	1. fix the intercept of the first statement in each item to the true value
#	2. fix intertrait correlation matrix to the true values
#set initial values of all item parameters to the true values:
 
tmp <- mirt(sim.data.3[[5]], Model.Triplet.2PL.5, 'Triplet_2PL', customItems=list(Triplet_2PL=Item.Triplet.2PL.5), pars = 'values')
par <- sim.data.3[[4]]
tmp1 <- NULL
for ( i in 1:nrow(par)){		
	tmp1 <- rbind(tmp1, data.frame(item=paste0("Item", i), name=paste0(c("a", "d"), rep(Dmatrix.3[i,], each=2)), value= par[i,paste0(c("a", "d"), rep(1:3, each=2))],est=c(T, F, T, T, T, T)))	
}

tmp3 <- merge (tmp, tmp1, by=c("item", "name"), all.x=T, sort=F)
tmp3 <- tmp3[order(tmp3$parnum),]

tmp3$value <- ifelse(is.na(tmp3$value.y), tmp3$value.x, tmp3$value.y)
tmp3$est <- ifelse(is.na(tmp3$est.y), tmp3$est.x, tmp3$est.y)
tmp4 <- tmp3[, colnames(tmp)]

#skip the following two lines if you do not want to fix the correlation matrix to the true values
tmp4[grepl("COV", tmp4$name, fixed=T), ]$value <- Trait.cor.3[lower.tri(Trait.cor.3,T)]
tmp4[grepl("COV", tmp4$name, fixed=T), ]$est <-F

sv.Triplet.2PL.5.fix <- tmp4

#run mirt() to estimate item parameters using "QMCEM" method; 'MHRM' method produces better estimates 
#but takes a much longer time
 
  
Triplet.2PL.5.est <-mirt(sim.data.3[[5]], Model.Triplet.2PL.5, 'Triplet_2PL', 
customItems=list(Triplet_2PL=Item.Triplet.2PL.5), pars =sv.Triplet.2PL.5.est, 
method='QMCEM',large=T, SE=T, customItemsData=user.data.3)

Triplet.2PL.5.fix <-mirt(sim.data.3[[5]], Model.Triplet.2PL.5, 'Triplet_2PL', 
customItems=list(Triplet_2PL=Item.Triplet.2PL.5), pars =sv.Triplet.2PL.5.fix, 
method='QMCEM',large=T, SE=F, customItemsData=user.data.3, TOL=NA)


##################################################################
#Function to extract item parameter estimates of Triplet_2PL 
#
#Description:
#Given a mirt object of Triplet_2PL estimation, extract item parameter estimates to a desired form. 
#
#Usage:
#extract.par.Triplet.2PL(mirt.obj, dim)
#
#Arguments:
#mirt.obj -- an object returned by the mirt function estimating Triplet_2PL
#dim -- a Nitem-row matrix including trait ID measured by each statement in each item, with column 
#names as "dim1", "dim2" and "dim3". 
#
#Return:
#a Nitem X 9 matrix including item parameter estimates and trait ID for each statement in each item
#with column names a1, b1, a2, b2, a3, b3, dim1, dim2, dim3. 
#
#Note:
#Assume that in "Triplet_2PL" item type, item parameters associated with a trait, for example, trait 1,
#is labeled as "a1" for discrimination and "d1" for intercept. 
#In the output file, intercept is labeled by "b", a different letter.   
#
##########################################################################

extract.par.Triplet.2PL <-
function (mirt.obj, dim) 
{
    tmp <- coef(mirt.obj, simplify = T)
    tmp1 <- round(extract.par.est3(tmp[[1]], dim[, c("dim1", 
        "dim2", "dim3")]), 2)
    colnames(tmp1) <- paste0(c("a", "b"), rep(1:3, each = 2))
    cbind(tmp1, dim[, c("dim1", "dim2", "dim3")])
}
extract.par.est3 <-
function (par, dim) 
{
    tmp1 <- NULL
    for (i in 1:nrow(par)) {
        tmp1 <- rbind(tmp1, par[i, paste0(c("a", "d"), rep(dim[i, 
            ], each = 2))])
    }
    tmp1
}

Triplet.2PL.5.par.ets <- extract.par.Triplet.2PL(Triplet.2PL.5.ets, sim.data.3[[1]])
Triplet.2PL.5.par.fix <- extract.par.Triplet.2PL(Triplet.2PL.5.fix, sim.data.3[[1]])


#estimate MAP trait scores
Triplet.2PL.5.MAP.est <- fscores(Triplet.2PL.5.est, method="MAP", QMC=T, full.scores.SE=T) 
Triplet.2PL.5.MAP.fix <- fscores(Triplet.2PL.5.fix, method="MAP", QMC=T, full.scores.SE=T) 

#item fit: Drasgow, Levine, & Williams (1985) Zh statistic
Triplet.2PL.5.ZH.est <- itemfit( Triplet.2PL.5.est, fit_stats ="Zh", na.rm=T, Theta=as.matrix( Triplet.2PL.5.MAP.est[,1:Ntrait.2]), QMC=TRUE)
Triplet.2PL.5.ZH.fix <- itemfit( Triplet.2PL.5.fix, fit_stats ="Zh", na.rm=T, Theta=as.matrix( Triplet.2PL.5.MAP.fix[,1:Ntrait.2]), QMC=TRUE)

#########################################################################
#Function to create Triplet-2PL ICC data that are used for drawing ICC plots 
#
#Description:
#Create observed and expected ICC data for one item.
#
#Usage:
#ICC.Triplet.2PL(item, model.est, dimcor, par.est, theta.up, theta.low, 
#    theta.interval = NULL, nquad = NULL, data, theta.score)
#
#Arguments:
#item -- numeric, item sequence number (e.g., 1) of the item whose ICC is calculated
#model.est -- an object returned by the mirt function estimating Triplet-2PL
#dimcor -- the covariance/correlation matrix of trait scores
#par.est -- a Nitem-row matrix containing the following colnames: 
#           discrimination and intercept parameters for Statements 1-3 in each item: "a1"   "b1"   "a2"   "b2"  "a3"   "b3" 
#		trait number(1:Ntrait) measured by each statement in each item: "dim1" "dim2"  "dim3"
#theta.up -- the upper bound of the theta range
#theta.low -- the lower bound of the theta range 
#theta.interval -- equal spaced interval in the theta range
#nquad -- number of theta points in the theta range; if both theta.interval and nquad are provided, theta.interval overwrites nquad
#data -- a matrix (Nexaminee X Nitem)containing item score (1=select the first statement,0=select the second statement)
#theta.score -- a matrix (Nexaminee X Ntrait) containing trait score(s) for each examinee
#
#Return:
#A list with 6 elements. 
#
#1. matrix  :  item true trait scores at a theta vector   
#  ..$ F1       : trait point measured by the first statement
#  ..$ F2       : trait point measured by the second statement
#  ..$ F3       : trait point measured by the third statement
#  ..$ weight   : weight (probability) of the theta vector (F1, F2, F3)
#  ..$ F1.subscore : expected item trait scores measured by the first statement at (F1, F2, F3)
#  ..$ F2.subscore : expected item trait scores measured by the second statement at (F1, F2, F3)
#  ..$ F3.subscore : expected item trait scores measured by the third statement at (F1, F2, F3)
#  ..$ F1.subscore.w : F1.subscore*weight
#  ..$ F2.subscore.w : F2.subscore*weight
#  ..$ F3.subscore.w : F3.subscore*weight
#
#2. matrix  : item true trait scores for each examinee 
#  ..$ score  : observed item score
#  ..$ F1.subscore  : observed item trait score measured by the first statement
#  ..$ F2.subscore  : observed item trait score measured by the second statement
#  ..$ F3.subscore  : observed item trait score measured by the third statement
#  ..$ F1.subscore.t  : expected item trait score measured by the first statement
#  ..$ F2.subscore.t  : expected item trait score measured by the second statement
#  ..$ F3.subscore.t  : expected item trait score measured by the third statement
#  ..$ F1      : trait score estimate by the first statement
#  ..$ F2      : trait score estimate by the second statement
#  ..$ F3      : trait score estimate by the third statement
#
#3. data.frame  : expected (marginal) item true trait score at a theta point 
#  ..$ Theta  : a theta point
#  ..$ F1.F1.subscore  : expected item true trait score measured by the first statement at a point of the trait measured by the first statement
#  ..$ F1.F2.subscore  : expected item true trait score measured by the second statement at a point of the trait measured by the first statement
#  ..$ F1.F3.subscore  : expected item true trait score measured by the third statement at a point of the trait measured by the first statement
#  ..$ F2.F1.subscore  : expected item true trait score measured by the first statement at a point of the trait measured by the second statement
#  ..$ F2.F2.subscore  : expected item true trait score measured by the second statement at a point of the trait measured by the second statement
#  ..$ F2.F3.subscore  : expected item true trait score measured by the third statement at a point of the trait measured by the second statement
#  ..$ F3.F1.subscore  : expected item true trait score measured by the first statement at a point of the trait measured by the third statement
#  ..$ F3.F2.subscore  : expected item true trait score measured by the second statement at a point of the trait measured by the third statement
#  ..$ F3.F3.subscore  : expected item true trait score measured by the third statement at a point of the trait measured by the third statement
#
#4. data.frame  :  item true trait score at a theta vector with negative statement(s) changed to positive one(s)    
#  ..$ F1       : trait point measured by the first statement
#  ..$ F2       : trait point measured by the second statement
#  ..$ F3       : trait point measured by the third statement
#  ..$ weight   : weight (probability) of the theta vector (F1, F2, F3)
#  ..$ F1.subscore : expected item trait scores measured by the first statement at (F1, F2, F3)
#  ..$ F2.subscore : expected item trait scores measured by the second statement at (F1, F2, F3)
#  ..$ F3.subscore : expected item trait scores measured by the third statement at (F1, F2, F3)
#  ..$ F1.subscore.w : F1.subscore*weight
#  ..$ F2.subscore.w : F2.subscore*weight
#  ..$ F3.subscore.w : F3.subscore*weight

#5. data.frame  : expected (marginal) item true trait score at a theta point with negative statement(s) changed to positive one(s)
#  ..$ Theta  : a theta point
#  ..$ F1.F1.subscore  : expected item true trait score measured by the first statement at a point of the trait measured by the first statement
#  ..$ F1.F2.subscore  : expected item true trait score measured by the second statement at a point of the trait measured by the first statement
#  ..$ F1.F3.subscore  : expected item true trait score measured by the third statement at a point of the trait measured by the first statement
#  ..$ F2.F1.subscore  : expected item true trait score measured by the first statement at a point of the trait measured by the second statement
#  ..$ F2.F2.subscore  : expected item true trait score measured by the second statement at a point of the trait measured by the second statement
#  ..$ F2.F3.subscore  : expected item true trait score measured by the third statement at a point of the trait measured by the second statement
#  ..$ F3.F1.subscore  : expected item true trait score measured by the first statement at a point of the trait measured by the third statement
#  ..$ F3.F2.subscore  : expected item true trait score measured by the second statement at a point of the trait measured by the third statement
#  ..$ F3.F3.subscore  : expected item true trait score measured by the third statement at a point of the trait measured by the third statement
#
#6. matrix  : adjusted observed item trait score for each examinee with negative statement(s) changed to positive one(s) 
#  ..$ score  : observed item score
#  ..$ F1.subscore  : observed item trait score measured by the first statement
#  ..$ F2.subscore  : observed item trait score measured by the second statement
#  ..$ F3.subscore  : observed item trait score measured by the third statement
#  ..$ F1      : trait score estimate by the first statement
#  ..$ F2      : trait score estimate by the second statement
#  ..$ F3      : trait score estimate by the third statement
#
##########################################################################

ICC.Triplet.2PL <-
function (item, model.est, dimcor, par.est, theta.up, theta.low, 
    theta.interval = NULL, nquad = NULL, data, theta.score) 
{
    block.size <- 3
    nf <- nrow(dimcor)
    if (is.null(theta.interval)) {
        tmp <- round(seq(theta.low, theta.up, length.out = nquad), 
            2)
    }
    else {
        tmp <- round(seq(theta.low, theta.up, by = theta.interval), 
            2)
    }
    if (block.size == 3) {
        tmp1 <- as.matrix(expand.grid(tmp, tmp, tmp))
    }
    else if (block.size == 2) {
        tmp1 <- as.matrix(expand.grid(tmp, tmp))
    }
    tmp1.1 <- as.matrix(expand.grid(tmp, tmp))
    tmp1.dim <- par.est[item, paste0("dim", 1:block.size)]
    colnames(tmp1) <- paste0("F", 1:length(tmp1.dim))
    tmp1.d <- dmvnorm(tmp1, sigma = dimcor[tmp1.dim, tmp1.dim])
    tmp1.d <- tmp1.d/sum(tmp1.d)
    tmp2.d <- dnorm(tmp)
    tmp2.d <- cbind(tmp, tmp2.d/sum(tmp2.d))
    colnames(tmp2.d) <- c("Group", "Weight.s")
    theta <- matrix(0, nc = nf, nr = nrow(tmp1), byrow = T)
    theta[, par.est[item, paste0("dim", 1:block.size)]] <- tmp1
    F1.score <- c(2, 2, 1, 0, 1, 0)
    F2.score <- c(1, 0, 2, 2, 0, 1)
    F3.score <- c(0, 1, 0, 1, 2, 2)
    i.trace <- probtrace(extract.item(model.est, item), Theta = theta.score[, 
        1:nf])
    i.exp <- cbind(apply(t(t(i.trace) * F1.score), 1, sum), apply(t(t(i.trace) * 
        F2.score), 1, sum), apply(t(t(i.trace) * F3.score), 1, 
        sum))
    icc.trace <- probtrace(extract.item(model.est, item), Theta = theta)
    icc.exp <- cbind(apply(t(t(icc.trace) * F1.score), 1, sum), 
        apply(t(t(icc.trace) * F2.score), 1, sum), apply(t(t(icc.trace) * 
            F3.score), 1, sum))
    icc.exp.t <- cbind(tmp1, tmp1.d, icc.exp, icc.exp * tmp1.d)
    colnames(icc.exp.t) <- c(paste0("F", 1:length(tmp1.dim)), 
        "weight", paste0("F", 1:length(tmp1.dim), ".subscore"), 
        paste0("F", 1:length(tmp1.dim), ".subscore.w"))
    i.exp.t <- cbind(data, F1.score[data], F2.score[data], F3.score[data], 
        i.exp, theta.score[, tmp1.dim])
    colnames(i.exp.t) <- c("score", paste0("F", 1:length(tmp1.dim), 
        ".subscore"), paste0("F", 1:length(tmp1.dim), ".subscore.t"), 
        paste0("F", 1:length(tmp1.dim)))
    for (i in 1:block.size) {
        tmp5 <- aggregate(icc.exp.t[, paste0("F", 1:length(tmp1.dim), 
            ".subscore.w")], by = list(icc.exp.t[, paste0("F", 
            i)]), FUN = sum)
        colnames(tmp5) <- c("Group", paste0("F", i, ".F", 1:length(tmp1.dim), 
            ".subscore"))
        if (i == 1) 
            icc.E <- tmp5
        else icc.E <- merge(icc.E, tmp5, by = "Group", sort = F)
    }
    icc.E1 <- merge(tmp2.d, icc.E, sort = F)
    icc.E2 <- cbind(Theta = icc.E1[, 1], icc.E1[, 3:ncol(icc.E1)]/icc.E1[, 
        2])
    sign.tmp <- sign(par.est[item, paste0("a", 1:block.size)])
    sign.tmp[sign.tmp == 0] <- 1
    icc.exp.pa.t <- icc.exp.t
    icc.E2.pa <- icc.E2
    if (sum(sign.tmp) == -3) {
        icc.trace.pa <- icc.trace[, c(6, 4, 5, 2, 3, 1)]
        icc.exp.pa <- cbind(apply(t(t(icc.trace.pa) * F1.score), 
            1, sum), apply(t(t(icc.trace.pa) * F2.score), 1, 
            sum), apply(t(t(icc.trace.pa) * F3.score), 1, sum))
        icc.exp.pa.t <- cbind(tmp1, tmp1.d, icc.exp.pa, icc.exp.pa * 
            tmp1.d)
        colnames(icc.exp.pa.t) <- c(paste0("F", 1:length(tmp1.dim)), 
            "weight", paste0("F", 1:length(tmp1.dim), ".subscore"), 
            paste0("F", 1:length(tmp1.dim), ".subscore.w"))
        for (i in 1:block.size) {
            tmp6 <- aggregate(icc.exp.pa.t[, paste0("F", 1:length(tmp1.dim), 
                ".subscore.w")], by = list(icc.exp.pa.t[, paste0("F", 
                i)]), FUN = sum)
            colnames(tmp6) <- c("Group", paste0("F", i, ".F", 
                1:length(tmp1.dim), ".subscore"))
            if (i == 1) 
                icc.E.pa <- tmp6
            else icc.E.pa <- merge(icc.E.pa, tmp6, by = "Group", 
                sort = F)
        }
        icc.E1.pa <- merge(tmp2.d, icc.E.pa, sort = F)
        icc.E2.pa <- cbind(Theta = icc.E1.pa[, 1], icc.E1.pa[, 
            3:ncol(icc.E1.pa)]/icc.E1.pa[, 2])
    }
    else if (sum(sign.tmp) %in% c(-1, 1)) {
        sign.n <- which(sign.tmp == -1)
        a <- abs(par.est[item, paste0("a", 1:block.size)])
        b <- par.est[item, paste0("b", 1:block.size)]
        b[sign.n] <- -b[sign.n]
        icc.trace.pa <- P.Triplet.2PL.s(c(a, b)[c(1, 4, 2, 5, 
            3, 6)], tmp1[, 1], tmp1[, 2], tmp1[, 3])
        icc.exp.pa <- cbind(apply(t(t(icc.trace.pa) * F1.score), 
            1, sum), apply(t(t(icc.trace.pa) * F2.score), 1, 
            sum), apply(t(t(icc.trace.pa) * F3.score), 1, sum))
        icc.exp.pa.t <- cbind(tmp1, tmp1.d, icc.exp.pa, icc.exp.pa * 
            tmp1.d)
        colnames(icc.exp.pa.t) <- c(paste0("F", 1:length(tmp1.dim)), 
            "weight", paste0("F", 1:length(tmp1.dim), ".subscore"), 
            paste0("F", 1:length(tmp1.dim), ".subscore.w"))
        for (i in 1:block.size) {
            tmp6 <- aggregate(icc.exp.pa.t[, paste0("F", 1:length(tmp1.dim), 
                ".subscore.w")], by = list(icc.exp.pa.t[, paste0("F", 
                i)]), FUN = sum)
            colnames(tmp6) <- c("Group", paste0("F", i, ".F", 
                1:length(tmp1.dim), ".subscore"))
            if (i == 1) 
                icc.E.pa <- tmp6
            else icc.E.pa <- merge(icc.E.pa, tmp6, by = "Group", 
                sort = F)
        }
        icc.E1.pa <- merge(tmp2.d, icc.E.pa, sort = F)
        icc.E2.pa <- cbind(Theta = icc.E1.pa[, 1], icc.E1.pa[, 
            3:ncol(icc.E1.pa)]/icc.E1.pa[, 2])
    }
    tmp <- i.exp.t[, c("score", paste0("F", 1:length(tmp1.dim)))]
    i.exp.pa <- i.exp.t[, c("score", paste0("F", 1:length(tmp1.dim), 
        ".subscore"), paste0("F", 1:length(tmp1.dim)))]
    if (sum(sign.tmp) == -3) {
        tmp1 <- tmp[, "score"]
        data1 <- ifelse(tmp1 == 1, 6, ifelse(tmp1 == 2, 4, ifelse(tmp1 == 
            3, 5, ifelse(tmp1 == 4, 2, ifelse(tmp1 == 5, 3, ifelse(tmp1 == 
            6, 1, NA))))))
        i.exp.pa <- cbind(data1, F1.score[data1], F2.score[data1], 
            F3.score[data1], tmp[, paste0("F", 1:length(tmp1.dim))])
    }
    else if (sum(sign.tmp) %in% c(-1, 1)) {
        sign.n <- which(sign.tmp == -1)
        a <- abs(par.est[item, paste0("a", 1:block.size)])
        b <- par.est[item, paste0("b", 1:block.size)]
        b[sign.n] <- -b[sign.n]
        i.trace.pa <- P.Triplet.2PL.s(c(a, b)[c(1, 4, 2, 5, 3, 
            6)], tmp[, "F1"], tmp[, "F2"], tmp[, "F3"])
        i.exp.pa <- cbind(apply(t(t(i.trace.pa) * F1.score), 
            1, sum), apply(t(t(i.trace.pa) * F2.score), 1, sum), 
            apply(t(t(i.trace.pa) * F3.score), 1, sum))
        i.exp.pa <- cbind(NA, i.exp.pa, tmp[, paste0("F", 1:length(tmp1.dim))])
    }
    colnames(i.exp.pa) <- c("score", paste0("F", 1:length(tmp1.dim), 
        ".subscore"), paste0("F", 1:length(tmp1.dim)))
    return(list(icc.exp.t, i.exp.t, icc.E2, icc.exp.pa.t, icc.E2.pa, 
        i.exp.pa))
}

#Create ICC data for Item 2 in a test

Triplet.2PL.5.fix.ICC.2 <- ICC.Triplet.2PL(item=2, 
	model.est=Triplet.2PL.5.fix, 
	dimcor=Trait.cor.3, 
	par.est=Triplet.2PL.5.par.fix, 
	theta.up=3, 
	theta.low=-3, 
	theta.interval=NULL, 
	nquad=31,  
	data= sim.data.3[[5]][, 2], 
	theta.score=Triplet.2PL.5.MAP.fix)

#Create ICC data for all items in a test

Triplet.2PL.5.fix.ICC <- vector('list', 10)
for ( i in 1:10){
Triplet.2PL.5.fix.ICC[[i]] <- ICC.Triplet.2PL(item=i, 
	model.est=Triplet.2PL.5.fix, 
	dimcor=Trait.cor.3, 
	par.est=Triplet.2PL.5.par.fix, 
	theta.up=3, 
	theta.low=-3, 
	theta.interval=NULL, 
	nquad=31,  
	data= sim.data.3[[5]][, i], 
	theta.score=Triplet.2PL.5.MAP.fix)
}

#########################################################################
#Function to draw Triplet-2PL ICC plots using the ICC data output from the ICC.Triplet.2PL function 
#
#Description:
#Draw one or multiple Triplet-2PL ICC plots using the ICC data output from the ICC.Triplet.2PL function
#
#Usage:
#ICC.plot.Triplet.2PL(ICC.data, plot.item, par.est, PageNo = T, fit) 
#
#Arguments:
#ICC.data -- a list with Nitem elements being the ICC data output of the ICC.Triplet.2PL function for all items in a test 
#plot.item -- a numeric vector containing item number(s) for which the ICC(s) is drawn 
#par.est -- a Nitem-row matrix containing the following colnames: 
#           discrimination and intercept parameters for Statements 1-3 in each item: "a1"   "b1"   "a2"   "b2"  "a3"   "b3" 
#		trait number(1:Ntrait) measured by each statement in each item: "dim1" "dim2"  "dim3"
#
#PageNo -- logic; T = print page number at the bottom center, F = no page number
#fit -- the output of Zh item fit statistics from the itemfit function in the mirt package
#
#Return:
#Output the ICC plot(s) for the input item(s)
##########################################################################
ICC.plot.Triplet.2PL <-
function (ICC.data, plot.item, par.est, PageNo = T, fit) 
{   
    block.size <- 3
    for (i in plot.item) {
        i.dim <- par.est[i, paste0("dim", 1:block.size)]
        par.item <- matrix(par.est[i, paste0(c("a", "b"), rep(1:block.size, 
            each = 2))], ncol = 2, byrow = T)
        state.text <- c("A", "B", "C")[1:block.size]
        i.dat <- ICC.data[[i]][[2]][, paste0("F", 1:3, ".subscore")]
        i.dat <- i.dat[apply(!is.na(i.dat), 1, all), ]
        i.dat.m <- round(apply(i.dat, 2, mean), 2)
        group.text <- paste0("F", i.dim, "-Statement ", state.text, 
            " (a= ", par.item[, 1], ", b= ", par.item[, 2], ") Mean_obs = ", 
            i.dat.m)
        group.text <- factor(group.text, levels = group.text)
        key.text <- paste0("Item", i, " N=", nrow(i.dat), " Zh=", 
            round(fit[i, 2], 2), " P(Zh)=", round(pnorm(fit[i, 
                2]), 2))
        i.exp.tmp <- ICC.data[[i]][[2]][apply(!is.na(i.dat), 
            1, all), c(paste0("F", 1:3, ".subscore"), paste0("F", 
            1:block.size))]
        tmp1 <- data.frame(group = group.text[1], i.exp.tmp[, 
            c("F1.subscore", "F1")])
        tmp2 <- data.frame(group = group.text[2], i.exp.tmp[, 
            c("F2.subscore", "F2")])
        tmp3 <- data.frame(group = group.text[3], i.exp.tmp[, 
            c("F3.subscore", "F3")])
        colnames(tmp1) <- colnames(tmp2) <- colnames(tmp3) <- c("group", 
            "Item", "Theta")
        i.exp.tmp1 <- rbind(tmp1, tmp2, tmp3)
        f.tmp <- as.formula(paste0("Item", "~Theta|group"))
        col.tmp <- sum(par.item[, 1] >= 0) + 1
        bg.col <- c(2, 3, 4, 1)[col.tmp]
        plot1 <- xyplot(f.tmp, data = i.exp.tmp1, panel = function(x, 
            y, ...) {
            panel.grid(h = -1, v = -1, col.line = grey(0.9))
            panel.smoother(x, y, ..., family = "g", n = 100, 
                span = 0.9, lty = 2)
        }, ylab = list(label = "Expected Item Trait Score", col = 1, 
            cex = 1), strip = strip.custom(strip.names = c(F, 
            TRUE), par.strip.text = list(cex = 0.7), bg = trellis.par.get("strip.background")$col[[bg.col]]), 
            key = list(text = list(c("Expected", "Observed(with 95% CI)"), 
                cex = 0.7, col = c(2, 4)), lines = list(lty = 1:2, 
                col = c(2, 4), cex = 0.7), space = "top", border = F, 
                between = 1, adj = 0, column = 2), layout = c(1, 
                3), main = list(label = key.text, cex = 1, fontface = 1), 
            as.table = T, sub = if (PageNo) 
                list(label = i, cex = 0.7)
            else "")
        tmp1 <- data.frame(group = group.text[1], ICC.data[[i]][[3]][, 
            c("F1.F1.subscore", "Theta")])
        tmp2 <- data.frame(group = group.text[2], ICC.data[[i]][[3]][, 
            c("F2.F2.subscore", "Theta")])
        tmp3 <- data.frame(group = group.text[3], ICC.data[[i]][[3]][, 
            c("F3.F3.subscore", "Theta")])
        colnames(tmp1) <- colnames(tmp2) <- colnames(tmp3) <- c("group", 
            "Item", "Theta")
        icc.E2.tmp <- rbind(tmp1, tmp2, tmp3)
        plot2 <- xyplot(f.tmp, icc.E2.tmp, type = "l", col = 2, 
            layout = c(1, 3), as.table = T)
        print(plot1 + as.layer(plot2))
    }
}
#draw the second item's ICC plot to screen
ICC.plot.Triplet.2PL(ICC.data=Triplet.2PL.5.fix.ICC, plot.item=2, par.est=Triplet.2PL.5.par.fix,PageNo=F, fit=Triplet.2PL.5.ZH.fix)

#draw the ICC plots of all items in test and save to a PDF file

pdf("Triplet_2PL_5_fix_ICC_plot.pdf")

ICC.plot.Triplet.2PL(ICC.data=Triplet.2PL.5.fix.ICC, plot.item=1:10, par.est=Triplet.2PL.5.par.fix,PageNo=T, fit=Triplet.2PL.5.ZH.fix)

dev.off()

#########################################################################
#Function to create TCC data for the Triplet-2PL
#
#Description:
#Create TCC data for the Triplet-2PL using the ICC data output from the ICC.Triplet.2PL function
#
#Usage:
#TCC.Triplet.2PL(ICC.data, par.est, theta.score)
#
#Arguments:
#ICC.data -- a list with Nitem elements being the ICC data output of the ICC.Triplet.2PL function for all items in a test  
#par.est -- a Nitem-row matrix containing the following colnames: 
#           discrimination and intercept parameters for Statements 1-3 in each item: "a1"   "b1"   "a2"   "b2"  "a3"   "b3" 
#		trait number(1:Ntrait) measured by each statement in each item: "dim1" "dim2"  "dim3"
#theta.score -- a matrix (Nexaminee X Ntrait) containing trait score(s) for each examinee
#
#Return: 
#A list with 4 elements. 
#
#1. matrix  :  expected trait test true score at a theta value with negative statement(s) changed to positive one(s)   
#  ..$ Theta    : trait score point defined by theta.up, theta.low, theta.interval, or nquad in the ICC.MUPP.2PL function  
#  ..$ F1.subscore : Trait 1's test scores
#  ..$ F2.subscore : Trait 2's test scores
#  ... 
#  ..$ F&Ntrait.subscore : last Trait's test scores
#  
#2. matrix  :  observed trait test scores for each examinee without changing the direction of negative statement(s)   
#  ..$ F1    : Trait 1's score estimate 
#  ..$ F2    : Trait 2's score estimate 
#  ... 
#  ..$ F&Ntrait : last Trait's score estimate
#  ..$ F1.subscore : Trait 1's observed test scores
#  ..$ F2.subscore : Trait 2's observed test scores
#  ... 
#  ..$ F&Ntrait.subscore : last Trait's observed test scores
#3. matrix  :  expected trait test true score at a theta value without changing the direction of negative statement(s)   
#  ..$ Theta    : trait score point defined by theta.up, theta.low, theta.interval, or nquad in the ICC.MUPP.2PL function  
#  ..$ F1.subscore : Trait 1's test scores
#  ..$ F2.subscore : Trait 2's test scores
#  ... 
#  ..$ F&Ntrait.subscore : last Trait's test scores
#  
#4. matrix  :  observed trait test scores for each examinee with negative statement(s) changed to positive one(s)   
#  ..$ F1    : Trait 1's score estimate 
#  ..$ F2    : Trait 2's score estimate 
#  ... 
#  ..$ F&Ntrait : last Trait's score estimate
#  ..$ F1.subscore : Trait 1's observed test scores
#  ..$ F2.subscore : Trait 2's observed test scores
#  ... 
#  ..$ F&Ntrait.subscore : last Trait's observed test scores
##########################################################################

TCC.Triplet.2PL <-
function (ICC.data, par.est, theta.score) 
{
    nf <- max(par.est[, c("dim1", "dim2", "dim3")])
    block.size <- 3
    subscore.t.adj <- matrix(0, nr = nrow(ICC.data[[1]][[5]]), 
        nc = nf)
    subscore.t <- matrix(0, nr = nrow(ICC.data[[1]][[3]]), nc = nf)
    subscore.obs <- subscore.obs.adj <- matrix(0, nr = nrow(ICC.data[[1]][[6]]), 
        nc = nf)
    for (i in 1:length(ICC.data)) {
        tmp1.dim <- par.est[i, paste0("dim", 1:block.size)]
        tmp2 <- as.matrix(ICC.data[[i]][[5]][, paste0("F", 1:3, 
            ".F", 1:3, ".subscore")])
        tmp5 <- as.matrix(ICC.data[[i]][[6]][, paste0("F", 1:3, 
            ".subscore")])
        tmp3 <- as.matrix(ICC.data[[i]][[2]][, paste0("F", 1:3, 
            ".subscore")])
        tmp4 <- as.matrix(ICC.data[[i]][[3]][, paste0("F", 1:3, 
            ".F", 1:3, ".subscore")])
        subscore.t.adj[, tmp1.dim] <- subscore.t.adj[, tmp1.dim] + 
            tmp2
        subscore.obs[, tmp1.dim] <- subscore.obs[, tmp1.dim] + 
            tmp3
        subscore.t[, tmp1.dim] <- subscore.t[, tmp1.dim] + tmp4
        subscore.obs.adj[, tmp1.dim] <- subscore.obs.adj[, tmp1.dim] + 
            tmp5
    }
    subscore.t.adj <- cbind(ICC.data[[1]][[5]][, 1], subscore.t.adj)
    subscore.obs <- cbind(theta.score[, 1:nf], subscore.obs)
    subscore.obs.adj <- cbind(theta.score[, 1:nf], subscore.obs.adj)
    subscore.t <- cbind(ICC.data[[1]][[3]][, 1], subscore.t)
    colnames(subscore.t) <- colnames(subscore.t.adj) <- c("Theta", 
        paste0("F", 1:nf, ".subscore"))
    colnames(subscore.obs) <- colnames(subscore.obs.adj) <- c(paste0("F", 
        1:nf), paste0("F", 1:nf, ".subscore"))
    return(list(subscore.t.adj, subscore.obs, subscore.t, subscore.obs.adj))
}

#create TCC data
Triplet.2PL.5.fix.TCC <- TCC.Triplet.2PL(ICC.data=Triplet.2PL.5.fix.ICC, par.est=Triplet.2PL.5.par.fix, theta.score=Triplet.2PL.5.MAP.fix)

#draw Trait 2's TCC plot to screen
TCC.plot(TCC.data=Triplet.2PL.5.fix.TCC, plot.Factors=2, par.est=Triplet.2PL.5.par.fix, trait=trait.label, PageNo = F)
TCC.plot(TCC.data=Triplet.2PL.5.fix.TCC, plot.Factors=2, par.est=Triplet.2PL.5.par.fix, trait=trait.label, PageNo = F, adj.a.sign=T)


#draw the TCC plots of all traits in a test and save to a PDF file

pdf("Triplet_2PL_5_fix_TCC_plot.pdf")

TCC.plot(TCC.data=Triplet.2PL.5.fix.TCC, plot.Factors=1:5, par.est=Triplet.2PL.5.par.fix, trait=trait.label, PageNo = T)

dev.off()

#draw the TCC plots of all traits in a test with negative statement(s) changed to positive one(s),
#and save to a PDF file

pdf("Triplet_2PL_5_fix_TCC_plot_adj.pdf")

TCC.plot(TCC.data=Triplet.2PL.5.fix.TCC, plot.Factors=1:5, par.est=Triplet.2PL.5.par.fix, trait=trait.label, PageNo = T, adj.a.sign=T)

dev.off()


#save workspace
save.image(".rdata")

