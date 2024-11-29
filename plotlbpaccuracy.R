# Reproduce "Accuracy of loopy BP" figures in main text (Fig 6) and supplementary
# material (Fig S3, S4), which compare clique tree and cluster graph estimates
# for 2 networks

library(scales)
library(SiPhyNetwork)
################################################################################
## Figure 6: lbp accuracy for regularization algorithm R1 vs R2
muller_net <- read.net("real_networks/muller_2022.phy")
getNetworkLevel(muller_net) # 358
# join-graph structuring loopy BP estimates using regularization algorithms R1 and R2
df_muller_R1 <- read.table("lbp_accuracy/muller2022_joingraph_R1.csv", header=T, sep=",")
df_muller_R2 <- read.table("lbp_accuracy/muller2022_joingraph_R2.csv", header=T, sep=",")
posmean_muller <- function(){
  plot(x=1:200, y=df_muller_R2$posmean[1:200], type="l", ann=F, cex.axis=1.3)
  lines(x=1:200, y=df_muller_R1$posmean[1:200], type="l", lty="dotted")
  abline(h=-15.696383011725299, col="red")}
posvar_muller <- function(){
  plot(x=1:50, y=df_muller_R2$posvar[1:50], type="l", ann=F,
       ylim=c(0,max(df_muller_R1$posvar, df_muller_R2$posvar)), cex.axis=1.3)
  lines(x=1:50, y=df_muller_R1$posvar[1:50], type="l", lty="dotted")
  abline(h=1425.651004065571, col="red")}
fenergy_muller <- function(){
  plot(x=1:50, y=df_muller_R2$fenergy[1:50], type="l", ann=F, cex.axis=1.3)
  lines(x=5:50, y=df_muller_R1$fenergy[5:50], type="l", lty="dotted")
  abline(h=-85.04744054045547, col="red")}

lipson_net <- read.net("real_networks/lipson_2020b.phy")
getNetworkLevel(lipson_net) # 12
# join-graph structuring loopy BP estimates using regularization algorithms R1 and R2
df_lipson_R1 <- read.table("lbp_accuracy/lipson2020b_joingraph_R1.csv", header=T, sep=",")
df_lipson_R2 <- read.table("lbp_accuracy/lipson2020b_joingraph_R2.csv", header=T, sep=",")
posmean_lipson <- function(){
  plot(x=1:20, y=df_lipson_R2$posmean, type="l", ann=F,
       ylim=range(c(df_lipson_R1$posmean, df_lipson_R2$posmean)), cex.axis=1.3)
  lines(x=1:20, y=df_lipson_R1$posmean, type="l", lty="dotted")
  abline(h=1.7010205811114112, col="red")
  legend(x=5, y=1, legend=c("clique tree", "join-graph str, R1",
                             "join-graph str, R2"),
         cex=1.6, x.intersp=.5, bty="n", lty=c("solid", "dotted","solid"),
         col=c("red", "black", "black"), seg.len=1, lwd=1.5)}
posvar_lipson <- function(){
  plot(x=1:20, y=df_lipson_R2$posvar, type="l", ann=F,
       ylim=c(0,max(df_lipson_R1$posvar, df_lipson_R2$posvar)), cex.axis=1.3)
  lines(x=1:20, y=df_lipson_R1$posvar, type="l", lty="dotted")
  abline(h=31.032000074578907, col="red")}
fenergy_lipson <- function(){
  plot(x=1:20, y=df_lipson_R2$fenergy, type="l", ann=F, cex.axis=1.3)
  lines(x=1:20, y=df_lipson_R1$fenergy, type="l", lty="dotted")
  abline(h=-35.89254355788221, col="red")}

combined_plot <- function(){
  par(mar=c(2.4, 2.2, 2.4, 2.2), oma=c(2,1,1,1), las=1)
  layout.matrix <- matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T)
  layout(mat=layout.matrix)
  
  posmean_lipson(); posvar_lipson(); fenergy_lipson()
  posmean_muller(); posvar_muller(); fenergy_muller()
  
  mtext(text="Number of iterations", side=1, outer=T, at=0.5, line=1, cex=1.1) # x-axis label
  mtext(text=expression(E * group("(",X[rho] ~ group("|",data,""),")")), side=3,
        outer=T, line=-17.3, at=0.15, cex=1.1) # y-axis label top-left
  mtext(text=expression(E * group("(",X[rho] ~ group("|",data,""),")")), side=3,
        outer=T, line=-38.5, at=0.15, cex=1.1) # y-axis label bottom-left
  mtext(text=expression(Var * group("(",X[rho] ~ group("|",data,""),")")), side=3,
        outer=T, line=-17.3, at=0.5, cex=1.1) # y-axis label top-mid
  mtext(text=expression(Var * group("(",X[rho] ~ group("|",data,""),")")), side=3,
        outer=T, line=-38.5, at=0.5, cex=1.1) # y-axis label bottom-mid
  mtext(text="Factored energy", side=3, outer=T, line=-16.8, at=0.83, cex=1.1) # y-axis label top-right
  mtext(text="Factored energy", side=3, outer=T, line=-38, at=0.83, cex=1.1) # y-axis label bottom-right
  mtext(text=expression("Lipson et al. (2020b): " *
                          list(italic(n)==12, "\u2113"==12, italic(h)==12, k==7, k^symbol("*")==3)),
        side=3, outer=T, line=-2, at=0, adj=0, cex=1.1) # top row label
  mtext(text=expression("Müller et al. (2022): " *
                          list(italic(n)==40, "\u2113"==358, italic(h)==361, k==54, k^symbol("*")==10)),
        side=3, outer=T, line=-23.2, at=0, adj=0, cex=1.1) # bottom row label
}

cairo_pdf(file="figures/lbpaccuracy.pdf", width=9, height=6)
combined_plot()
dev.off()

################################################################################
## Figure S3: lbp accuracy for join-graph structuring vs factor graph
# applies modulus transformation with transformation exponent set to 0
# (https://www.jstor.org/stable/2986305)
logmodulus_trans <- function(x) {
  modulus_trans(0)$transform(x)
}

# factor graph loopy BP estimates using regularization algorithm R2
df_muller_fg_R2 <- read.table("lbp_accuracy/muller2022_factorgraph_R2.csv",
                              header=T, sep=",")
# range(df_muller_R1$posmean[1:200],df_muller_R2$posmean[1:200]) # -77.64321  24.45600
posmean_muller_fg <- function(){
  plot(x=1:200, y=df_muller_fg_R2$posmean[1:200], type="l", ann=F,
       ylim=c(-77.64321,24.45600), col="purple", cex.axis=1.2)
  lines(x=1:200, y=df_muller_R2$posmean[1:200], type="l")
  abline(h=-15.696383011725299, col="red")
}
posmean_muller_fg_logmodulus <- function(){
  plot(x=1:200,
       y=logmodulus_trans(df_muller_fg_R2$posmean[1:200]+15.696383011725299),
       type="l", col="purple", ann=F, yaxt="n", cex.axis=1.2)
  
  axis(side=2, at=logmodulus_trans(c(-10^12,-10^8,-10^4,0,10^4,10^8,10^12)),
       cex.axis=1.2,
       labels=c(expression("-" * 10^12), expression("-" * 10^8),
                expression("-" * 10^4), 0, expression(10^4), expression(10^8),
                expression(10^12)), hadj=0.9)
  lines(x=1:200,
        y=logmodulus_trans(df_muller_R2$posmean[1:200]+15.696383011725299),
        type="l")
  abline(h=0, col="red")
}
# max(df_muller_R1$posvar[1:50],df_muller_R2$posvar[1:50]) # 2841.798
posvar_muller_fg <- function(){
  plot(x=1:50, y=df_muller_fg_R2$posvar[1:50], type="l", ann=F, ylim=c(0,2841.798),
       col="purple", cex.axis=1.2)
  lines(x=1:50, y=df_muller_R2$posvar[1:50], type="l")
  abline(h=1425.651004065571, col="red")
}
posvar_muller_fg_logmodulus <- function(){
  plot(x=1:50, y=logmodulus_trans(df_muller_fg_R2$posvar[1:50]-1425.651004065571),
       type="l", col="purple", ann=F, yaxt="n", cex.axis=1.2)
  axis(side=2, at=logmodulus_trans(c(-10^2,-10^1,0,10^1,10^2,10^4)), cex.axis=1.2,
       labels=c(expression("-" * 10^2), -10, 0, 10, expression(10^2), expression(10^4)), hadj=0.9)
  lines(x=1:50, y=logmodulus_trans(df_muller_R2$posvar[1:50]-1425.651004065571),
        type="l")
  abline(h=0, col="red")
}
# range(df_muller_R1$fenergy[5:50],df_muller_R2$fenergy[1:50]) # -589.68140  -82.07795
fenergy_muller_fg <- function(){
  plot(x=1:50, y=df_muller_fg_R2$fenergy[1:50], type="l", ann=F, col="purple",
       ylim=c(-589.68140,-82.07795), cex.axis=1.2)
  lines(x=1:50, y=df_muller_R2$fenergy[1:50], type="l")
  abline(h=-85.04744054045547, col="red")}
fenergy_muller_fg_logmodulus <- function(){
  plot(x=1:50, y=logmodulus_trans(df_muller_fg_R2$fenergy[1:50]+85.04744054045547),
       type="l", ann=F, col="purple", ylim=c(-13,2.5), yaxt="n", cex.axis=1.2)
  axis(side=2, at=logmodulus_trans(c(-10^4,-10^2,-10^1, -2, 0, 2, 10)), cex.axis=1.2,
       labels=c(expression("-" * 10^4), expression("-" * 10^2), -10, -2, 0, 2, 10), hadj=0.9)
  lines(x=1:50, y=logmodulus_trans(df_muller_R2$fenergy[1:50]+85.04744054045547),
        type="l")
  abline(h=0, col="red")}

# factor graph loopy BP estimates using regularization algorithm R2
df_lipson_fg_R2 <- read.table("lbp_accuracy/lipson2020b_factorgraph_R2.csv",
                           header=T, sep=",")
posmean_lipson_fg <- function(){
  plot(x=1:20, y=df_lipson_fg_R2$posmean, type="l", ann=F,
       ylim=range(df_lipson_R1$posmean,df_lipson_R2$posmean), col="purple", cex.axis=1.2)
  lines(x=1:20, y=df_lipson_R2$posmean, type="l")
  abline(h=1.7010205811114112, col="red")
  legend(x=8.5, y=1,
         legend=c("clique tree", "join-graph str", "factor graph"),
         cex=1.5, bty="n", col=c("red", "black","purple"), x.intersp=.5,
         seg.len=1, lwd=1.5)}
posvar_lipson_fg <- function(){
  plot(x=1:20, y=df_lipson_fg_R2$posvar, type="l", ann=F,
       ylim=c(0,max(df_lipson_R1$posvar,df_lipson_R2$posvar)), col="purple", cex.axis=1.2)
  lines(x=1:20, y=df_lipson_R2$posvar, type="l")
  abline(h=31.032000074578907, col="red")}
fenergy_lipson_fg <- function(){
  plot(x=1:20, y=df_lipson_fg_R2$fenergy, type="l", ann=F,
       ylim=range(df_lipson_R1$fenergy,df_lipson_R2$fenergy), col="purple", cex.axis=1.2)
  lines(x=1:20, y=df_lipson_R2$fenergy, type="l")
  abline(h=-35.89254355788221, col="red")}

combined_plot_fg <- function(){
  par(mar=c(2.4, 2.2, 2.4, 2.8), oma=c(2,2.5,1,1), las=1)
  layout.matrix <- matrix(c(1,2,3,4,5,6,7,8,9),nrow=3,ncol=3,byrow=T)
  layout(mat=layout.matrix)
  
  posmean_lipson_fg(); posvar_lipson_fg(); fenergy_lipson_fg()
  posmean_muller_fg(); posvar_muller_fg(); fenergy_muller_fg()
  posmean_muller_fg_logmodulus(); posvar_muller_fg_logmodulus(); fenergy_muller_fg_logmodulus()
  
  mtext(text="Number of iterations", side=1, outer=T, at=0.5, line=1) # x-axis label
  mtext(text=expression(E * group("(",X[rho] ~ group("|",data,""),")")),
        side=3, outer=T, line=-17.3, at=0.15) # y-axis label top-left
  mtext(text=expression(E * group("(",X[rho] ~ group("|",data,""),")")), side=3,
        outer=T, line=-38.5, at=0.15) # y-axis label mid-left
  mtext(text=expression(Delta == E * group("(",X[rho] ~ group("|",data,""),")") - 15.7),
        side=3, outer=T, line=-62.7, at=0.15) # y-axis label bottom-left
  mtext(text=expression(Var * group("(",X[rho] ~ group("|",data,""),")")), side=3,
        outer=T, line=-17.3, at=0.5) # y-axis label top-mid
  mtext(text=expression(Var * group("(",X[rho] ~ group("|",data,""),")")), side=3,
        outer=T, line=-38.5, at=0.5) # y-axis label mid-mid
  mtext(text=expression(Delta == Var * group("(",X[rho] ~ group("|",data,""),")") - 1425.7),
        side=3, outer=T, line=-62.7, at=0.5) # y-axis label bottom-mid
  mtext(text="Factored energy", side=3,
        outer=T, line=-16.8, at=0.83) # y-axis label top-right
  mtext(text="Factored energy", side=3, outer=T, line=-38, at=0.83) # y-axis label mid-right
  mtext(text=expression(Delta == "Factored energy" - (-85)), side=3, outer=T,
        line=-62.7, at=0.83) # y-axis label bottom-right
  mtext(text=expression("Lipson et al. (2020b): " * italic(n)==12),
        side=3, outer=T, line=-2, at=0, adj=0) # row label top
  mtext(text=expression("Müller et al. (2022): " * italic(n)==40),
        side=3, outer=T, line=-23.7, at=0, adj=0) # row label mid
  mtext(text=expression("Müller et al. (2022): " * italic(n)==40 *
                          ", difference from clique tree estimate shown on log-modulus transformation scale"),
        side=3, outer=T, line=-45.5, at=0, adj=0) # row label mid
}

cairo_pdf(file="figures/lbpaccuracy_fg.pdf", width=9.5, height=9)
combined_plot_fg()
dev.off()

################################################################################
## Figure S4: cluster sizes for clique tree and join-graph structuring
df_clustersizes <- read.table("lbp_accuracy/lbpaccuracy_clustersizes.csv", header=T, sep=",")
clustersize_boxplots <- function(){
  col1 <- rgb(0,0,255,alpha=120,max=255); col2 <- rgb(238,130,238,alpha=120,max=255)
  boxplot(clustersize ~ grp, data=df_clustersizes, range=0, log="x",
          names=c(11, 8.3, 5.8, 4.9), horizontal=T, las=1, ylab="", xlab="",
          at=c(11, 8.3, 5.8, 4.9), boxfill=c(col1, col1, col2, col2), cex.lab=1.5)
  points(x=c(mean(df_clustersizes$clustersize[df_clustersizes$grp == 1]),
             mean(df_clustersizes$clustersize[df_clustersizes$grp == 2]),
             mean(df_clustersizes$clustersize[df_clustersizes$grp == 3]),
             mean(df_clustersizes$clustersize[df_clustersizes$grp == 4])),
         y=c(11, 8.3, 5.8, 4.9))
  title(ylab=expression("Mean time ("*mu*"s) per message"),
        xlab="Cluster size (log-scale)", cex.lab=1.5, line=2.5)
  legend(x=8, y=6.4,
         legend=c(expression("Müller (" * italic(n)==40 * ")"),
                  expression("Lipson (" * italic(n)==12 * ")")),
         x.intersp=.5, cex=1.5, bty="n",
         fill=c(rgb(0,0,255,alpha=120,max=255),rgb(238,130,238,alpha=120,max=255)))
  text(x=1.5,y=11,labels=expression(italic(U)))
  text(x=13,y=8.3,labels=expression(paste(italic(U),"*")))
  text(x=1.5,y=5.8,labels=expression(italic(U)))
  text(x=4,y=4.9,labels=expression(paste(italic(U),"*")))
}

pdf(file="figures/lbpaccuracy_clustersizeboxplots.pdf", width=6, height=6)
clustersize_boxplots()
dev.off()