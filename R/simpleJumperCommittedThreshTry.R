cols <- c(committed    = "#BD8E17FF",
          horizonFree = "#1F78B4E6",
          universal    = "#E31A1CE6")
labs <- c(latex2exp::TeX(r"($b_{N,\alpha}=11.20)"),
          latex2exp::TeX(r"($b_{\alpha}=18.69)"),
          latex2exp::TeX(r"($1/\alpha=20)"))

# N=200 -----
kip <- uniroot(function(x){
  processNullSamplePaths(eValueMatrixNull[, 1:200], b=x)$pB-0.05
}, c(5, 15))


chim <- processNullSamplePaths(eValueMatrixNull[, 1:200], b=11.19918)
chim$b*chim$pB+chim$vB

chimS <- processNullSamplePaths(eValueMatrixNull[, 1:200], b=20)
chimS$b*chimS$pB+chimS$vB

bunny <- uniroot(function(x){
  iets <- processNullSamplePaths(eValueMatrixNull[, 1:200], b=x)
  iets$pB+iets$rB/iets$b-0.05
}, c(12, 20))

chimH <- processNullSamplePaths(eValueMatrixNull[, 1:200], b=18.68956)
chimH$b*chimH$pB+chimH$vB
chimH$pB+chimH$rB/chimH$b


# Paths -----
cPaths <- processTimeSamplePaths(eValueMatrixNull, b=11.1992)
hfPaths <- processTimeSamplePaths(eValueMatrixNull, b=18.68956)
uPaths <- processTimeSamplePaths(eValueMatrixNull, b=20)
round(11.19918, 5)
round(uPaths$pB[200], 5)
round(uPaths$rB[2000], 5)
round(uPaths$vB[2000], 5)

uPaths$pB[200]
uPaths$pB[200]+uPaths$vB[200]/20

18.68956
# Fpr plot ----
pdf("conformalFpr.pdf", width = 14, height = 7)
par(mai = c(1.2, 1.5, 0.2, 1.7), cex.axis = 4/3, las = 1)

plot(NA, xlim = c(10, 2000), ylim = c(0, 0.08), log = "x",
     axes = FALSE, xlab = "", ylab = "")
axis(1, at = c(10, 20, 50, 100, 200, 500, 1000, 2000), cex.axis=cexFactor*myCex)
axis(2, at = seq(0, 0.08, 0.01), labels = 0:8, cex.axis=cexFactor*myCex)
mtext("Time horizon N", side = 1, line = 3.8, cex = cexFactor*myCex)

mtext(latex2exp::TeX(r"($P(T_{b} \leq N)$ in %)"),
      side = 2, line = 4, las=0,
      cex = cexFactor*myCex,
      adj=0.5, padj=-0.4)


lines(c(1, 2000), c(alpha, alpha),
      lwd=myLwd, col="lightgrey", lty=2)

lines(c(N, N), c(0, 0.08),
      lwd=myLwd, col="grey", lty=2)


# draw universal first so committed ends up on top, as in the original
for (k in 1:3) {
  if (k==1) {
    lines(10:200, cPaths$pB[10:200], col=cols[k], lwd=myLwd)
    lines(201:2000, cPaths$pB[201:2000], col="darkorchid3",
          lwd=myLwd, lty=2)
  } else {
    yLine <- if (k==2) hfPaths$pB[10:2000] else uPaths$pB[10:2000]
    lines(10:2000, yLine, col=cols[k], lwd=myLwd)
  }
}

par(xpd = NA)
for (k in 1:3) {
  yCor <- if (k==1) cPaths$pB[2000] else if (k==2) hfPaths$pB[2000] else uPaths$pB[2000]

  text(2000, yCor,
       labs[k], pos = 4, offset = 0.85, col = cols[k],
       cex = 0.85*myCex)
}


text(N, par("usr")[4] + 0.0001, paste("N =", N), col = "gray40", cex = myCex)
text(11, alpha + 0.0025, bquote(alpha == .(alpha)), adj = c(0, 0.5),
     col = "gray40", cex = myCex)
dev.off()


# "power" plot ----

adjustTextHeight <- 0.04

pdf("conformalPower.pdf", width = 14, height = 7)
par(mai = c(1.2, 1.5, 0.2, 1.7), cex.axis = 4/3, las = 1)

plot(NA, xlim = c(10, 2000), ylim = c(0, 1), log = "x",
     axes = FALSE, xlab = "", ylab = "")
axis(1, at = c(10, 20, 50, 100, 200, 500, 1000, 2000), cex.axis=cexFactor*myCex)
axis(2, at = seq(0, 1, 0.1), labels = seq(0, 100, 10), cex.axis=cexFactor*myCex)
mtext("Time horizon N", side = 1, line = 3.8, cex = cexFactor*myCex)

mtext(latex2exp::TeX(r"($b p_{b,N} + V_{b,N}$ in %)"),
      side = 2, line = 4, las=0,
      cex = cexFactor*myCex,
      adj=0.5, padj=-0.4)


# lines(c(1, 2000), c(alpha, alpha),
#       lwd=myLwd, col="lightgrey", lty=2)

lines(c(N, N), c(0, 1),
      lwd=myLwd, col="grey", lty=2)


# draw universal first so committed ends up on top, as in the original
for (k in 1:3) {
  if (k==1) {
    lines(10:200, cPaths$power[10:200], col=cols[k], lwd=myLwd)
    lines(201:2000, cPaths$power[201:2000], col="darkorchid3",
          lwd=myLwd, lty=2)
  } else {
    yLine <- if (k==2) hfPaths$power[10:2000] else uPaths$power[10:2000]

    tempLwd <- if (k==2) 1.5*myLwd else myLwd

    lines(10:2000, yLine, col=cols[k], lwd=tempLwd)
  }
}

par(xpd = NA)
for (k in 1:3) {
  yCor <- if (k==1) cPaths$power[2000]+adjustTextHeight else if (k==2) hfPaths$power[2000] else uPaths$power[2000]-adjustTextHeight

  text(2000, yCor,
       labs[k], pos = 4, offset = 0.85, col = cols[k],
       cex = 0.85*myCex)
}


text(N, par("usr")[4] + 0.0001, paste("N =", N), col = "gray40", cex = myCex)
# text(11, alpha + 0.0025, bquote(alpha == .(alpha)), adj = c(0, 0.5),
#      col = "gray40", cex = myCex)
dev.off()


# N=150 -----
chop <- processNullSamplePaths(eValueMatrixNull[, 1:150], b=9.918542)
chop$b*chop$pB+chop$vB

chopS <- processNullSamplePaths(eValueMatrixNull[, 1:150], b=20)
chopS$b*chopS$pB+chopS$vB

processNullSamplePaths(eValueMatrixNull, b=15)$pB
processNullSamplePaths(eValueMatrixNull, b=5)$pB

# N=100 -----

chan <- processNullSamplePaths(eValueMatrixNull[, 1:100], b=8.116551)
chanS <- processNullSamplePaths(eValueMatrixNull[, 1:100], b=20)

chan$b*chan$pB+chan$vB
chanS$b*chanS$pB+chanS$vB



# Logan 500  ------

aap <- processNullSamplePaths(eValueMatrixNull, b=15.2, nEnd=500)
aapF <- processNullSamplePaths(eValueMatrixNull, b=18.1, nEnd=500)
aapS <- processNullSamplePaths(eValueMatrixNull, b=20, nEnd=500)

aap$pB
aapF$pB
aapS$pB

aap$b*aap$pB+aap$vB
aapF$b*aapF$pB+aapF$vB
aapS$b*aapS$pB+aapS$vB

# Logan 2000  ------
brie <- processNullSamplePaths(eValueMatrixNull, b=15.2)
brieF <- processNullSamplePaths(eValueMatrixNull, b=18.1)
brieS <- processNullSamplePaths(eValueMatrixNull, b=20)

brie$pB
brieF$pB
brieS$pB

# N=2000 -------
bNAlpha2K <- uniroot(function(x){
  processNullSamplePaths(eValueMatrixNull, b=x)$pB-0.05
}, c(10, 20))

bAlpha2K <- uniroot(function(x){
  kop <- processNullSamplePaths(eValueMatrixNull, b=x)
  kop$pB+kop$rB/kop$b-0.05
}, c(10, 20))


bNAlpha2K$root
bAlpha2K$root

cTerms2K <- processNullSamplePaths(eValueMatrixNull, b=17.6908)
round(cTerms2K$pB, 5)
round(cTerms2K$pB+cTerms2K$rB/cTerms2K$b, 5)
round(cTerms2K$b*cTerms2K$pB+cTerms2K$vB, 5)

hfTerms2K <- processNullSamplePaths(eValueMatrixNull, b=17.75044)

round(hfTerms2K$pB, 5)
round(hfTerms2K$pB+hfTerms2K$rB/hfTerms2K$b, 5)
round(hfTerms2K$b*hfTerms2K$pB+hfTerms2K$vB, 5)
