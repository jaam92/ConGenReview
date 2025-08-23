setwd("~/Documents/USC/projects/jazlyn_dfe/")
a = read.table("results/20_0.05_0.01.txt")
b = read.table("results/8_0.05_0.01.txt")

big_a1 = mean(abs(a$V1 / a$V5))
sma_a1 = mean(abs(b$V1 / b$V5))
big_a2 = mean(abs(a$V2 / a$V6))
sma_a2 = mean(abs(b$V2 / b$V6))

big_d1 = mean(abs(a$V7 / a$V11))
sma_d1 = mean(abs(b$V7 / b$V11))
big_d2 = mean(abs(a$V8 / a$V12))
sma_d2 = mean(abs(b$V8 / b$V12))

pdf("error.pdf", 10, 5)
par(cex=1.3, mfrow=c(1,2))

barplot(matrix(c(big_a1,sma_a1,big_a2, sma_a2), 2, 2), 
        beside=T, ylab="Relative Mean Absolute Error",
        names=c("alpha", "beta"), col=c("grey80", "grey20"))
mtext("A", side = 3, line = 1, adj = -0.1, cex = 1.5, font = 2)
legend("topleft", c("n=20", "n=8"), fill=c("grey80", "grey20"), bty='n')
barplot(matrix(c(big_d1,sma_d1,big_d2, sma_d2), 2, 2), 
        beside=T, ylab="Relative Mean Absolute Error",
        names=c("Contraction size", "Contraction time"), col=c("grey80", "grey20"))
mtext("B", side = 3, line = 1, adj = -0.1, cex = 1.5, font = 2)

dev.off()