## use a separate file for forward reads and reverse reads
## n is the number of positions to plot

plotmutfreq_forward <- function(file=NULL, n=100) {
	freq <- read.table(file,head=T)
	pos <- c(1:n)
	plot(pos,freq$CT1,type="l",col=2,lty=2, ylab="Substitution freq.",xlab="5'-3'",ylim=c(0,0.01))
	points(pos,freq$CA1,type="l", col=7)
	points(pos,freq$CG1,type="l",col=8)
	points(pos,freq$GA1,type="l",col=3,lty=2)
	points(pos,freq$GC1,type="l",col=1)
	points(pos,freq$GT1,type="l",col=2)
	points(pos,freq$TA1,type="l",col=5)
	points(pos,freq$TC1,type="l",col=1,lty=2)
	points(pos,freq$TG1,type="l",col=4)
	points(pos,freq$AT1,type="l",col=3)
	points(pos,freq$AG1,type="l",col=5,lty=2)
	points(pos,freq$AC1,type="l",col=6)
	legend("top",c("C>T","C>A","C>G","G>A","G>C","G>T","T>A","T>C","T>G","A>T","A>G","A>C"),
	col=c(2,7,8,3,1,2,5,1,4,3,5,6),lty=c(2,1,1,2,1,1,1,2,1,1,2,1),lwd=2,ncol=6,bty="n")
}
	

plotmutfreq_reverse <- function(file=NULL, n=100) {
	freq <- read.table(file,head=T)
	pos <- c(1:n)
	plot(pos,freq$CT2,type="l",col=2,lty=2, ylab="Substitution freq.",xlab="3'-5'",ylim=c(0,0.01))
	points(pos,freq$CA2,type="l", col=7)
	points(pos,freq$CG2,type="l",col=8)
	points(pos,freq$GA2,type="l",col=3,lty=2)
	points(pos,freq$GC2,type="l",col=1)
	points(pos,freq$GT2,type="l",col=2)
	points(pos,freq$TA2,type="l",col=5)
	points(pos,freq$TC2,type="l",col=1,lty=2)
	points(pos,freq$TG2,type="l",col=4)
	points(pos,freq$AT2,type="l",col=3)
	points(pos,freq$AG2,type="l",col=5,lty=2)
	points(pos,freq$AC2,type="l",col=6)
	legend("top",c("C>T","C>A","C>G","G>A","G>C","G>T","T>A","T>C","T>G","A>T","A>G","A>C"),
	col=c(2,7,8,3,1,2,5,1,4,3,5,6),lty=c(2,1,1,2,1,1,1,2,1,1,2,1),lwd=2,ncol=6,bty="n")
}
