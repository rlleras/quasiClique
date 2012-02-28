for (e in commandArgs()) {
    ta <- strsplit(e, '=');
    if (ta[[1]][1] == '--fam') {
        fam <- ta[[1]][2];
    }
	if (ta[[1]][1] == '--input') {
		input_filename <- ta[[1]][2];
	}
	if (ta[[1]][1] == '--scorename') {
		scorename <- ta[[1]][2];
	}
}
# fam = 'purine' or something
# scorename = 'pscore' or 'rank'
# input = eval file sorted by score

x<-read.table(input_filename, sep='\t', row.names=1)
N <- length( row.names(x) );

# we round the scores to integers and repeat it twice (once for tp, once for fp)
round0rep2 <- function(x) { rep(round(as.double(x), 0), 2) }
row_names_alt_list <- c()

scores <- c();
# use boxplots, already sorted by increasing pscore/rank
V <- list();
for (i in 1:N) {
  tp <- as.numeric( strsplit( as.character(x[i,6]), split="\\," )[[1]]);
  fp <- as.numeric( strsplit( as.character(x[i,7]), split="\\," )[[1]]);
  if (scorename=='pscore') { scores <- c( scores, round0rep2(x[i,2])); }
  else if (scorename=='pscore-RAxML') { scores <- c( scores, round0rep2(x[i,3])); }
  else { scores <- c( scores, round0rep2(x[i,1])); }
  new_name <- sub("\\.fna", "", row.names(x)[i]);
  new_name <- sub("\\.motif", "", new_name);
  row_names_alt_list <- c(row_names_alt_list, paste(new_name,'(',scores[2*i],')',sep=''), ''); 
  V[[i*2-1]] <- tp;
  V[[i*2]] <- fp;
}
png(paste(fam,'.',scorename,'_vs_CMhits.png',sep=''),width=40*N);
par(mar=par()$mar+c(12,0,0,0));
boxplot(V,names=row_names_alt_list,col=c("blue","yellow"),main=paste(fam,' (',scorename,')',sep='') ,ylab="CM hit score",las=2);
legend("topright", fill=c("blue","yellow"), legend=list("TPs","FPs"));
dev.off()

#min_y <- Inf;
#max_y <- -Inf;
#pscores <- c();
#for (i in 1:N) {
#    tp <- as.numeric( strsplit( as.character(x[i,4]), split="\\," )[[1]]);
#    fp <- as.numeric( strsplit( as.character(x[i,5]), split="\\," )[[1]]);
#    min_y <- min(min_y, min(tp), min(fp));
#    max_y <- max(max_y, max(tp), max(fp));
#    pscores <- c( pscores, as.numeric(x[i,2]) );
#}
#
#i <- 1;
#tp <- as.numeric( strsplit( as.character(x[i,4]), split="\\," )[[1]]);
#fp <- as.numeric( strsplit( as.character(x[i,5]), split="\\," )[[1]]);
#png('test.png')
#plot( rep(pscores[1], length(tp)), tp, xlim=c(min(pscores),max(pscores)), ylim=c(min_y,max_y), col='blue' );
#points( rep(pscores[1], length(fp)), fp, pch=17, col='black' );
#for (i in 2:N) {
#    tp <- as.numeric( strsplit( as.character(x[i,4]), split="\\," )[[1]]);
#    fp <- as.numeric( strsplit( as.character(x[i,5]), split="\\," )[[1]]);
#    points( rep(pscores[i], length(tp)), tp, col='blue');
#    points( rep(pscores[i], length(fp)), fp, pch=17, col='black' );
#}
#dev.off()
