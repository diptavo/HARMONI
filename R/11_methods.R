###############################################################################
# MS-BPWAS  -  11_methods.R
# S3 methods: print, summary, plot
###############################################################################

#' @export
print.msbpwas <- function(x, ...) {
 cat("MS-BPWAS results\n")
 cat(sprintf("  Subtypes : %d (%s)\n", x$S, toString(x$subtype_names)))
 cat(sprintf("  Proteins : %d\n", x$P))
 cat(sprintf("  Method   : %s | Prior: %s | Overlap: %s\n",
     x$control$method, x$control$prior_type, x$control$overlap_correction))
 cat(sprintf("  pi0      : %.4f\n", x$hyperparams$pi0))
 cat(sprintf("  xi       : %s\n", paste(sprintf("%.3f",x$hyperparams$xi),collapse=", ")))
 cat(sprintf("  BFDR<%.2f: %d discoveries\n",
     x$control$bfdr_threshold, sum(x$bfdr < x$control$bfdr_threshold)))
 invisible(x)
}

#' @export
summary.msbpwas <- function(object, top_n=20, ...) {
 x <- object
 cat("======== MS-BPWAS Summary ========\n\n")
 cat(sprintf("Learned null rate: %.1f%%\n", 100*x$hyperparams$pi0))
 cat("Per-subtype inclusion:\n")
 for(s in seq_along(x$subtype_names))
   cat(sprintf("  %-15s %.1f%%\n", x$subtype_names[s], 100*x$hyperparams$xi[s]))

 if(!is.null(x$hyperparams$pi_mix)) {
   cat("\nMixture weights:\n")
   for(l in seq_along(x$control$sigma2_grid))
     cat(sprintf("  sigma=%.4f  w=%.3f\n",
         sqrt(x$control$sigma2_grid[l]), x$hyperparams$pi_mix[l]))
 }

 cat("\nDiscoveries:\n")
 for(t in c(0.01,0.05,0.10,0.20))
   cat(sprintf("  BFDR < %.2f : %d\n", t, sum(x$bfdr < t)))

 cat("\nMAP config distribution:\n")
 tb <- sort(table(x$map_config), decreasing=TRUE)
 for(i in seq_len(min(10,length(tb))))
   cat(sprintf("  %-30s %d\n", names(tb)[i], tb[i]))

 cat(sprintf("\nTop %d proteins:\n", top_n))
 ord <- order(x$p_null); n <- min(top_n, x$P)
 cat(sprintf("  %-15s P(null) BFDR   MAP%-26s p_cond     PIPs\n","Protein",""))
 for(i in seq_len(n)) { p<-ord[i]
   cat(sprintf("  %-15s %.4f  %.4f %-28s %.2e  %s\n",
     x$protein_ids[p], x$p_null[p], x$bfdr[p], x$map_config[p],
     x$p_conditional[p], paste(sprintf("%.2f",x$pip[p,]),collapse=" ")))
 }
 invisible(x)
}

#' @export
plot.msbpwas <- function(x, type=c("manhattan","pip_heatmap",
                                     "convergence","config_bar"), ...) {
 type <- match.arg(type)

 if (type=="manhattan") {
   y <- -log10(pmax(x$p_null, 1e-300))
   plot(seq_along(y), y, pch=16, cex=0.4,
        col=ifelse(x$bfdr<x$control$bfdr_threshold,"red","grey60"),
        xlab="Protein index", ylab=expression(-log[10](P(null))),
        main="MS-BPWAS Manhattan", ...)
   abline(h=-log10(0.10), lty=2, col="blue")
   top <- order(y,decreasing=TRUE)[1:min(5,length(y))]
   text(top, y[top], x$protein_ids[top], pos=3, cex=0.65, col="red")
 }
 if (type=="pip_heatmap") {
   n <- min(30, sum(x$p_null<0.50))
   if(n<2){ message("Too few hits"); return(invisible()) }
   o <- order(x$p_null)[1:n]; mat <- x$pip[o,,drop=FALSE]
   image(t(mat[nrow(mat):1,]), axes=FALSE,
         col=colorRampPalette(c("white","orange","red3"))(100),
         main="PIP heatmap")
   axis(1, seq(0,1,length=ncol(mat)), x$subtype_names, las=2)
   axis(2, seq(0,1,length=nrow(mat)), rev(x$protein_ids[o]), las=1, cex.axis=.55)
 }
 if (type=="convergence") {
   h <- x$convergence$elbo_hist %||% x$convergence$Q_hist
   if(is.null(h)){ message("No history for MCMC"); return(invisible()) }
   plot(h, type="l", lwd=2, col="steelblue",
        xlab="Iteration", ylab=if(!is.null(x$convergence$elbo_hist)) "ELBO" else "Q",
        main=paste(toupper(x$control$method),"convergence"))
 }
 if (type=="config_bar") {
   tb <- sort(table(x$map_config[x$p_null<0.50]), decreasing=TRUE)
   if(!length(tb)){ message("No non-null hits"); return(invisible()) }
   par(mar=c(9,4,3,1))
   barplot(tb[seq_len(min(15,length(tb)))], las=2, col="steelblue",
           main="MAP configs (non-null)", ylab="Count")
 }
 invisible(x)
}
