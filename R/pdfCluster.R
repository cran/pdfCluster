# setClasses
setClass("qdelaunay", representation(nb="list", new="matrix", call="call"))

setClass("summary.qdelaunay", representation(nregions="numeric", nnb="numeric", nblist="list"))

setClass("kepdf", representation(x="matrix", eval.points="matrix", estimate="numeric", h="numeric", h.weights="numeric", call="call"))

setClass("summary.kepdf",representation(obj.class="character", mode="numeric", props="numeric", indices="list"))

setOldClass("dendrogram")   

setClass("pdfCluster", representation(call="call", x="matrix", estimate="numeric", h="vector", hmult="numeric", nc="list",  
              cluster.cores="ANY", tree="ANY", ctrl="logical", stages="ANY", clusters="ANY"))
			  
setClass("pdfSilhouette", representation(x="matrix", prior="numeric", dbs="numeric", clusters="numeric", nc="numeric", stage="ANY", call="call"))

setClass("summary.pdfCluster", representation(obj.class="character", cluster.cores="numeric", clusters="ANY", tree="ANY"))

setClass("summary.pdfSilhouette", representation(obj.class="character", dbs.groups = "list"))	


##############################################################
################################################################
################################################################	
# functions

kepdf <- function(x, h=h.norm(x), eval.points=x, h.weights=1){
	x <- as.matrix(x) 
	eval.points <- matrix(eval.points,ncol=ncol(x))
	nx <- as.integer(nrow(x))
	ndim <- as.integer(ncol(x))
	y <- matrix(as.double(x), nx, ndim)
	neval <- as.integer(nrow(eval.points))
	eval.points <- matrix(as.double(eval.points), neval, ndim)
	f <- as.double(rep(0,neval))
	if(length(h.weights==1)) h.weights <- rep(h.weights, nx)
	hm <- outer(h.weights, h)
	hm <- matrix(as.double(hm), nx, ndim)
	pdf <- .Fortran("kepdf", y, hm, eval.points, nx, ndim, neval, f, PACKAGE="pdfCluster")
	new("kepdf", x=x, eval.points=eval.points, estimate=pdf[[7]], h=h, h.weights=h.weights, call=match.call())
}

summary.kepdf <- function (object, ..., props=c(75,50,25)) {
	if(dim(object@eval.points)[1]==dim(object@x)[1] && all(object@x==object@eval.points)) pdf.obs <- object@estimate else pdf.obs <- kepdf(object@x, h=object@h, eval.points=object@x, h.weights=object@h.weights)@estimate
		hts <- quantile(rep(pdf.obs, object@h.weights), prob = (100 - props)/100)
	modepos=which.max(pdf.obs)				
	indices <- list()
	for(i in 1:length(props)){
		indices[[i]] <- (which(pdf.obs>hts[i]))
	}				
	names(indices) <- paste(as.character(props),"%",sep="")				
    new("summary.kepdf", obj.class=class(object)[1], mode=modepos, props=props, indices=indices)
}

plot.kepdf1d <- function(x, y, eval.points=NULL, n.grid=25, data=NULL, add=FALSE, main	= NULL, xlab=NULL, ylab=NULL, col=NULL, col.data=2, type="l",...){
	##########################################
	##########################################
	# comments by GM 03/02/2011
	# plot of univariate kernel density estimate 
	# DEFAULT: density is evaluated on a grid defined on the range of evaluation points of the density
	##########################################
	##########################################
	if(is.null(main))	main="Kernel density estimate of data"
	if(is.null(xlab))	xlab="eval.points"
	if(is.null(ylab))	ylab="density function"
	if(is.null(eval.points)) {
	datarange <- apply(x@eval.points, 2, extendrange)
	eval.points <- seq(datarange[1, 1], datarange[2, 1], length=n.grid)
	} 
	pdf <- kepdf(as.vector(x@x), h=x@h, eval.points=eval.points, h.weights=x@h.weights)
	pdf.est.sort <- pdf@estimate[order(eval.points)]
	eval.points.sort=sort(eval.points)
	if(is.null(col)) col=1
	if (add){lines(eval.points.sort, pdf.est.sort, col=col,...)} else {
	plot(eval.points.sort, pdf.est.sort, type=type, main=main, xlab=xlab, ylab=ylab, col=col, ...)}
	if (!is.null(data)) rug(data, col=col.data)
	invisible(list(eval.points=as.vector(eval.points), estimate=pdf@estimate))
}

# lines.kepdf <- function(x, eval.points= NULL, ...){
	# if(is.null(eval.points)) {
	# datarange <- apply(x@eval.points, 2, extendrange)
	# eval.points <- seq(datarange[1, 1], datarange[2, 1], length=n.grid)
	# } else {
	# pdf <- kepdf(as.vector(x@x), h=x@h, eval.points=eval.points, h.weights=x@h.weights)}
	# pdf.est.sort <- pdf@estimate[order(eval.points)]
	# eval.points.sort=sort(eval.points)
	# lines(eval.points.sort, pdf.est.sort, ...)
	# invisible(list(eval.points=as.vector(eval.points), estimate=pdf@estimate))
# }

plot.kepdf2d <- function (x, y, eval.points=NULL, n.grid=c(25,25), 
	data=NULL, add=FALSE, main=NULL, xlab=NULL, ylab=NULL, zlab=NULL, col=NULL, col.data=2, props=c(75,50,25), method="contour", ticktype = "detailed", ...) {
	if(is.null(main)) main="Kernel density estimate of data"
	if(is.null(eval.points)) {  
	datarange <- apply(x@eval.points, 2, extendrange)
	xgrid <- seq(datarange[1, 1], datarange[2, 1], length=n.grid[1]) # sort(x@eval.points[,1])
	ygrid <- seq(datarange[1, 2], datarange[2, 2], length=n.grid[2]) # sort(x@eval.points[,2])
	} else { 
	xgrid=sort(eval.points[, 1])
	ygrid=sort(eval.points[, 2])
	}
	nx <- length(xgrid)
	ny <- length(ygrid)
	grid.points <- as.matrix(expand.grid(xgrid,ygrid))
	pdf <- kepdf(x@x, h=x@h, eval.points=grid.points, h.weights=x@h.weights)
	lab <- colnames(x@x) 
	if(is.null(xlab)) if(is.null(lab)) xlab <- "Var 1" else xlab <- lab[1]
	if(is.null(ylab)) if(is.null(lab)) ylab <- "Var 2" else ylab <- lab[2]
	if(method=="image"){
		if(is.null(col)) col=heat.colors(12)		
		image(xgrid, ygrid, matrix(pdf@estimate, nx,ny), add = add, main=main, col=col, xlab=xlab, ylab=ylab, ...)
	} else if(method=="perspective"){
		if (is.null(zlab)) zlab="density function"
		if(is.null(col)) col="lightblue"		
		persp(xgrid,ygrid, matrix(pdf@estimate, nx,ny), main=main, xlab = xlab, ylab = ylab, zlab=zlab, col=col, ticktype = "detailed", ...)
	} else {
	if(method!="contour") warning("\"method\" should be one of \"contour\", \"image\", \"perspective\". A \"contour\" plot has been created.")
	if(!add) plot(grid.points, type = "n", main=main, xlab = xlab, ylab = ylab, ...)
	if(dim(x@eval.points)[1]==dim(x@x)[1] && all(x@x==x@eval.points)) pdf.obs <- x@estimate else pdf.obs <- kepdf(x@x, h=x@h, eval.points=x@x, h.weights=x@h.weights)@estimate
	hts <- quantile(rep(pdf.obs, x@h.weights), prob = (100 - props)/100)
	if(is.null(col)) col=1		
	for (i in 1:length(props)) {
	scale <- props[i]/hts[i]
	contour(xgrid, ygrid, matrix(pdf@estimate, nx,ny) * scale, level = hts[i] *scale, add = TRUE, main=main, col=col, ...)
    }
	if(!is.null(data)) points(data, col=col.data, ...)#aggiunge i punti	
	}
	invisible(list(eval.points = pdf@eval.points, estimate = pdf@estimate))
	}		
		
plot.kepdfdd <- function (x, y, ..., main=NULL, method="contour", text.diag.panel = NULL, gap = 0.5) {
	textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, y, txt, cex = cex, font = font)
	localAxis <- function(side, x, y, xpd, bg, col = NULL, main,  oma, ...) {
    if (side%%2 == 1) 
        Axis(x, side = side, xpd = NA, ...)
        else Axis(y, side = side, xpd = NA, ...)
    }
	localPlot <- function(..., oma, font.main, cex.main) plot.kepdf2d(...)
	dots <- list(...)
    nmdots <- names(dots)
	if (!is(x, "kepdf")) stop("argument to 'pairs.kepdf must be of class \"kepdf\"")
	nc <- ncol(x@x)
    if (is.null(text.diag.panel)) if(!is.null(lab <- colnames(x@x))) text.diag.panel <- lab else text.diag.panel <- paste(rep("V", nc),1:nc,sep="")
	oma <- if ("oma" %in% nmdots) 
    dots$oma else NULL
    #main <- if ("main" %in% nmdots) 
    #    dots$main  else "Pairwise marginal kernel density estimates of data"
    if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main)) 
    oma[3L] <- 6
    }    
	#if (!is.null(main)) main="Pairwise marginal kernel density estimates of data"
    opar <- par(mfrow = c(nc, nc), mar = rep(c(gap,gap/2),each=2), oma=oma)
    on.exit(par(opar))
	out <- list()
	count <- 1
    for (i in 1L:nc){ 
		for (j in (1L:nc)) 
			if(i==j) {
		plot(1,type="n",xlab="",ylab="",axes=FALSE)
		text(1,1,text.diag.panel[i],cex=2)
		box()} else {
        marg <- kepdf(x@x[,c(j,i)], h=x@h[c(i,j)])
		out[[count]] <- localPlot(x=marg, xlab = "", ylab = "", main="", axes=FALSE, method=method, ...)
        count <- count + 1     
			if (method!="perspective"){
			if(i==1) {axis(3); if(j==nc) axis(4)}
			if(j==1) {axis(2); if(i==nc) axis(1)}
			box()} #else if (i == j) {j=j+1}
                }
		count <- count + 1}			
    par(new = FALSE)
	if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots) dots$font.main
        else par("font.main") 
        cex.main <- if ("cex.main" %in% nmdots) dots$cex.main
        else par("cex.main") 
        mtext(main, side=3, TRUE, line = 3, outer = TRUE, at = NA, cex = cex.main, font = font.main, adj= 0.5)
    }	
    invisible(out)
	}

plot.kepdf <- function(x, y,...){
	if(is.vector(x@x)) d <- 1 else d <- dim(x@x)[2]
	if (d==1) plot.kepdf1d(x, y, ...) else {
	if (d==2) plot.kepdf2d(x, y, ...) else 
	plot.kepdfdd(x, y, ...) }
	}

require(geometry)
   
qdelaunay <- function(x, Q="QJ", mult=1.1) {
   # ----------------------------------------------------------------------
   # Crea la triangolazione e elenca i VICINI PER OGNI PUNTO
 
   # REQUIRES geometry
   # Options implementate per qdelaunay.exe:
   # Qt  - tutte le regioni sono simplicial (in 2D triangoli)
   # Qj  - joggle input tutte le regioni sono simplicial (in 2D triangoli)
   # Fv  - l'output e' la lista dei punti per ogni triangolo
   #       (i punti sono numerati da 0 a n-1 come nel data set)
 
   # OUTPUT - per ogni punto lista dei vicini, punti per cui esiste un arco
   # ----------------------------------------------------------------------
   x <- as.matrix(x)
   x0 <- x
   sd <- sqrt(diag(var(x)))
   m  <- apply(x,2,mean)
	#GM 17/01/2011 standardizzo i dati
   x.scaled <- t((t(x)-m)/sd)
   ncx <- dim(x.scaled)[2]
   xrange <- apply(x.scaled,2,range)
	#GM 17/01/2011 inizializzo matrice che contiene nuovi punti, vertici dell'involucro convesso (credo)
   xnew <- matrix(0,ncol=ncx,nrow=2^ncx)
   for (i in 1:ncx)
      xnew[,i] <- rep(xrange[,i],each=(2^(i-1)),times= (2^(ncx-i)))
	xnew.id <- rep(c(FALSE,TRUE),c(dim(x.scaled)[1],dim(xnew)[1]))
    #GM 17/01/2011 xnew contiene i vertici, corrispondenti a max e min di x molt per 1.1
   xnew <- xnew*mult
   x.scaled <- rbind(as.matrix(x.scaled),xnew)
   nrx <- dim(x.scaled)[1]
   if (Q=="QT") xt <- delaunayn(x.scaled, options="Qt")
   if (Q=="QJ") xt <- delaunayn(x.scaled, options="QJ")
   on.exit(unlink("qhull_out.txt"))
   x.nb <- vector(mode="list", length=nrx)
   for (i in 1:nrx) {x.nb[[i]] <- 
as.integer(sort(unique(as.vector(xt[(xt==i)%*%rep(1,ncx+1)>0,]), FALSE))) 
   }
rownames <- as.character(1:length(x.nb))
#if (ncx > 2 & plot) { plot <- FALSE; warning("Plot of Delaunay triangulation is available for 2-d data only")}
x.out <- new("qdelaunay", nb=x.nb, new=xnew, call=match.call())
#  if (plot)  plot(x=x.out, data=x)
  x.out
  }

	summary.qdelaunay <- function (object,...) {   
nregions <- length(object@nb)
nnb <- unlist(lapply(object@nb,length))
nblist <- object@nb
new("summary.qdelaunay", nregions=nregions, nnb=nnb, nblist=nblist)
}

#	plot.qdelaunay <- function(x, y, data){
#   # ----------------------------------------------------------------------
#   # Plot di una triangolazione - solo in due dimensioni
#   # 
#   # INPUT
#   #  data    - matrice dei dati
#   #  y       - controllo, solo per creare metodo plot
#   #  x       - triangolazione dei dati
#   # ----------------------------------------------------------------------
#   x.nb <- x@nb
#   x1 <- apply(as.matrix(data),2,scale)
#   x1 <- rbind(x1,x@new)
#   if (dim(x1)[2]>2) plot <- FALSE else plot=TRUE
#   if (plot) {
#      if (is.null(colnames(data))) colnames(x1) <- c("x","y")
#      plot(x1 , main="Delaunay triangulation",cex=0.5)
#      for (j in 1:length(x.nb)) {
#        for (i in 1:length(x.nb[[j]])) 
#                    lines(rbind(x1[j,],x1[x.nb[[j]],][i,]))
#      }
#   }
#}  


num.con.multi <- function(x, xd, pn=0.9, profile=TRUE, profile.grid=25, 
                  correct=FALSE, profile2=FALSE, qnv=0, plot.it=FALSE, Q="QJ", mult=1.1){
   # ----------------------------------------------------------------------
   # computes number of connected regions by using the depth-first search algorithm

      # INPUT
   #  x         data matrix --> used as input for the Delaunay triangulation, internally performed.  
   #  			The considered output of executing qdelaunay is x.nb, the structure of the neightbors of x
   #  xd        density estimate 
   #  pn        proportion of data to be included
   #  profile   evaluates connected regions on a grid; pn is ignored

   # OUTPUT   - list with elements
   #  nc      - number of connected for each value of pn
   #  p       - probab. pn [AA]
   #  id      - matrix of group labels at each point
   #  pq      - (CTRL) for each pn gives the corresponding quantile
   # ----------------------------------------------------------------------

   nb <- qdelaunay(x, Q=Q, mult=mult)

   x.nb <- nb@nb

   #require(prabclus)
   n   <- length(xd@estimate)
   ngrid <- profile.grid
   if (profile)   pn <- seq(0,1,length=ngrid+2)[-c(1,ngrid+2)]
   # p <- n * (1-pn)      								#GM 18/01/2011 non mi sembra sia usato dopo
   # p.trunc <- trunc(p) 								#GM 18/01/2011 non mi sembra sia usato dopo
   fest <- c(xd@estimate, rep(0,dim(nb@new)[1]))
   # fhat <- sort(xd@estimate)  
   # qn <- fhat[p.trunc] + (p - p.trunc) * (fhat[p.trunc + 1] - fhat[p.trunc])
   #
   # 30-06-2004: ripristiamo istruzione commentata che segue (ma modificata)
   #  qn <- quantile(xd@estimate, pn)
   qn <- as.numeric(quantile(xd@estimate, 1-pn))
   if (qnv != 0) qn <- qnv
        colnames(nb@new) <- colnames(xd@eval.points)
   x <- rbind(xd@eval.points,nb@new)
   n1  <- n + dim(nb@new)[1] 
   n0  <- length(x.nb) 
   #GM 18/01/2011 xd.id is a logical matrix.  is the density below or above the threshold?  
   xd.id <- (matrix(fest , ncol=1)%*%(1/qn) ) < 1

   nc.nc                  <- rep(0,length(qn))
   names(nc.nc)           <- format(pn, digit=3)
   nc.id                  <- matrix(0,nrow=n0,ncol=length(qn))
   nc.id1                 <- matrix(0,nrow=n1,ncol=length(qn))
   nc.id2                 <- matrix(0,nrow=n,ncol=length(qn))

   colnames(nc.id)  <- pn
   #browser()
    ##############
			
	##############
    for (i in 1:length(qn)) {							    #GM 18/01/2011 at each section "i" of the density
	ni    <- sum(xd.id[,i]) 								#GM 18/01/2011 ni: how many values are below the threshold
      if (ni<n0) { 											#GM 18/01/2011 if alle the values are below the threshold
         x.nbc <- x.nb
			ind <- which(!xd.id[,i])
			for(j in ind)
			x.nbc[[j]] <- intersect(x.nb[[j]],ind)
			x.nbc[xd.id[,i]] <- as.integer(0)
		 nc                   <- find.nc(x.nbc)
         nc.nc[i]             <- nc$nc - ni
         # nc.id[,i]          <- nc@comp.id
         nc.id[xd.id[,i],i]   <- -1
         nc.id[!xd.id[,i], i] <- unclass(factor(nc$comp.id[!xd.id[,i]]))
         # (AA, prima era) : codes(factor(nc@comp.id[!xd.id[,i] ]  ))
      }
      if (ni==n0) {
         nc.nc[i]  <- 0
         nc.id[,i] <- -1
      }
      if (correct){
         tt  <- table(nc.id[nc.id[,i] != -1, i])
         ttc <- as.integer(names(tt)[tt<=1])
         ttp <- as.integer(names(tt)[tt>1])
         if (length(ttc)>0) 
             for (l in 1:length(ttc)) nc.id[nc.id[,i]==ttc[l],i] <- -1
         nc.nc[i] <- length(ttp)
      }
             nc.id2[,i] <- nc.id[1:n,i]
     if (correct){
         tt <- table(nc.id2[nc.id2[,i]!=-1,i])
         ttc <- as.integer(names(tt)[tt<=1])
         ttp <- as.integer(names(tt)[tt>1])
         if (length(ttc)>0) 
             for (l in 1:length(ttc)) nc.id2[nc.id2[,i]==ttc[l],i] <- -1
         nc.nc[i] <- length(ttp)
         nc.nc[i] <- length(tt)   #GM 20/01/2011 a cosa serve? e soprattutto, a cosa serve questa correzione dal momento che poi viene annullata una riga sotto?
      }
      nc.nc[i] <- length(unique(nc.id2[,i]))-1 # 09-07-2004
   }
   if (profile & plot.it) {
      plot(c(0, pn, 1, 1), c(0, nc.nc, nc.nc[length(nc.nc)],0), 
           ylim=c(0,max(nc.nc)), yaxt="n", type="S", cex=0.7,
           xlab="Fraction of data points included",
           ylab="m(p)")
      axis(2,at=c(0:max(nc.nc)))
      #  lines(pn,nc.nc)
   }
   if (profile2 & plot.it) {
      par(mfrow=c(1,2))
      plot(pn,nc.nc,xlim=c(0,1),ylim=c(0,max(nc.nc)),cex=0.7, pch=3,
                xlab="Percentuale oss. incluse",ylab="Numero di gruppi")
      lines(pn,nc.nc)
      plot(qn/max(fest),nc.nc,xlim=c(0,1),ylim=c(0,max(nc.nc)),cex=0.7,
                 pch=3, xlab="Density level", ylab="Number of groups")
      lines(qn/max(fest),nc.nc)
      par(mfrow=c(1,1))
   }

   list(nc=nc.nc, p=pn, id=nc.id2, pq=cbind(p=pn,q=qn))
}

num.con.uni <- function(x, xd, pn = 0.9, profile = TRUE, profile.grid=25,
						correct = FALSE, profile2 = FALSE, qnv = 0, plot.it = FALSE){
	x1 <- sort(x)
	f <- xd@estimate[order(x)]
	riordino <- rank(x)
	ngrid <- profile.grid
		if (profile)   pn <- seq(0,1,length=ngrid+2)[-c(1,ngrid+2)]
	qn <- as.numeric(quantile(xd@estimate, 1-pn))
	if (qnv != 0) qn <- qnv
	xd.id <- (matrix(f , ncol=1)%*%(1/qn) ) > 1
	xd.id1 <- rbind(FALSE,xd.id,FALSE)
	jumps <- apply(xd.id1,2,diff)
		which.jumps.ini <- which(rbind(jumps)>0,arr.ind=TRUE)
		which.jumps.fin <- which(rbind(jumps)<0,arr.ind=TRUE)
	num.con <- as.vector(table(which.jumps.ini[,2]))
	length.con <- split(which.jumps.fin[,1]-which.jumps.ini[,1],which.jumps.ini[,2])
	lab.con.fun <- function(num.con,length.con)rep(1:num.con,times=length.con)
	lab.con <- mapply(lab.con.fun,num.con,length.con)
	y <- as.matrix(xd.id)
	y[y] <- unlist(lab.con)
	y[y==0] <- -1
	nc.id <- as.matrix(y[riordino,])
	if (correct){
		correct.id <- function(vec){
		t <- table(vec)
		uni <- as.integer(names(t)[t==1])
		new <- vec
		new[match(uni,vec)] <- -1
		new}	
        nc.id <- apply(nc.id,2,correct.id)   #GM 20/01/2011 giusto?
      }
	#browser()  
	nc <- (apply(nc.id,2,max))
	nc[nc<0] <- 0
	names(nc) <- format(pn, digit=3)
		out <- list(nc=nc,p=pn,id=nc.id,pq=cbind(pn,qn))
	if (profile & plot.it) {
      plot(c(0, pn, 1, 1), c(0, nc, nc[length(nc)],0), 
           ylim=c(0,max(nc)), yaxt="n", type="S", cex=0.7,
           xlab="Fraction of included data points ",
           ylab="m(p)")
      axis(2,at=c(0:max(nc)))
      #  lines(pn,nc.nc)
	}
	 if (profile2 & plot.it) {
       par(mfrow=c(1,2))
       plot(pn,nc,xlim=c(0,1),ylim=c(0,max(nc)),cex=0.7, pch=3,
                 xlab="Fraction of included data points ",ylab="Number of groups")
       lines(pn,nc)
       plot(qn/max(xd@estimate),nc,xlim=c(0,1),ylim=c(0,max(nc)),cex=0.7,
                  pch=3, xlab="Density level",ylab="Number of groups")
       lines(qn/max(xd@estimate),nc)
       par(mfrow=c(1,1))
    }
	out
}


con2tree <- function(object, f, debug=FALSE){
  # object is the output of num.con()
  # f density estimate
  # 
  nc <- object$nc
  p <- object$p
  index <- which(diff(nc) != 0)
  ns <- length(index) # numero di punti di salto
  ps <- p[index]      # frazione di dati inclusi ai vari punti di salto
  lista <- list()     # lista degli insiemi connessi ai punti di salto
  if(ns > 0)
     for(j in 1:ns) lista[[j]] <-  as.vector(object$id[,index[j]])
  else
     lista[[1]] <-  as.vector(object$id[,ncol(object$id)])
  if(debug){cat("con2tree\n"); browser()}
  #...inizio della ex funzione groups(), qui in parte modificata
  K <- length(lista)
  gruppi <- rep(0,length(lista[[1]]))      
  M <- 0          
  insiemi <- list()
  # passo con k=1
  k <- 1
  sets <- lista[[k]]
  sdiff <- setdiff(unique(sets), -1)
  for (m in sdiff)
      gruppi[which(sets==m)] <- M <- M+1
  allocated <- (gruppi>0)
  insiemi[[k]] <- setdiff(unique(gruppi),0)
  if(debug){cat(">>> prima allocazione (k=1):\n")
    print(table(gruppi))
    cat(">>> insiemi[[1]] :", insiemi[[k]],"\n")         
    # browser()
  }
  # ciclo k=2,...,K
  while(k < K) {
    k <- k+1
    sets  <- lista[[k]]
    insieme <- list()
    for(m in 1:max(sets)) {
      set <- which(sets==m)
      new <- setdiff(set, which(allocated))
      if(length(new)>0){
        g.new <- unique(gruppi[intersect(set, which(allocated))])
        if(debug)  { cat(">>> g.new:", g.new,"\n");  browser() }
        if(length(g.new) == 0)  gruppi[set] <- M <- M+1           
        if(length(g.new) == 1)  gruppi[set] <- g.new 
        allocated <- (gruppi>0)
        }
      gg <- sort(setdiff(unique(gruppi[set]),0))
      if(length(gg)>0) insieme[[length(insieme)+1]] <- gg 
      if(debug){ print(table(gruppi));    browser()}
    }
    insiemi[[k]] <- insieme
    if(debug) { 
        cat("****** fine ciclo k-mo, k=", k,"\n"); 
        print(insieme)
        cat("***************************\n")
     }
  }
  if(!missing(f)){
    if(debug) cat("re-assign labels, according to pdf highest value\n")
    u <- unique(gruppi[rev(order(f))])
    g <- rep(0,length(f))
    u0 <- u[u>0]
    for(i in 1:max(gruppi)) g[gruppi==u0[i]] <- i
    gruppi <- g
    if(debug){
      cat("densita` massime dei gruppi 1...M,0:\n")
      for(m in c(1:M,0)) print(c(m, max(f[gruppi==m])))
     }
  } #..fine ex funzione groups()
  if(debug) cat("--------fine fase 1, inizia albero\n")
  salti <- diff(c(0,nc))
  salta.giu <- rev(which(diff(nc)[index]<0))
  altezza <- numeric(M)
  m <- 0
  salti.su <- salti[salti>0]
  while(m <M){
    m <- m+1
    # altezza[m] <- p[1] * (which(cumsum(salti.su) >= m)[1] - 1)
    r <- min(which(cumsum(salti.su) >= m))
    altezza[m] <- p[salti>0][r]
  }
  if(debug){ cat("altezza:",altezza,"\n");   browser()}
  bad.grid <- any(is.na(altezza))
  sotto.albero <- function(tree, k, set,  debug=FALSE){
      insieme <- insiemi[[salta.giu[k]]]      
      r <- 0
      branch <- list()
      if(debug){
        cat("*** inizia sotto.albero, livello k:", k,"\n")
        print(insieme)
        browser()
      }
     for(item0 in insieme){
        item <- intersect(set, unlist(item0))
        if(debug){
          cat("**** k, r:", k, r,"\n")
          cat("item:", item,"\n") 
          cat("set:", set,"\n") 
        }
        if(length(item)==1){
          r <- r+1
          u <- item
          attr(u,"members") <- 1
          attr(u,"height") <- altezza[item]
          attr(u,"label")  <- paste(as.character(item), " ", sep="")
          attr(u,"leaf")   <- TRUE
          branch[[r]] <- u

        }     
      if(length(item) > 1) {

        r <- r+1
        u <- sotto.albero(list(), k+1, item)
        attr(u,"members") <- length(unlist(item))
        attr(u,"height") <- max(ps[salta.giu[k+1]])
        attr(u,"label")  <- paste("{", paste(item, collapse=","),"}", sep="")
        attr(u,"leaf")   <- FALSE   
        branch[[r]] <- u
        # browser()
      }
      }
      if(debug){
        cat("***branch[k]:", k,"\n")
        print(branch)
        browser()
      }
      branch
    }        
  tree <- sotto.albero(list(), 1, 1:M,  debug=debug)
  attr(tree, "members") <- M
  attr(tree, "height") <-  max(ps[salta.giu[1]])
  attr(tree, "label") <- paste("{", paste(1:M, collapse=","),"}", sep="")
  tree <- list(tree)
  attr(tree, "members") <- M
  attr(tree, "height") <- 1
  if (M > 1) {attr(tree, "class") <- "dendrogram"
  ctrl <- FALSE } else ctrl <- TRUE 
  if(debug) browser()
  invisible(list(g=gruppi, sets=insiemi, tree=tree, bad=bad.grid, ctrl=ctrl))
 }

pdfCluster <- function(x, h=h.norm(x), hmult=.75, n.grid=min(50,nrow(as.matrix(x))), n.stage=5,  debug=FALSE){
  x <- data.matrix(x)
  pdf <- kepdf(x=x, h=h*hmult, eval.points=x)
  estimate <- pdf@estimate
  n <- nrow(x)
  if(n.grid > n){
    warning("n.grid too large, set equal to n")
    n.grid <- min(n.grid, nrow(x))
  }
  pn <- seq(0,1,length=n.grid+2)[-c(1,n.grid+2)]
  nc <- num.con(x, pdf, profile.grid=n.grid-1, correct=TRUE, plot.it=FALSE)
  struct <- con2tree(nc, estimate,  debug=debug)
  if(struct$bad)
    stop("The grid is too coarse: re-run with larger value of 'n.grid'")
  g <- struct$g
  g[struct$g==0] <- NA
  out <- new("pdfCluster", call=match.call(), x=data.matrix(x), estimate=estimate, h=h, hmult=hmult, nc=nc,  
              cluster.cores=g, tree=struct$tree, ctrl= struct$ctrl, stages=NULL, clusters=NULL)
  if(n.stage> 0) 
    {out <- pdfClassification(out, n.stage=n.stage)
	}  
  out
}
#summary.pdfCluster <- function(object, ...){
summary.pdfCluster <- function(object,...){
 new("summary.pdfCluster", obj.class=class(object)[1], cluster.cores=object@cluster.cores, clusters=object@clusters, tree=object@tree)
}

pdfClassification <- function(obj, h.funct="h.norm", n.stage=5){
  x <- obj@x
  if (missing(h.funct)) h0 <- obj@h else {
  h.fn <-  get(h.funct, inherit = TRUE)
  h0 <- h.fn(x)}
  hmult <- obj@hmult
  lista <- list()
  g <- obj@cluster.cores
  n <- length(g)
  g[is.na(g)] <- 0
  M <- max(g)
	if (obj@ctrl) {
	obj@cluster.cores[g==0] <- 1
	for (stage in 1:n.stage) lista[[stage]] <- rep(1, nrow(x))} else{
  for(stage in 1:n.stage){
  unallocated <- which(g==0)
  f <- matrix(NA, nrow=sum(g==0), ncol=M)
  for(m in 1:M){
    h <- h0 * hmult * (1 + (stage)/n.stage)
	f[,m] <- kepdf(x=x[g==m,], h=h, eval.points=x[g==0,])@estimate
  }
  ind.max <- t(apply(f, 1,order)[rev(1:M),])#at each line decreasing sort
  i1 <- ind.max[,1]
  i2 <- ind.max[,2]
  f1 <- diag(f[,i1])
  f2 <- diag(f[,i2])
  log.ratio <- log(f1)- log(f2)
  if(stage < n.stage){
    alti <- (log.ratio >= quantile(log.ratio, (n.stage-stage)/n.stage))
    # alti <- alti & (pmax(f1,f2)> 0.33 * rep(max(c(f1,f2)), length(alti)))
    }
  else
    alti <- rep(TRUE, sum(g==0))
  g.alti <-  ind.max[alti,1]
  for(m in 1:M){
    nuovi <- unallocated[alti][g.alti==m]
    g[nuovi] <- m
    }
  lista[[stage]] <- g
  }}
  #browser()
  unlista <- unlist(lista)
  unlista[unlista==0] <- NA
  obj@stages <- split(unlista,rep(1:n.stage,each=nrow(x)))
  obj@clusters <- lista[[n.stage]]
  obj
}

plot.pdfCluster <- function(x, y, which=1:3, stage=Inf, pch=NULL, col=NULL, ...){
  if (is.element(2, which) & x@ctrl) which <- setdiff(which, 2)
  w12 <- sort(which)[1:min(2,length(which))]
  if(setequal(1:2,w12)) par(mfrow=c(1,2))
  if(is.element(1, which)){
     nc <- x@nc$nc
     p <- x@nc$p
     plot(c(0,p,1,1), c(0,nc,nc[length(nc)],0), 
           ylim=c(0,max(nc)), yaxt="n", type="S", cex=0.7,
           xlab="fraction of data points",
           ylab="mode function", ...)
     axis(2,at=c(0:max(nc)))
  }
  if(is.element(2, which)) {
    plot(x@tree, center=TRUE, ...)
    title(sub="groups tree")
    }
  if(setequal(1:2,w12)) par(mfrow=c(1,1))
  if(length(setdiff(w12,3)) > 0 ) {
     cat("press <enter> to continue...")
     readline()
   }
  if(is.element(3, which))    {
    #browser()  
    stage <- min(stage, length(x@stages))
    if(stage==0) {g <- x@cluster.cores; g[is.na(x@cluster.cores)] <- 0}
    else {g <- x@stages[[stage]]; g[is.na(x@stages[[stage]])] <- 0} 
	x <- x@x
    M <-  length(table(g))
    if (is.null(pch)) #|(length(pch) < M & length(pch) > 1)) {
	pch <- as.character(g)
	if (length(pch) < M & length(pch) > 1) warning("pch argument should have length 1 or equal to the number of clusters")
	if(length(pch) > 1){ 
	newpch <- as.factor(g)
	levels(newpch) <- rep(levels(as.factor(pch)), length(levels(as.factor(pch))))[1:M]
	if (is.numeric(pch)) pch <- as.numeric(as.character(newpch)) else pch <- as.character(newpch)
	} 
	if (is.null(col)){#|(length(col) < M & length(col) > 1)) {
	col=g 
	if (0 %in% g) col=g+1} 
	if (length(col) < M & length(col) > 1) warning("col argument should have length 1 or equal to the number of clusters")
	if(length(col) > 1){ 
	newc <- g
	newc <- as.factor(newc)
	levels(newc) <- rep(levels(as.factor(col)), length(levels(as.factor(col))))[1:M]
	if (is.numeric(col)) col <- as.numeric(as.character(newc)) else col <- as.character(newc)
	}
	if(ncol(x)==1){
	if(is.null(list(...)$add))  plot(x[,1], pch=pch, col=col, cex=.75, ylab="x", ...)
#	if(is.null(list(...)$add))  plot(x[,1], type="n", ...)
#    for(m in 0:M) {
#      gr <- (g==m)
#      if(m==0) points(seq(along=x)[gr], x[gr,],  pch=0, cex=0.2)
#      else points(seq(along=x)[gr], x[gr,], pch=m, col=m+1, cex=0.75)
#      } 
	} else if (ncol(x)==2){
if(is.null(lab <- colnames(x))) {lab <- paste(rep("V", ncol(x)), 1:ncol(x), sep="")}	
if(is.null(list(...)$add))  plot(x, pch=pch, col=col, xlab=lab[1], ylab=lab[2], cex=.75, ...)
#	if(is.null(list(...)$add))  plot(x[,coord], type="n", ...)
#    for(m in 0:M) {
#      gr <- (g==m)
#      if(m==0) {pch=pch[1]; points(x[gr,coord],  pch=pch , cex=0.2)}
#      else points(x[gr,coord], pch=pch[m], col=m, cex=0.75)
#      }
	  } else pairs(x, pch=pch, col=col, ...)
  }
  #browser()
  invisible()
}

dbs.cluster <- function(x, clusters, h.funct="h.norm", hmult=0.75, prior, ...){
x <- as.matrix(x)
h.fn <- get(h.funct, inherit = TRUE)
M <- max(clusters)
if(M==1) stop("dbs can be computed for 2 groups at least")
lik.fine <- matrix(NA,nrow=nrow(x),ncol=M)
if(missing(prior)) prior <- rep(1/M,M) 
for(m in 1:M){
	h=(h.fn(x=x[clusters==m,]))
	lik.fine[,m] <- kepdf(x=x[clusters==m,],h=h*hmult, eval.points=x, ...)@estimate
	}
ind <- t(t(lik.fine)*prior)
den <- apply(ind,1,sum)
tau_m <- ind/den
ordered <- t(apply(tau_m,1,order,decreasing = T))
ind.ug <- which(clusters==ordered[,1])
dbs.ind <- diag(log(tau_m[,clusters]/tau_m[,ordered[,1]]))
dbs.ind[ind.ug] <- diag(log(tau_m[ind.ug,ordered[ind.ug,1]]/ind[ind.ug,ordered[ind.ug,2]]))
new("pdfSilhouette", x=x, prior=prior, dbs=dbs.ind/max(dbs.ind), clusters=clusters, nc=M, stage=NULL, call=match.call())
}

dbs.pdfCluster <- function(x, h.funct="h.norm", hmult=0.75, prior, stage=NULL, ...){
if(is.null(stage)) stage <- length(x@stages)
if (stage == 0) {
	ind.all <- which(!is.na(x@cluster.cores))
	clusters <- x@cluster.cores[ind.all] 
	} else {
	ind.all <- which(!is.na(x@stages[[stage]]))
	clusters <- x@stages[[stage]][ind.all] 
	}
y <- as.matrix(x@x)[ind.all, ]
h.fn <- get(h.funct, inherit = TRUE)
M <- max(clusters)
if(M==1) stop("dbs can be computed for 2 groups at least")
lik.fine <- matrix(NA,nrow=nrow(y),ncol=M)
if(missing(prior)) prior <- rep(1/M,M) 
for(m in 1:M){
	h=(h.fn(x=y[clusters==m,]))
	lik.fine[, m] <- kepdf(x=y[clusters==m, ],h=h*hmult, eval.points=y, ...)@estimate
	}
ind <- t(t(lik.fine)*prior)
den <- apply(ind,1,sum)
tau_m <- ind/den
ordered <- t(apply(tau_m,1,order,decreasing = T))
ind.ug <- which(clusters==ordered[,1])
dbs.ind <- diag(log(tau_m[,clusters]/tau_m[,ordered[,1]]))
dbs.ind[ind.ug] <- diag(log(tau_m[ind.ug,ordered[ind.ug,1]]/ind[ind.ug,ordered[ind.ug,2]]))
dbs.all <- rep(NA, length(x@cluster.cores))
dbs.all[ind.all] <- dbs.ind/max(dbs.ind)
cl <- rep(NA, length(x@cluster.cores))
cl[ind.all] <- clusters
new("pdfSilhouette", x=x@x, prior=prior, dbs=dbs.all, clusters=cl, nc=M, stage=stage, call=match.call())
}

plot.pdfSilhouette <- function(x, y, xlab="", ylab="", col= NULL, lwd =3, cex=.9, cex.axis=.5, main=NULL, labels=FALSE,...){
dbs=x@dbs
nc=x@nc
#par(mar=c(2,2.1,.1,2))
if (is.null(col)) cl <- paste(rep("grey",nc),round(seq(1,70,length=nc)),sep="") else cl=rep(col,nc)
if (is.null(main)) main="dbs plot" 
plot(0,axes=F, xlim=c(min(0, min(dbs, na.rm=T)),1), ylim=c(0,sum(!is.na(dbs))),type="n",ylab=ylab, xlab=xlab, main=main, ...)
ini <- 0
#lb <- 0
pos <- numeric(nc)
lab <- numeric(0)
ind <- numeric(0)
for(m in 1:nc){
	lab.gr <- which(x@clusters==m)
	ind.gr <- dbs[lab.gr]
	lab <- c(lab,lab.gr[order(ind.gr)])
	ind.gr <- c(sort(ind.gr))
	ind <- c(ind,ind.gr)
	median.gr <- median(ind.gr)
	fin <- length(lab)
	ini <- (fin-(length(lab.gr)))+1
	segments(rep(0,(length(lab.gr))), ini:fin,ind.gr,ini:fin,lwd=lwd,col=cl[m])
	text(0.65,(ini+fin)/2, paste("cluster median dbs = ", round(median.gr,2)), cex=cex,pos=4,xpd=T)
	#segments(0,ini,0,fin,lwd=4,col=cl[m])
	}
	axis(1, seq(round(min(0, min(dbs, na.rm=T)), 1), 1, by=.1), seq(round(min(0, min(dbs, na.rm=T)), 1), 1, by=.1), cex.axis=cex.axis)
	if (labels) axis(2, at=(1:sum(!is.na(dbs))), lab, cex.axis=cex.axis)
	}

summary.pdfSilhouette <- function(object, ...){
dbs.groups <- list()
for(m in 1:max(object@clusters,na.rm=T)){
ind <- which(object@clusters==m)
dbs.groups[[m]] <- object@dbs[ind]
names(dbs.groups[[m]]) <- ind
}
names(dbs.groups) <- paste("cluster", 1:length(dbs.groups), sep="")
new("summary.pdfSilhouette", obj.class=class(object)[1], dbs.groups = dbs.groups)
}

###########################################
###########################################
###########################################
###########################################
#setMethods

#new methods

setGeneric("num.con", function(x, xd, pn=0.9, profile=TRUE, profile.grid=25, correct=FALSE, profile2=FALSE, qnv=0, plot.it=FALSE, ...) standardGeneric("num.con"))
setMethod("num.con", signature(x="matrix"),  function(x, xd, pn=0.9, profile=TRUE, profile.grid=25, correct=FALSE, profile2=FALSE, qnv=0, plot.it=FALSE, ...) {
d <- dim(x)[2]
if(d==1) out <- num.con.uni(x, xd, pn=pn, profile=profile, profile.grid=profile.grid, correct=correct, profile2=profile2, qnv=qnv, plot.it=plot.it) else out <- num.con.multi(x, xd, pn=pn, profile=profile, profile.grid=profile.grid, correct=correct, profile2=profile2, qnv=qnv, plot.it=plot.it, ...)
return(out)})

h.norm <- function(x){
   x <- as.matrix(x)
   nd <- ncol(x)
   sd <- sqrt(diag(var(x)))
   sd*(4/((nd+2)*nrow(x)))^(1/(nd+4))
 }
 
# setMethod("h.norm", "matrix", function(x){
   # x <- as.matrix(x)
   # nd <- ncol(x)
   # sd <- sqrt(diag(var(x)))
   # sd*(4/((nd+2)*nrow(x)))^(1/(nd+4))
 # })

# setMethod("h.norm", "numeric", function(x){
   # sd <- sd(x)
   # sd*(4/(3*length(x)))^(1/5)
 # })
 
setGeneric("dbs", function(x, clusters, h.funct="h.norm", hmult=1, prior, ...) standardGeneric("dbs"))
setMethod("dbs", signature(x="matrix", clusters="numeric"), dbs.cluster)
setMethod("dbs", signature(x="pdfCluster",clusters="missing"), dbs.pdfCluster)

#show methods
setMethod("show","kepdf", function(object){
cat("An S4 object of class \"",class(object),"\"","\n","\n")
cat("Call: ")
print(object@call)
cat("\n")
if(ncol(object@x)==1) 
	cat("Smothing parameter: ", object@h, "\n")
	else 
	cat("Diagonal elements of the smothing matrix: ", object@h, "\n")
	cat("\n")
	cat("Density estimate at evaluation points: ", "\n")
	print(object@estimate)
	cat("\n")
 }	)

setMethod("show","summary.kepdf", function(object){
cat("An S4 object of class \"",object@obj.class,"\"","\n","\n",sep="")
cat("The highest density data point has position", object@mode, "in the sample data", "\n","\n",sep=" ")
for(i in 1:length(object@props)){
cat("Rows of ", object@props[i], "% top density data points:", object@indices[[i]],"\n","\n",sep=" ")
}		
})

#
#setMethod("show",signature("qdelaunay"), function (object) {   
#	cat("An object of class \"qdelaunay\"", "\n","\n")
#	cat("Neighbour list object:\n")
#    nnew <- nrow(object@new)
#	nr <- length(object@nb)-nnew
#	cat("Indices greater than", nr, "are vertices of the convex hull of the data",":\n")
#    print (object@nb[1:nr])
#	invisible(object)
#}) 

#setMethod("show",signature("summary.qdelaunay"), function(object){
#cat("Number of regions:",object@nregions,"\n") 
#cat("Average number of links:", mean(object@nnb),"\n") 
#cat("Link number distribution:\n")
#print(table(object@nnb, deparse.level = 0))
#})

#
setMethod("show",signature("pdfCluster"), function(object){
	cat("Clustering via nonparametric density estimation", "\n","\n")
	h  <-  object@hmult*object@h
	if(ncol(object@x)==1) 
	cat("Smothing parameter: ", h, "\n", "\n")
	else 
	cat("Diagonal elements of the smothing matrix: ", h, "\n","\n")
	cat("Density estimate at data points: ", "\n")
	print(object@estimate)
	cat("\n")
    if (class(object@tree)=="dendrogram") {cat("Groups tree:\n")
	str(object@tree)
	cat("\n")}
  	cat("Initial groupings:", "\n")
      print(object@cluster.cores)
    cat("\n")
	if(!is.null(object@clusters)){
    for(i in seq(1:(length(object@stages)-1))){
	cat("Stage", i, "groupings:","\n")
    print(object@stages[[i]])
    cat("\n")}
	cat("Final groupings:", "\n")
      print(object@clusters)
    cat("\n")
		}
   })

setMethod("show",signature("summary.pdfCluster"), function(object){
 cat("An S4 object of class \"", object@obj.class, "\"","\n","\n",sep="")
 cat("Initial groupings:")
      print(table(object@cluster.cores,exclude=NULL))
  cat("\n")
  cat("Final groupings:");
      print(table(object@clusters))
  cat("\n")
  cat("Groups tree:\n")
  str(object@tree)
  invisible(object)})

#  
setMethod("show", "pdfSilhouette", function(object){
cat("Density-based silhouette ", "\n", "\n")
cat("dbs index: ", "\n")
print(round(object@dbs,digits=4))
cat("\n")
cat("clusters: ", "\n")
print(object@clusters)
	})

setMethod("show", "summary.pdfSilhouette", function(object){
cat("An S4 object of class", object@obj.class, "\n", "\n")
cat("Density based silhouette summary of cluster", "\n")
for(m in 1: length(object@dbs.groups)){
cat(m, ":", "\n")
print(summary(object@dbs.groups[[m]]))
cat("\n")
}
cat("Density based silhouette summary of data", "\n")
print(summary(unlist(object@dbs.groups)))
})

#plot methods

#setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setMethod("plot", signature(x="kepdf", y="missing"), plot.kepdf)

#if (!isGeneric("lines")) setGeneric("lines", function(x,...) standardGeneric("lines"))
#setMethod("lines", signature("kepdf"),  lines.kepdf)

#setMethod("plot", signature(x="qdelaunay", y="missing"), plot.qdelaunay)

setMethod("plot",signature(x="pdfCluster", y="missing"), plot.pdfCluster)

setMethod("plot", signature(x="pdfSilhouette",y="missing"),  plot.pdfSilhouette)


find.nc <- function (nb.obj) 
{
    comp <- integer(length(nb.obj))
    comp <- .Call("g_components", nb.obj, as.integer(comp), PACKAGE="pdfCluster")
    answ <- list(nc = length(unique(comp)), comp.id = comp)
    answ
}

