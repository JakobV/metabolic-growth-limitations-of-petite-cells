#general settings
library(xlsx)
dodge <- position_dodge(width=0.9)
logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n)) 

# ggbiplot function----

ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, 
                     groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                     labels = NULL, labels.size = 3, alpha = 1, custom.var = F,
                     custom.var.data = NA, custom.var.palette = rainbow(10),
                     var.axes = TRUE, var.alpha=1, custom.var.linetype=1, custom.var.size=1, custom.varname.size=3,
                     circle = FALSE, circle.prob = 0.69, ellipse.names = T, ellipse.names.custom = NA,
                     varname.size = 3, varname.adjust = 1.5,select = 1:nrow(df.v), 
                     varname.abbrev = FALSE, ...)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    if(custom.var==F){
      v <- pcobj$rotation
    } else {
      v <- pcobj$rotation
      v1 <- as.matrix(custom.var.data)
    }
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
  df.v <- as.data.frame(v[, choices])
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  if(custom.var == T){
    v1 <- sweep(v1, 2, d^var.scale, FUN='*')
    df.v1 <- as.data.frame(v1[, choices])
  }
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  if(custom.var == T){
    names(df.v1) <- names(df.u)
  }
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  if(custom.var == T){
    v1.scale <- rowSums(v1^2)
    df.v1 <- r * df.v1 / sqrt(max(v1.scale))
  }
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
    if(custom.var == T){
      df.v1$varname <- rownames(v1)
    }
  } else {
    df.v$varname <- rownames(v)
    if(custom.var == T){
      df.v1$varname <- rownames(v1)
    }
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  if(custom.var == T){
    df.v1$angle <- with(df.v1, (180/pi) * atan(yvar / xvar))
    df.v1$hjust = with(df.v1, (1 - varname.adjust * sign(xvar)) / 2)
  }
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()
  
  if(var.axes) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v[select,],
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas')), 
                   alpha=var.alpha,
                   color = custom.var.palette[1])
    if(custom.var == T){
      g <- g +
        geom_segment(data = df.v1,
                     aes(x = 0, y = 0, xend = xvar, yend = yvar),
                     arrow = arrow(length = unit(1/2, 'picas')), linetype=custom.var.linetype, size=custom.var.size,
                     color = custom.var.palette[2])
    }
  }
  
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha, size = labels.size)
    } else {
      g <- g + geom_point(alpha = alpha, size = labels.size)      
    }
  }
  
  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    df.u_temp <- df.u
    for(z in unique(df.u$groups)){
      if(nrow(df.u[df.u$groups == z,]) == 2) {
        temp <- aggregate(df.u[df.u$groups == z,1:2], by=list(df.u$groups[df.u$groups == z], df.u$labels[df.u$groups == z]), FUN=mean)
        temp$xvar <- temp$xvar+(temp$xvar/10)
        temp <- temp[,c(3,4,2,1)]
        names(temp)[3:4] <- c("labels","groups")
        df.u_temp <- rbind(df.u_temp, temp)
      }
    }
    
    ell <- ddply(df.u_temp, 'groups', function(x) {
      
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    group.positions <- aggregate(ell[,1:2], by=list(ell$groups), FUN=mean)
    names(group.positions) <- c( "label", "xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
    
    if(ellipse.names == T){
      
      g <- g + geom_text(data=group.positions, aes(x = xvar, y = yvar, label=label), colour="black")
      
    }}
  
  # Label the variable axes
  if(var.axes) {
    if(is.na(select)){
      g <- g + 
        geom_text(data = df.v[,], 
                  aes(label = varname, x = xvar, y = yvar, 
                      angle = angle, hjust = hjust), alpha=var.alpha,
                  color = custom.var.palette[1], size = varname.size)
    }
    if(!is.na(select)){
      g <- g + 
        geom_text(data = df.v[select,], 
                  aes(label = varname, x = xvar, y = yvar, 
                      angle = angle, hjust = hjust), alpha=var.alpha,
                  color = custom.var.palette[1], size = varname.size)
    }
    if(custom.var == T){
      g <- g + 
        geom_text(data = df.v1, 
                  aes(label = varname, x = xvar, y = yvar,
                      angle = angle, hjust = hjust), 
                  color = custom.var.palette[2], size = custom.varname.size)
    }
  }
  # Change the name of the legend for groups
  # if(!is.null(groups)) {
  #   g <- g + scale_color_brewer(name = deparse(substitute(groups)), 
  #                               palette = 'Dark2')
  # }
  
  # TODO: Add a second set of axes
  
  return(g)
}



###----

CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}
capitalize_str <- function(charcter_string){
  sapply(charcter_string, CapStr)
}

se <- function(x) sqrt(var(x)/length(x))

stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

require(scales)
mylog_trans <- function(base=exp(1), from=0) 
{
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), 
            domain = c(base^from, Inf))
}

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#color settings
palette <- "Paired"
palette2 <- "Dark2"
wt <- "black"
rho <- "red"
aco <- "orange"
mybreakslog=seq(-3,3,length.out=501)
mycollog <- colorpanel(length(mybreakslog)-1,low="blue",mid="white",high="red")
wtnoIron <- "darkgrey"
rhonoIron <- "#FF9999"
evonoIron <- "#D2E6F1"



#legends
lb1 <- bquote('wild type' ~ 'ancestor')
lb2 <- bquote('wild type' ~ 'evolved')
lb3 <- bquote(rho^0 ~ 'ancestor')
lb4 <- bquote(rho^0 ~ 'evolved')
lb5 <- bquote(rho^0)
lb6 <- bquote('wild type')
lb6b <- bquote('wt')
lb7 <- bquote(ATP3^A307P~'[ATP3-7]')
lb8 <- bquote(ATP3^I129T~'[ATP3-8]')
lb9 <- bquote(ATP3^I304N~'[ATP3-6]')
#lb10 <- bquote('wild type' ~ '')
lb10 <- expression(atop("wild type", ))
lb10b <- "wild type "
lb10b1 <- "wild type  + Iron"
lb10b2 <- "wild type  - Iron"
lb10c <- expression(atop("wild type", 'ATP3-6'))
lb10d <- "wild type ATP3-6"
lb11 <- expression(atop(rho^0, ''))
lb11b <- bquote(rho^0 ~ '')
lb11b1 <- bquote(rho^0 ~ ' + Iron')
lb11b2 <- bquote(rho^0 ~ ' - Iron')
lb12 <- expression(atop(rho^0, 'ATP3'))
lb12b <- bquote(rho^0 ~ 'ATP3')
lb13 <- expression(atop(rho^0, 'ATP3-6' ))
lb13b <- bquote(rho^0 ~ 'ATP3-6')
lb13b1 <- bquote(rho^0 ~ 'ATP3-6 + Iron')
lb13b2 <- bquote(rho^0 ~ 'ATP3-6 - Iron')
lb14 <- expression(atop(rho^0, 'ATP3-7' ))
lb14b <- bquote(rho^0 ~ 'ATP3-7')
lb15 <- expression(atop(rho^0, 'ATP3-67' ))
# lb12 <- bquote(rho^0 ~ 'ATP3')
# lb13 <- bquote(rho^0 ~ 'ATP3-6')
# lb14 <- bquote(rho^0 ~ 'ATP3-7')
# lb15 <- bquote(rho^0 ~ 'ATP3-67')
lb16 <- bquote('')
lb17 <- bquote(' \n ATP3')
lb18 <- bquote(' \n ATP3-6')
lb19 <- bquote(' \n ATP3-7')
lb20 <- bquote(' \n ATP3-67')
lb21 <- bquote(rho^0 ~ '1')
lb22 <- bquote(rho^0 ~ '2')
lb23 <- bquote(rho^0 ~ '3')
lb24 <- bquote('wild type' ~ '1')
lb25 <- bquote('wild type' ~ '2')
lb26 <- bquote('wild type' ~ '3')
lb27 <- bquote(rho^0 ~ '\n ATP3-6')
lb27b <-expression(atop(rho^0,'ATP3-6'))
lb28 <- bquote(rho^0 ~ '\n ATP3-7')
lb29 <- bquote(rho^0 ~ '\n ATP3-8')
lb30 <- bquote(rho^0 ~ '\n ATP3-9')
lb31 <- expression(atop(wt, '+ tert-BOOH' ))
lb32 <- expression(atop('Isocitrate dehydrogenase', '(NAD)' ))
lb33 <- expression(atop('Isocitrate dehydrogenase', '(NADP)' ))
lb34 <- bquote(Delta ~ 'aco1')
lb35 <- expression(atop('100 ?M', 'Iron [II] chloride' ))
lb36 <- expression(atop(Delta~"aco1", '' ))
lb37 <- expression(atop(Delta~"aco1", 'ATP3-6' ))
lb38 <- expression(atop('wild type', '+ H'[2]*'O'[2] ))
lb39 <- bquote(Delta~'aco1 ')
lb40 <- bquote(Delta~'aco1 ATP3-6')


