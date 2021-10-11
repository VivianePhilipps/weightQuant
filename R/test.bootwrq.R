test.bootwrq <-
function(x,m)
    {
        if(missing(m)) stop("Please specify the model on initial data in argument m")
        if(length(m$tau)<2) stop("At least 2 quantile regressions should be estimated in model m")
        
        noms <- rownames(x)[-nrow(x)]
        tabnoms <- matrix(unlist(str_split(noms, pattern="_", n=3)),ncol=3,byrow=TRUE)
        
        tau <- unique(as.numeric(tabnoms[,2]))
        if(any(tau != m$tau)) stop("Quantiles should be the same in model m as in the bootstrap results")

        if(length(tau==1))
            {
                if(any(tabnoms[,3] != names(m$coef))) stop("Coefficients should be in the same order in the model m than in the bootstrap results")
            }
        else
            {
                if(any(tabnoms[,3] != rep(rownames(m$coef),length(tau)))) stop("Coefficients should be in the same order in the model m than in the bootstrap results")
            }
        
        res0 <- NULL
        if(length(grep("calc0",noms)))
            {
                cat("\n")
                cat(" Without computation of the weights in each bootstrap sample : \n \n")
          
                for(j in 1:(length(tau)-1))
                    {
                        m0tau1 <- coef(m)[,j]
                        m0tau2 <- coef(m)[,j+1]

                        x0tau1 <- x[grep(paste("calc0",tau[j],sep="_"),noms),,drop=FALSE]
                        x0tau2 <- x[grep(paste("calc0",tau[j+1],sep="_"),noms),,drop=FALSE]

                        s0diff <- apply(x0tau1-x0tau2,1,sd)
                        p0diff <- 2*pnorm(abs((m0tau1-m0tau2)/s0diff),lower.tail=FALSE)
                        cat(paste("tau =",tau[j],"versus tau =",tau[j+1],": \n"))
                        print(cbind(coef=m0tau1-m0tau2,se=s0diff,pvalue=p0diff))
                        cat("\n")
                        
                        res0 <- cbind(res0,m0tau1-m0tau2,s0diff,p0diff)
                    }
            }

        res1 <- NULL
        if(length(grep("calc1",noms)))
            {
                cat("\n")
                cat(" With computation of the weights in each bootstrap sample : \n \n")
          
                for(j in 1:(length(tau)-1))
                    {
                        m1tau1 <- coef(m)[,j]
                        m1tau2 <- coef(m)[,j+1]

                        x1tau1 <- x[grep(paste("calc1",tau[j],sep="_"),noms),,drop=FALSE]
                        x1tau2 <- x[grep(paste("calc1",tau[j+1],sep="_"),noms),,drop=FALSE]

                        s1diff <- apply(x1tau1-x1tau2,1,sd)
                        p1diff <- 2*pnorm(abs((m1tau1-m1tau2)/s1diff),lower.tail=FALSE)
                        cat(paste("tau =",tau[j],"versus tau =",tau[j+1],": \n"))
                        print(cbind(coef=m1tau1-m1tau2,se=s1diff,pvalue=p1diff))
                        cat("\n")
                        
                        res1 <- cbind(res1,m1tau1-m1tau2,s1diff,p1diff)
                    }
            }
              
        return(invisible(list(results0=res0,results1=res1)))
    }
