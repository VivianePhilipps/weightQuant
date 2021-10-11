summary.bootwrq <-
function(object,...)
    {
        noms <- rownames(object)[-nrow(object)]
        tabnoms <- matrix(unlist(str_split(rownames(object)[-nrow(object)], pattern="_", n=3)),ncol=3,byrow=TRUE)

        tau <- unique(as.numeric(tabnoms[,2]))

        dots <- list(...)
        isrq <- sapply(dots,function(x) class(x) %in% c("rq","rqs"))
        m <- NULL
        if(any(isrq))
            {
                j <- which(isrq==TRUE)[1]
                m <- dots[[j]]
            }

        if(!is.null(m))
            {
                if(any(tau != m$tau)) stop("Quantiles should be the same in model m as in the bootstrap results")

                if(length(tau==1))
                    {
                        if(any(tabnoms[,3] != names(m$coef))) stop("Coefficients should be in the same order in the model m than in the bootstrap results")
                    }
                else
                    {
                        if(any(tabnoms[,3] != rep(rownames(m$coef),length(tau)))) stop("Coefficients should be in the same order in the model m than in the bootstrap results")
                    }
            }

        
 
        res0 <- NULL
        if(length(grep("calc0",noms))) #poids non recalcules
            {
                for(j in 1:length(tau))
                    {
                        x0tau <- object[grep(paste("calc0",tau[j],sep="_"),noms),,drop=FALSE]
                        if(is.null(m))
                            {
                                m0tau <- apply(x0tau,1,mean)
                            }
                        else
                            {
                                if(length(tau)>1)
                                    {
                                        m0tau <- coef(m)[,j]
                                    }
                                else
                                    {
                                        m0tau <- coef(m)
                                    }
                            }
                        s0tau <- apply(x0tau,1,sd)
                        p0tau <- 2*pnorm(abs(m0tau/s0tau),lower.tail=FALSE)

                        res0 <- rbind(res0,cbind(m0tau,s0tau,p0tau))
                    }

                colnames(res0) <- c("coef","se","p-value")
                rownames(res0) <- tabnoms[1:nrow(res0),3]

                cat(" Without computation of the weights in each bootstrap sample :\n")
                cat("\n")

                k <- 0
                for(j in 1:length(tau))
                    {
                        cat("Quantile regression estimates for tau =",tau[j]," :\n")

                        print(res0[k+1:(nrow(res0)/length(tau)),])
                        k <- k+nrow(res0)/length(tau)
                        cat("\n")
                    }
            }
        
        res1 <- NULL
        if(length(grep("calc1",noms))) #poids recalcules
            {
                for(j in 1:length(tau))
                    {
                        x1tau <- object[grep(paste("calc1",tau[j],sep="_"),noms),,drop=FALSE]
                        if(is.null(m))
                            {
                                m1tau <- apply(x1tau,1,mean)
                            }
                        else
                            {
                                if(length(tau)>1)
                                    {
                                        m1tau <- coef(m)[,j]
                                    }
                                else
                                    {
                                        m1tau <- coef(m)
                                    }
                            }
                        s1tau <- apply(x1tau,1,sd)
                        p1tau <- 2*pnorm(abs(m1tau/s1tau),lower.tail=FALSE)

                        res1 <- rbind(res1,cbind(m1tau,s1tau,p1tau))
                    }
                
                colnames(res1) <- c("coef","se","p-value")
                rownames(res1) <- tabnoms[length(tabnoms[,3])-(nrow(res1)-1):0,3]
                
                cat(" With computation of the weights in each bootstrap sample :\n")
                cat("\n")

                k <- 0
                for(j in 1:length(tau))
                    {
                        cat("Quantile regression estimates for tau =",tau[j]," :\n")

                        print(res1[k+1:(nrow(res1)/length(tau)),])
                        k <- k+nrow(res1)/length(tau)
                        cat("\n")
                    }
            }

        return(invisible(list(results0=res0,results1=res1)))
    }
