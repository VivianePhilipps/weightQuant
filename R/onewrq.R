.onewrq <-
function(form, tau, data, Y, X1, X2, subject, death, time, interval.death, impute, weight, wcompute, seed, intermittent)
    {
        ## graine
        set.seed(seed)
        
        ## sujet dans data
        numeros <- unique(data[,subject])
        n <- length(numeros)

        ## poids dans l'echantillon de depart
        poidsechdepart <- data[,weight]
        if(wcompute!=1)
            {
                data$poidsechdepart <- data[,weight]
                data <- data[,-which(colnames(data)==weight)]
            }
        
        ## echantillon boot
        num_b <- sample(numeros, size=n, replace=TRUE)
        j_b <- sapply(num_b,function(i) which(data[,subject]==i,useNames=FALSE))
        j_b <- unlist(j_b,use.names=FALSE)
        nbmes_b <-  sapply(num_b,function(i) length(which(data[,subject]==i)),USE.NAMES=FALSE)
        ech_b <- data[j_b,]
        ech_b[,subject] <- rep(1:n,nbmes_b)
                
        ## estimation des modeles
        if(wcompute==0) ## on ne recalcule pas
            {        
                ## modeles si on ne recalcule pas les poids
                mold <- rq(formula=form,tau=tau,data=ech_b,weights=poidsechdepart)
            }
        else
            {
                if(wcompute==1) ## on recalcule
                    {
                        ## ajout des nouveaux poids
                        if(intermittent==FALSE)
                            {
                                dataw <- weightsMMD(data=ech_b,Y=Y,X1=X1,X2=X2,subject=subject,death=death,time=time,interval.death=interval.death)$data
                            }

                        if(intermittent==TRUE)
                            {
                                dataw <- weightsIMD(data=ech_b,Y=Y,X1=X1,X2=X2,subject=subject,death=death,time=time,impute=impute)$data
                            }
                        
                        ## modeles
                        mnew <- rq(formula=form,tau=tau,data=dataw,weights=weight)
                    }
                else ## on fait les 2
                    {
                        ## modeles si on ne recalcule pas les poids
                         mold <- rq(formula=form,tau=tau,data=ech_b,weights=poidsechdepart)
                 
                        ## ajout des poids
                        if(intermittent==FALSE)
                            {
                                dataw <- weightsMMD(data=ech_b,Y=Y,X1=X1,X2=X2,subject=subject,death=death,time=time,interval.death=interval.death)$data
                            }

                        if(intermittent==TRUE)
                            {
                                dataw <- weightsIMD(data=ech_b,Y=Y,X1=X1,X2=X2,subject=subject,death=death,time=time,impute=impute)$data
                            }
                      
                        ## modeles
                        mnew <- rq(formula=form,tau=tau,data=dataw,weights=weight)
                    }
            }

        

        ## garder les coef
        coef_b0 <- NULL
        nbcoef0 <- 0
        nomcoef0 <- NULL
        if(exists("mold"))
            {
                coef_b0 <- mold$coefficients
                if(length(tau)==1)
                    {
                        nbcoef0 <- length(coef_b0)
                        nomcoef0 <- paste("calc0",rep(tau,each=nbcoef0),names(coef_b0),sep="_")
                    }
                else
                    {
                        nbcoef0 <- nrow(coef_b0)
                        nomcoef0 <- paste("calc0",rep(tau,each=nbcoef0),rownames(coef_b0),sep="_")
                    }
            }

        
        coef_b1 <- NULL
        nbcoef1 <- 0
        nomcoef1 <- NULL
        if(exists("mnew"))
            {
                coef_b1 <- mnew$coefficients
                if(length(tau)==1)
                    {
                        nbcoef1 <- length(coef_b1)
                        nomcoef1 <- paste("calc1",rep(tau,each=nbcoef1),names(coef_b1),sep="_")
                    }
                else
                    {
                        nbcoef1 <- nrow(coef_b1)
                        nomcoef1 <- paste("calc1",rep(tau,each=nbcoef1),rownames(coef_b1),sep="_")
                    }
            }


        coef_b <- c(coef_b0,coef_b1)
        nomcoef <- c(nomcoef0,nomcoef1)
                
        res <- c(as.vector(coef_b),seed)
        names(res) <- c(nomcoef,"seed")
         
        return(res)
    }
