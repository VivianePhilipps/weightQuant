weightsMMD <-
function(data,Y,X1,X2,subject,death,time,interval.death=0,name="weight")
    {
        ## verif des arguments
        if(missing(data)) stop("Please specify the dataset in argument data")
        if(missing(Y)) stop("Please specify the outcome in argument Y")
        if(missing(subject)) stop("Please specify the group variable in argument subject")
        if(missing(death)) stop("Please specify death time in argument death")
        if(missing(time)) stop("Please specify time variable in argument time")

        if(!is.data.frame(data)) stop("data should be a data frame")
        if(!is.character(Y)) stop("Y should be a character")
        if(!(Y %in% colnames(data))) stop("data should contain Y")
        if(!is.null(X1))
            {
                if(!all(is.character(X1))) stop("X1 should only contain characters")
                if(!all((X1 %in% colnames(data)))) stop("data should contain X1")
            }        
        if(!is.null(X2))
            {
                if(!all(is.character(X2))) stop("X2 should only contain characters")
                if(!all((X2 %in% colnames(data)))) stop("data should contain X2")
            }
        if(!is.character(subject)) stop("subject should be a character")
        if(!(subject %in% colnames(data))) stop("data should contain subject")    
        if(!is.character(death)) stop("death should be a character")
        if(!(death %in% colnames(data))) stop("data should contain death")
        if(!is.character(time)) stop("time should be a character")
        if(!(time %in% colnames(data))) stop("data should contain time")
        if(!all(is.numeric(interval.death))) stop("interval.death should only contain numeric values")
        if(!is.character(name)) stop("name should be a character")
        
        ## une colonne par suivi, un sujet par ligne :
        data2 <- data[which(!is.na(data[,Y])),c(subject,time,X1,X2,Y,death)]   
        wide <- reshape(data2, v.names=Y, idvar=subject, timevar=time, direction = "wide")
        
        y.t <- paste(Y,unique(data[,time]),sep=".")
        nt <- length(y.t)

        ## indicateur obs chez les vivants
        matobs <- matrix(-1,nrow(wide),nt)
        colnames(matobs) <- paste("obs_t",1:nt,sep="")
        for(j in 1:nt)
            {
                matobs[,j] <- ifelse(!is.na(wide[,y.t[j]]),1,ifelse((is.na(wide[,death])) | (wide[,death]>j),0,NA))
            }
        
        wide_avecobs <- data.frame(wide,matobs)


        ## fonction pour un delta
        prob <- function(wide,Y,X1,X2,subject,death,suivi,delta)
            {
                matobs <- wide[,(ncol(wide)-nt+1):ncol(wide)]
                
                ## selectionner les sujets vivants pour chaque suivi
                sample <- vector("list",nt-1)
                for(j in 2:nt)
                    {
                        if((j+delta)<=nt)
                            {
                                sample[[j-1]] <- subset(wide,!is.na(matobs[,j]) & !is.na(matobs[,j+delta]))
                            }
                    }

                ## ech poole avec toutes les var
                reponse.var <- paste("obs_t",1:nt,sep="")
                nmes <- rep(NA,nt-1)
                poole <- NULL
                for(j in 1:(nt-1-delta))
                    {
                        dat <- sample[[j]][,c(subject,reponse.var[j+1],X1,X2,y.t[j])]
                        
                        dat$suivi <- j+1
                        
                        nmes[j] <- nrow(dat)
                        
                        colnames(dat) <- c(subject,"R",X1,X2,"Yavt","suivi")
                        poole <- rbind(poole,dat)
                    }
                
                colnames(poole) <- c(subject,"R",X1,X2,"Yavt","suivi")

                ## enlever les pas obs visite precedente
                poole <- poole[which(!is.na(poole$Yavt)),]

                ## regression logistique numerateur
                if(delta==0)
                    {
                        covar1 <- c(X1,X2)
                        form1 <- formula(paste("R~-1+factor(suivi)+",paste(covar1,collapse="+")))
                        reg1 <- glm(form1,family="binomial",data=poole)
                        
                        if(reg1$converged==TRUE)
                            {
                                coefnum <- reg1$coefficients
                                senum <- sqrt(diag(vcov(reg1)))
                            }
                        else
                            {
                                coefnum <- NA
                                senum <- NA
                            }
                    }

                ## regression logistique denominateur
                form2 <- "R~-1+factor(suivi)"
                if(length(X1))
                    {
                        form2 <- paste(form2,"+(",paste(X1,collapse="+"),")*Yavt",sep="")
                    }
                if(length(X2))
                    {
                        form2 <- paste(form2,"+",paste(X2,collapse="+"),sep="")
                    }
                reg2 <- glm(form2,family="binomial",data=poole)
                
                if(reg2$converged==TRUE)
                    {
                        coefden <- reg2$coefficients
                        seden <- sqrt(diag(vcov(reg2)))
                    }
                else
                    {
                        coefden <- NA
                        seden <- NA
                    }

                
                ## calcul denominateur
                dtmp <- poole[,c(subject,"R",X1,X2,"Yavt","suivi")]
                dpred <- dtmp[order(dtmp$suivi,dtmp[,subject]),]
                
                pred <- predict(reg2,newdata=dpred)
                d_avecpred <- data.frame(dpred,pden=1/(1+exp(-pred)))

                ## d_avecpred <- d_avecpred[order(d_avecpred[,subject],d_avecpred$suivi),]

                ## d_avecpred$pden <- NA
                ## for(i in 1:nrow(d_avecpred))
                ##     {
                ##         if(d_avecpred[i,"suivi"]==2)
                ##             {
                ##                 d_avecpred$pden[i] <- d_avecpred$pred[i]
                ##             }
                ##         else
                ##             {
                ##                 d_avecpred$pden[i] <- d_avecpred$pred[i]*d_avecpred$pden[i-1]
                ##             }
                ##     }
                
                ## calcul numerateur
                if(delta==0)
                    {
                        d_avecpred <- d_avecpred[order(d_avecpred[,subject],d_avecpred[,"suivi"]),]
                        pred <- predict(reg1,newdata=d_avecpred)
                        d_avecpred$pnum <- 1/(1+exp(-pred))
                        d_avecpred$pnum <- unlist(tapply(d_avecpred$pnum,d_avecpred[,subject],cumprod))
                    }

                colnames(d_avecpred)[which(colnames(d_avecpred)=="pden")] <- paste("pden",delta,sep="")
 
                if(delta==0)
                    {
                        res <- list(d_avecpred,coefden,seden,coefnum,senum)
                    }
                else
                    {
                        res <- list(d_avecpred,coefden,seden)                     
                    }
                
                return(res)
            }

        ## calculer les probas pour tous les deltas
        for(delta in interval.death)
            {
                res <- prob(wide=wide_avecobs,Y=Y,X1=X1,
                            X2=X2,subject=subject,death=death,
                            suivi=time,delta=delta)

                if(delta==0)
                    {
                        data3 <- merge(data2,res[[1]][,c(subject,"suivi","Yavt","pnum","pden0")],by.x=c(subject,time),by.y=c(subject,"suivi"),all.x=TRUE)

                        coef <- list(res[[4]],res[[2]])
                        se <- list(res[[5]],res[[3]])
                    }
                else
                    {
                        data3 <- merge(data3,res[[1]][,c(subject,"suivi",paste("pden",delta,sep=""))],by.x=c(subject,time),by.y=c(subject,"suivi"),all.x=TRUE)

                        coef <- c(coef,list(res[[2]]))
                        se <- c(se,list(res[[3]]))
                    }
            }

        ## calculer les poids
        data3 <- data3[order(data3[,subject],data3[,time]),]
  
        data3[which(data3[,time]==1),paste("pden",interval.death,sep="")] <- 1

        
        for(l in 1:nrow(data3))
            {
                if(data3[l,time]==1)
                    {
                        data3$pden[l] <- 1
                        data3$pnum[l] <- 1
                    }
                else
                    {
                        j <- data3[l,time]

                        kk <- cut(0:(j-2),breaks=c(interval.death,nt+1),labels=interval.death,right=FALSE)
                        kk <- as.numeric(as.character(kk))
                        ll <- l-0:(j-2)

                        p <- apply(cbind(ll,kk),1,function(x,d) d[x[1],paste("pden",x[2],sep="")],d=data3)
                        data3$pden[l] <- prod(p)
                        

                        # pden[l] = pden0[l] * pden1[l-1] * .. *pdendelta[l-delta]
                    }
            }

        data3$poids <- data3$pnum/data3$pden
        
        ## ajouter au data initial
        ajout <- data3[,c(subject,time,"poids")]
        colnames(ajout) <- c(subject,time,name)

        data_poids <- merge(data,ajout,all.x=TRUE,sort=FALSE)


        return(list(data=data_poids,coef=coef,se=se))
    }
