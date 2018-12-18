weightsIMD <-
function(data,Y,X1,X2,subject,death,time,impute=0,name="weight")
{
    ##verif
    if(missing(data)) stop("Please specify the dataset in argument data")
    if(missing(Y)) stop("Please specify the outcome in argument Y")
    if(missing(subject)) stop("Please specify the group variable in argument subject")
    if(missing(death)) stop("Please specify death time in argument death")
    if(missing(time)) stop("Please specify time variable in argument time")
    
    if(!is.data.frame(data)) stop("data should be a data frame")
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
    if(!is.numeric(impute)) stop("impute should be numeric")
    if(!is.character(name)) stop("name should be a character")
    
    ## une colonne par suivi, un sujet par ligne :
    data2 <- data[which(!is.na(data[,Y])),c(subject,time,X1,X2,Y,death)]   
    wide <- reshape(data2, v.names=Y, idvar=subject, timevar=time, direction = "wide")

    y.t <- paste(Y,sort(unique(data[,time])),sep=".")
    nt <- length(y.t)

    ## indicateur obs chez les vivants
    matobs <- matrix(-1,nrow(wide),nt)
    colnames(matobs) <- paste("obs_t",1:nt,sep="")
    for(j in 1:nt)
        {
            matobs[,j] <- ifelse(!is.na(wide[,y.t[j]]),1,ifelse((is.na(wide[,death])) | (wide[,death]>j),0,NA))
        }

    wide_avecobs <- cbind(wide,matobs)

    ## selectionner les sujets viviants pour chaque suivi
    sample <- vector("list",nt-1)
    for(j in 2:nt)
        {
            sample[[j-1]] <- subset(wide_avecobs,!is.na(matobs[,j]))
        }

    ## ech poole avec toutes les var
    reponse.var <- paste("obs_t",1:nt,sep="")
    nmes <- rep(NA,nt-1)
    poole <- NULL
    for(j in 1:(nt-1))
        {
            dat <- sample[[j]][,c(subject,reponse.var[j],reponse.var[j+1],X1,X2,y.t[j])]
            dat[which(is.na(dat[,y.t[j]])),y.t[j]] <- impute

            dat$suivi <- j+1

            nmes[j] <- nrow(dat)

            colnames(dat) <- c(subject,"Ravt","R",X1,X2,"Yavt","suivi")
            poole <- rbind(poole,dat)
        }
    
    colnames(poole) <- c(subject,"Ravt","R",X1,X2,"Yavt","suivi")

    ## regression logistique numerateur
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
    form2 <- formula(paste(form2,"+I(1-Ravt)"))
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
    dtmp <- dtmp[order(dtmp$suivi,dtmp[,subject]),]
    dpred0 <- cbind(dtmp,Ravt=0)
    dpred1 <- cbind(dtmp,Ravt=1)
    dpred <- rbind(dpred0,dpred1)

    pred <- predict(reg2,newdata=dpred)
    np <- length(pred)/2
    d_avecpred <- cbind(dtmp,pred0=1/(1+exp(-pred[1:np])),pred1=1/(1+exp(-pred[np+1:np])))

    d_avecpred <- d_avecpred[order(d_avecpred[,subject],d_avecpred$suivi),]

    d_avecpred$pden <- NA
    for(i in 1:nrow(d_avecpred))
        {
            if(d_avecpred[i,"suivi"]==2)
                {
                    d_avecpred$pden[i] <- d_avecpred$pred1[i]
                }
            else
                {
                    d_avecpred$pden[i] <- d_avecpred$pred0[i]*(1-d_avecpred$pden[i-1])+
                        d_avecpred$pred1[i]*d_avecpred$pden[i-1]
                }
        }

    ## calcul numerateur
    pred <- predict(reg1,newdata=d_avecpred)
    d_avecpred$pnum <- 1/(1+exp(-pred))

    ## poids
    d_avecpred$poids <- d_avecpred$pnum/d_avecpred$pden
    d_poids <- d_avecpred[,c(subject,"suivi","poids")]
    colnames(d_poids) <- c(subject,time,name)

    ## ajout aux donnees initiales
    data_poids <- merge(data,d_poids,all.x=TRUE,sort=FALSE)
    data_poids[which(data_poids[,time]==1),name] <- 1

    return(list(data=data_poids,coef=list(coefnum,coefden),se=list(senum,seden)))
}
