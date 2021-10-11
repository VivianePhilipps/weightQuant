bootwrq <-
function(B, form, tau, data, Y, X1=NULL, X2=NULL, subject, death, time, interval.death=NULL, impute=NULL, weight=NULL, wcompute=2, seed=NULL, intermittent, file=NULL, nproc=1, MPI=FALSE)
{
    ## verif arguments
    if(missing(B)) stop("Please specify the number of bootstrap sample in argument B")
    if(missing(form)) stop("Please specify the model forula in argument form")
    if(missing(tau)) stop("Please specify quantiles in argument tau")
    if(missing(data)) stop("Please specify the dataset in argument data")
    if(missing(Y)) stop("Please specify the outcome in argument Y")
    if(missing(subject)) stop("Please specify the group variable in argument subject")
    if(missing(death)) stop("Please specify death time in argument death")
    if(missing(time)) stop("Please specify time variable in argument time")
    if(missing(intermittent)) stop("Please specify if there are intermittent missing data in argument intermittent")

    if(!is.numeric(B)) stop("B should be numeric")
    if(class(form)!="formula") stop("form should be a formula")
    if(!all(is.numeric(tau))) stop("tau should contain numeric values")
    if(!(all((tau>0) & (tau<1)))) stop("tau should contain values between 0 and 1")
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
    if(is.null(interval.death)) interval.death <- 0
    if(!is.null(interval.death)){ if(!all(is.numeric(interval.death))) stop("interval.death should only contain numeric values")}
    if(!is.null(impute)){ if(!is.numeric(impute)) stop("impute should be numeric")}
    if(!is.null(weight))
        {
            if(!is.character(weight)) stop("weight should be a character")
            if(!(weight %in% colnames(data))) stop("data should contain weight")
        }
    if(!(wcompute %in% c(0,1,2))) stop("wcompute should be 0, 1 or 2")
    if(!is.null(seed)){ if(!is.numeric(seed)) stop("seed should be numeric")}
    if(!(intermittent %in% c(TRUE,FALSE))) stop("wcompute should be TRUE or FALSE")
    if(!is.numeric(nproc)) stop("nproc should be numeric")
    if(!is.null(file)) {if(!is.character(file)) stop("file should be a character")}
    
    ## parallele
    if(nproc>1)
        {
            if(MPI==TRUE)
                {
                    cl <- makeCluster(nproc, type = "MPI") 
                }
            else
                {
                    cl <- makeCluster(nproc, type = "SOCK")
                }
            registerDoParallel(cl)
        }

    ## graines
    if(missing(seed))
        {
            seed <- floor(runif(B)*1000000)
        }


    ## calcul des poids sur data si pas deja fait
    if(wcompute!=1 & missing(weight))
        {
            if(intermittent==FALSE)
                {
                    data <- weightsMMD(data=data,Y=Y,X1=X1,
                                        X2=X2,subject=subject,death=death,
                                        time=time,
                                        interval.death=interval.death)$data
                    weight <- "weight"
                }

            if(intermittent==TRUE)
                {
                    data <- weightsIMD(data=data,Y=Y,X1=X1,
                                        X2=X2,subject=subject,death=death,
                                        time=time,impute=impute)$data
                    weight <- "weight"
                }
        }

    ## B reech
    b <- NULL
    if(nproc>1)
        {
            res <- foreach(b=1:B, .combine=cbind, .errorhandling="remove") %dopar%
            {
                .onewrq(form=form,tau=tau,data=data,Y=Y,X1=X1,
                       X2=X2,subject=subject,death=death,time=time,
                       interval.death=interval.death,impute=impute,
                       weight=weight,wcompute=wcompute,seed=seed[b],
                       intermittent=intermittent)
            }
        }
    else
        {
            res <- foreach(b=1:B, .combine=cbind, .errorhandling="remove") %do%
            {
                .onewrq(form=form,tau=tau,data=data,Y=Y,X1=X1,
                       X2=X2,subject=subject,death=death,time=time,
                       interval.death=interval.death,impute=impute,
                       weight=weight,wcompute=wcompute,seed=seed[b],
                       intermittent=intermittent)
            }
        }

    if(nproc>1)
        {
            stopCluster(cl)
        }

    if(!is.null(file))
        {
            write.table(res,file=file,sep="\t")
        }

    class(res) <- "bootwrq"
    return(invisible(res))
}
