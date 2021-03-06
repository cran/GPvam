rGP.un <-                                            
function (Z_mat, fixed_effects, control) 
{
    gr.eta <- function(eta, X, Y, Z, ybetas, R_inv, G.inv) {
        gr.y <- crossprod(Z, R_inv) %*% (Y - X %*% ybetas - Z %*% 
            eta)
        gr.p.eta <- -G.inv %*% eta
        -as.vector(gr.y + gr.p.eta)
    }
    H.eta <- function(G.inv, Z, R_inv) {
        h.eta <- G.inv
        h.y <- crossprod(Z, R_inv) %*% Z
       forceSymmetric(h.eta + h.y)
    }
    ltriangle <- function(x) {
        if (!is.null(dim(x)[2])) {
            resA <- as.vector(x[lower.tri(x, diag = TRUE)])
            resA
        }
        else {
            nx <- length(x)
            d <- 0.5 * (-1 + sqrt(1 + 8 * nx))
            resB <- .symDiagonal(d)
            resB[lower.tri(resB, diag = TRUE)] <- x
            if (nx > 1) {
                resB <- resB + t(resB) - diag(diag(resB))
            }
            as(resB, "sparseMatrix")
        }
    }
    reduce.G <- function(G, nyear, nteacher, Kg) {
        if (!is.null(dim(G)[2])) {
            temp_mat <- G
            index1 <- 0
            resA <- c(NULL)
            for (j in 1:nyear) {
                temp_mat_j <- temp_mat[(index1 + 1):(index1 + 
                  Kg[j]), (index1 + 1):(index1 + Kg[j])]
                resA <- c(resA, ltriangle(as.matrix(temp_mat_j)))
                index1 <- index1 + nteacher[j] * Kg[j]
            }
            resA
        }
        else {
            resB <- Matrix(0, 0, 0,doDiag=FALSE)
            index <- 1
            for (j in 1:nyear) {
                ne <- (Kg[j] * (Kg[j] + 1))/2
                resB <- bdiag(resB, suppressMessages(kronecker(suppressMessages(.symDiagonal(nteacher[j])), 
                  ltriangle(G[index:(index + ne - 1)]))))
                index <- index + ne
            }
            rm(j)
            resB
        }
    }
    update.eta <- function(X, Y, Z, R_inv, ybetas, G, nyear, 
        cons.logLik, Ny, nstudent, n_eta) {
        G.chol <- chol(G)
        G.inv <- chol2inv(G.chol)
        H <- H.eta(G.inv, Z, R_inv)
        chol.H <- chol(H)
        var.eta <- as.matrix(chol2inv(chol.H))
        eta<- var.eta%*%as.vector(crossprod(Z, R_inv) %*% (Y - X %*% ybetas))
        log.p.eta <- -(n_eta/2) * log(2 * pi) - sum(log(diag(G.chol))) - 
            0.5 * crossprod(eta, G.inv) %*% eta
        log.p.y <- -(Ny/2) * log(2 * pi) + sum(log(diag(chol(R_inv)))) - 
            0.5 * crossprod(Y - X %*% ybetas - Z %*% eta, R_inv) %*% 
                (Y - X %*% ybetas - Z %*% eta)
        res <- var.eta
        attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + 
            log.p.y - sum(log(diag(chol.H))))
        attr(res, "eta") <- eta
        res
    }
    Score <- function(thetas, eta = eta.hat, ybetas, X, Y, Z, 
        year.count, n_ybeta, nyear, n_eta, nstudent, nteacher, 
        Kg, cons.logLik, con = control, mis.list, pattern.parmlist2, 
        pattern.count, pattern.length, pattern.Rtemplate, pattern.diag, 
        pattern.key, Ny, pattern.sum = pattern.sum) {
        n_Rparm <- nyear * (nyear + 1)/2
        G <- thetas[seq(n_Rparm + 1, length(thetas))]
        G <- reduce.G(G = G, nyear = nyear, nteacher = nteacher, 
            Kg = Kg)
        R_i <- ltriangle(as.vector(thetas[1:n_Rparm]))
        R_i.parm <- as.vector(thetas[1:n_Rparm])
        if (length(mis.list) > 0) {
            R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                R_i)[-mis.list, -mis.list]))
            R_inv <- solve(R)
        }
        else {
            R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                R_i)))
            R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                chol2inv(chol(R_i)))))
        }
        new.eta <- update.eta(X = X, Y = Y, Z = Z, 
            R_inv = R_inv, ybetas = ybetas, G = G, nyear = nyear, 
            cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, 
            n_eta = n_eta)
        eta.hat <- attr(new.eta, "eta")
        var.eta.hat <- new.eta
        temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
                        pattern.sum <- list()
for (p in unique(Z_mat$pat)) {
pattern.sum[[p]]<-R_mstep2(invsqrtW_=as.matrix(rep(1,Ny)),JYp_=as.matrix(Y.p[[p]]),loopsize_=pattern.count[[p]]/pattern.length[[p]],
  patternlength_=pattern.length[[p]],rownumber_=as.matrix(Y.p.rownumber[[p]]),ybetas_=as.matrix(ybetas),
  etahat_=as.matrix(eta.hat),tempmatR_=as.matrix(temp_mat),
  JXpi_=as.matrix(X.p[[p]]@i),JXpp_=as.matrix(X.p[[p]]@p),JXpx_=as.matrix(X.p[[p]]@x),JXpdim_=as.matrix(X.p[[p]]@Dim),
  JZpi_=as.matrix(Z.p[[p]]@i),JZpp_=as.matrix(Z.p[[p]]@p),JZpx_=as.matrix(Z.p[[p]]@x),JZpdim_=as.matrix(Z.p[[p]]@Dim))
}  
 
        score.R <- -pattern.f.score(R_i.parm, nyear, pattern.parmlist2, 
            pattern.count, pattern.length, pattern.Rtemplate, 
            pattern.diag, pattern.key, pattern.sum)
        temp_mat <- G
        gam_t <- list()
        sv_gam_t <- list()
        index1 <- 0
        for (j in 1:nyear) {
            gam_t[[j]] <- matrix(0, Kg[j], Kg[j])
            temp_mat_j <- temp_mat[(index1 + 1):(index1 + nteacher[j] * 
                Kg[j]), (index1 + 1):(index1 + nteacher[j] * 
                Kg[j])]
            gam_t[[j]] <- temp_mat_j[1:Kg[j], 1:Kg[j]]
            sv_gam_t[[j]] <- chol2inv(chol(gam_t[[j]]))
            index1 <- index1 + nteacher[j] * Kg[j]
        }
        rm(j)
        score_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
        gam_t_sc <- list()
        index1 <- 0
        score.G <- Matrix(0, 0, 0,doDiag=FALSE)
        for (j in 1:nyear) {
            gam_t_sc[[j]] <- matrix(0, Kg[j], Kg[j])
            score_mat_j <- score_mat[(index1 + 1):(index1 + nteacher[j] * 
                Kg[j]), (index1 + 1):(index1 + nteacher[j] * 
                Kg[j])]
            index2 <- c(1)
            for (k in 1:nteacher[j]) {
                gam_t_sc[[j]] <- gam_t_sc[[j]] + score_mat_j[(index2):(index2 + 
                  Kg[j] - 1), (index2):(index2 + Kg[j] - 1)]
                index2 <- index2 + Kg[j]
            }
            index1 <- index1 + nteacher[j] * Kg[j]
            der <- -0.5 * (nteacher[j] * sv_gam_t[[j]] - sv_gam_t[[j]] %*% 
                gam_t_sc[[j]] %*% sv_gam_t[[j]])
            if (is.numeric(drop(sv_gam_t[[j]]))) {
                score.eta.t <- der
            }
            else {
                score.eta.t <- 2 * der - diag(diag(der))
            }
            for (k in 1:nteacher[j]) {
                score.G <- bdiag(score.G, score.eta.t)
            }
        }
        rm(j, k)
        -c(score.R, reduce.G(G = score.G, nyear = nyear, nteacher = nteacher, 
            Kg = Kg))
    }
    update.ybeta <- function(X, Y, Z, R_inv, eta.hat) {
        A.ybeta <- crossprod(X, R_inv) %*% X
        B.ybeta <- crossprod(X, R_inv) %*% (Y - Z %*% eta.hat)
        as.vector(solve(A.ybeta, B.ybeta))
    }
    bin2dec <- function(s) sum(s * 2^(rev(seq_along(s)) - 1))
    dec2bin <- function(s) {
        L <- length(s)
        maxs <- max(s)
        digits <- floor(logb(maxs, base = 2)) + 1
        res <- array(NA, dim = c(L, digits))
        for (i in 1:digits) {
            res[, digits - i + 1] <- (s%%2)
            s <- (s%/%2)
        }
        if (L == 1) 
            res[1, ]
        else res
    }
   R_mstep2 <- function(invsqrtW_,JYp_,loopsize_, patternlength_,rownumber_,ybetas_,etahat_,tempmatR_,JXpi_,JXpp_,JXpx_,JXpdim_,JZpi_,JZpp_,JZpx_,JZpdim_){
.Call( "R_mstep_cpp",invsqrtW_,JYp_,loopsize_, patternlength_,rownumber_,ybetas_,etahat_,tempmatR_,JXpi_,JXpp_,JXpx_,JXpdim_,JZpi_,JZpp_,JZpx_,JZpdim_)
}

pattern.f.score <- function(R_i.parm, nyear, pattern.parmlist2, pattern.count, pattern.length, pattern.Rtemplate, pattern.diag, pattern.key, pattern.sum) {
     R_i <- as.matrix(ltriangle(as.vector(R_i.parm)))
    pattern.score <- numeric(nyear/2 * (nyear + 1) )
    for (p in nonempty.patterns) {
    pattern.y <- solve(pattern.f.R(R_i, p, nyear, pattern.key))
    YSY<-pattern.y %*% pattern.sum[[p]] %*% pattern.y
    PCL<-pattern.countoverlength[[p]]
        for (r.parm in 1:(nyear/2 * (nyear + 1))) {
            if(is.null(pat.coord <- pat.coord.guide[[r.parm]][[p]])) next
            pattern.yc <- pattern.y[pat.coord]
            pattern.score[r.parm] <- pattern.score[r.parm] - ( PCL* pattern.yc - YSY[pat.coord])

        }

    }
    pattern.score[1:(nyear/2 * (nyear + 1) ) %in% pattern.diag] <- 0.5 * pattern.score[1:(nyear/2 * (nyear + 1) ) %in% pattern.diag]
    -pattern.score
}
    pattern.f.R <- function(R, p, nyear, pattern.key) {
        R[pattern.key[p, ] * (1:nyear), pattern.key[p, ] * (1:nyear), 
            drop = FALSE]
    }
    Z_mat$year <- as.numeric(Z_mat$year)
    nyear <- length(unique(Z_mat$year))
    Z_mat$mis <- rep(0, dim(Z_mat)[1])
    student_list <- unique(Z_mat$student)
    Z_mat[is.na(Z_mat$y), ]$mis <- rep(1, dim(Z_mat[is.na(Z_mat$y), 
        ])[1])
    for (g in 1:nyear) {
        mis_stu <- student_list[!(student_list %in% unique(Z_mat[Z_mat$year == 
            g, ]$student))]
        le <- length(mis_stu)
        if (le > 0) {
            temp.exp <- Z_mat[1:le, ]
            temp.exp$year <- rep(g, le)
            temp.exp$mis <- rep(1, le)
            temp.exp$student <- mis_stu
            temp.exp$teacher <- rep(NA, le)
            temp.exp$y <- rep(NA, le)
            Z_mat <- rbind(Z_mat, temp.exp)
        }
    }
    rm(g, le)
    Z_mat.full <- Z_mat
    Z_mat <- Z_mat[!is.na(Z_mat$y), ]
    Z_mat.full <- Z_mat.full[which((Z_mat.full$student %in% Z_mat$student)), 
        ]
    Ny <- sum(Z_mat$mis == 0)
    nstudent <- length(unique(Z_mat$student))
    year.count <- numeric(nyear)
    for (j in 1:nyear) {
        year.count[j] <- sum(Z_mat[Z_mat$year == j, ]$mis == 
            0)
    }
    rm(j)
    RE_s_start_pos <- 1
    Kg <- c(rep(2,nyear-1),1)
    Z_mat <- Z_mat[order(Z_mat$year, Z_mat$teacher), ]
    Z_mat.full <- Z_mat.full[order(Z_mat.full$year, Z_mat.full$teacher), 
        ]
    na_list <- grep("^NA", Z_mat$teacher)
    if (length(na_list) > 0) {
        teachyearcomb <- unique(cbind(Z_mat[-na_list, ]$year, 
            Z_mat[-na_list, ]$teacher))
    }else {
        teachyearcomb <- unique(cbind(Z_mat$year, Z_mat$teacher))
    }
        nteacher <- as.vector(tapply(teachyearcomb[, 2], teachyearcomb[, 
        1], length))
    nteach_effects <- sum(nteacher*c(rep(2,nyear-1),1))
    teacheffZ_mat <- Matrix(0, nrow = nrow(Z_mat), ncol = nteach_effects,doDiag=FALSE)
    t_effects <- rep(NA, nteach_effects)
    indx <- 1
    eblup.tracker <- matrix(0, 0, 3)
    for (k in 1:nrow(teachyearcomb)) {
        student_subset <- Z_mat.full$student[Z_mat.full$year == 
            teachyearcomb[k, 1] & Z_mat.full$teacher == teachyearcomb[k, 
            2]]
        for (yr in teachyearcomb[k, 1]:nyear) {
            if (sum(is.element(Z_mat$student, student_subset) & 
                Z_mat$year == yr) != 0) {
                teacheffZ_mat[is.element(Z_mat$student, student_subset) & 
                  Z_mat$year == yr & !is.na(Z_mat$y), indx] <- 1
            }
            if (yr==as.numeric(teachyearcomb[k,1])){
             t_effects[indx] <- paste(teachyearcomb[k, 1], "_", 
                teachyearcomb[k, 2], "_current", sep = "")
                            eblup.tracker <- rbind(eblup.tracker, c(teachyearcomb[k, 
                ], yr))
                }
            if (yr==(as.numeric(teachyearcomb[k,1])+1)){
             t_effects[indx] <- paste(teachyearcomb[k, 1], "_", 
                teachyearcomb[k, 2], "_future", sep = "")
                   eblup.tracker <- rbind(eblup.tracker,c(teachyearcomb[k, 
                ], yr))
                }
            if (yr==(as.numeric(teachyearcomb[k,1]))|yr==nyear) indx <- indx + 1

        }
    }
    rm(k, yr, indx)
    Z_mat.full <- Z_mat.full[order(Z_mat.full$student, Z_mat.full$year, 
        Z_mat.full$teacher), , drop = FALSE]
    mis.list <- which(Z_mat.full$mis == 1)
    colnames(teacheffZ_mat) <- t_effects
    RE_mat <- teacheffZ_mat
    X_mat <- sparse.model.matrix(fixed_effects, Z_mat, drop.unused.levels = TRUE)
    X_mat <- X_mat[, !(colSums(abs(X_mat)) == 0), drop = FALSE]
    if (rankMatrix(X_mat,method = 'qrLINPACK')[1] != dim(X_mat)[2]) {
        stop("WARNING: Fixed-effects design matrix not full-rank")

                              
    }
    RE_mat <- RE_mat[order(Z_mat$student, Z_mat$year, Z_mat$teacher), 
        , drop = FALSE]
    X_mat <- X_mat[order(Z_mat$student, Z_mat$year, Z_mat$teacher), 
        , drop = FALSE]
    Z_mat <- Z_mat[order(Z_mat$student, Z_mat$year, Z_mat$teacher), 
        , drop = FALSE]
    Z_mat.full <- Z_mat.full[order(Z_mat.full$student, Z_mat.full$year, 
        Z_mat.full$teacher), ]
    n_eta <- nteach_effects
    n_ybeta <- dim(X_mat)[2]
    Z <- Matrix(RE_mat,doDiag=FALSE)
    huge.flag<-TRUE
    Y <- as.vector(Z_mat$y)
    X <- Matrix(X_mat,doDiag=FALSE)
    if(!huge.flag){
    Z.dense <- as.matrix(Z)
    X.dense <- as.matrix(X)
    }
    Z_mat.full$r <- 1 - Z_mat.full$mis
    pattern.student <- matrix(Z_mat.full$r, nstudent, nyear, 
        byrow = TRUE)
    Z_mat.full$pat <- rep(apply(pattern.student, 1, bin2dec), 
        each = nyear)
    Z_mat$pat <- Z_mat.full[Z_mat.full$r == 1, ]$pat
    pat <- list()
    pattern.count <- list()
    pattern.length <- list()
    X.p <- list()
    Y.p <- list()
    Z.p <- list()
        Y.p.rownumber <- list()
    pattern.countoverlength<-list()
    rownumber <- 1:Ny
    pattern.key <- dec2bin(1:(2^nyear - 1))
        X<-as(X,"sparseMatrix")
   for (p in unique(Z_mat$pat)) {
        pat[[p]] <- which(Z_mat$pat == p)
        if(!huge.flag){
        X.p[[p]] <- X.dense[pat[[p]], , drop = FALSE]     
        Z.p[[p]] <- Z.dense[pat[[p]], , drop = FALSE]
        }else{
        X.p[[p]] <- X[pat[[p]], , drop = FALSE]     
        Z.p[[p]] <- Z[pat[[p]], , drop = FALSE]
        }
        Y.p[[p]] <- Y[pat[[p]]]
        Y.p.rownumber[[p]] <- rownumber[pat[[p]]]
        pattern.count[[p]] <- length(Y.p[[p]])
        pattern.length[[p]] <- sum(pattern.key[p, ])
        pattern.countoverlength[[p]]<-pattern.count[[p]]/pattern.length[[p]]
    }
 
    pattern.yguide <- list()
    for (g in 1:nyear) {
        pattern.yguide[[g]] <- which(pattern.key[, g] == 1)
    }
    pattern.Rtemplate <- ltriangle(1:(nyear/2 * (nyear + 1)))
    pattern.diag <- diag(pattern.Rtemplate)
    pattern.Rtemplate <- ltriangle(1:(nyear/2 * (nyear + 1)))
    pattern.parmlist1 <- list()
    pattern.parmlist2 <- list()
    for (p in unique(Z_mat$pat)) {
        pattern.parmlist1[[p]] <- sort(unique(as.vector(pattern.f.R(pattern.Rtemplate, 
            p, nyear, pattern.key))))
    }
    for (r.parm in 1:(nyear/2 * (nyear + 1))) {
        pattern.parmlist2[[r.parm]] <- which(sapply(pattern.parmlist1, 
            f <- function(x) r.parm %in% x))
    }
    eta.hat <- numeric(n_eta)
    var.eta.hat <- Matrix(0, n_eta, n_eta,doDiag=FALSE)
    R.temp.comp <- numeric(nyear)
    for (g in 1:nyear) {
        R.temp.comp[g] <- var(Z_mat[Z_mat$year == g, ]$y)/4
    }
    R_i <- diag(R.temp.comp)
    if (length(mis.list) > 0) {
        R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
            R_i)[-mis.list, -mis.list]))
        R_inv <- solve(R)
    }else {
        R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
            R_i)))
        R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
            chol2inv(chol(R_i)))))
    }
    
      pat.coord.guide<-list()
for (r.parm in 1:(nyear/2 * (nyear + 1) )) {
pat.coord.guide[[r.parm]]<-list()
        for (p in pattern.parmlist2[[r.parm]]) {
            pat.coord.guide[[r.parm]][[p]] <- which(tril(pattern.f.R(pattern.Rtemplate, p, nyear, pattern.key)) == r.parm)

         }
}
    nonempty.patterns<-NULL
for(p in 1:length(pat.coord.guide[[1]])){
if(is.null(retrieve.parm<-pattern.parmlist1[[p]])) next
nonempty.patterns<-c(nonempty.patterns,p)
}
    ybetas <- update.ybeta(X, Y, Z, R_inv, eta.hat)
    names(ybetas) <- colnames(X_mat)
    G <- 100*.symDiagonal(n_eta)
    cons.logLik <- 0.5 * n_eta * log(2 * pi)
    iter <- control$max.iter.EM
    Y.mat <- Matrix(0, iter, n_ybeta,doDiag=FALSE)
    G.mat <- Matrix(0, iter, length(reduce.G(G = G, nyear = nyear, 
        nteacher = nteacher, Kg = Kg)),doDiag=FALSE)
    R.mat <- Matrix(0, iter, nyear * (nyear + 1)/2,doDiag=FALSE)
    lgLik <- numeric(iter)
    conv <- FALSE
   if (control$verbose) cat("Beginning EM algorithm\n")
    flush.console()
    for (it in 1:iter) {
        ptm <- proc.time()[3]
        suppressWarnings(rm(var.eta.hat,temp_mat))
        new.eta <- update.eta(X = X, Y = Y, Z = Z, 
            R_inv = R_inv, ybetas = ybetas, G = G, nyear = nyear, 
            cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, 
            n_eta = n_eta)
        Y.mat[it, ] <- c(ybetas)
        R.mat[it, ] <- ltriangle(as.matrix(R_i))
        G.mat[it, ] <- reduce.G(G = G, nyear = nyear, nteacher = nteacher, 
            Kg = Kg)
        eta.hat <- attr(new.eta, "eta")
        var.eta.hat <- new.eta
        temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
        lgLik[it] <- attr(new.eta, "likelihood")
        rm(new.eta)
        thets1 <- c(Y.mat[it - 1, ], R.mat[it - 1, ], G.mat[it - 
            1, ])
        thets2 <- c(Y.mat[it, ], R.mat[it, ], G.mat[it, ])
        if (it > 5) {
            check.lik <- abs(lgLik[it] - lgLik[it - 1])/abs(lgLik[it] + 
                control$tol1) < control$tol1
            if (check.lik) {
                conv <- TRUE
                if (control$verbose) {
                  cat("\n\n Algorithm converged.\n")
                  cat("\n\niter:", it, "\n")
                  cat("log-likelihood:", sprintf("%.7f", lgLik[it]), 
                    "\n")
                  cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                    lgLik[it - 1]), "\n")
                  cat("fixed effects:", round(ybetas, 4), "\n")
                  cat("R_i:\n")
                  print(round(as.matrix(R_i), 4))
                  cat("\n")
                  print(round(cov2cor(as.matrix(R_i)), 4))
                  cat("\n")
                  for (j in 1:nyear) {
                    cat("\ngamma_teach_year", j, "\n")
                    print(round(as.matrix(gam_t[[j]]), 4))
                    cat("\n")
                    print(try(round(cov2cor(as.matrix(gam_t[[j]])), 
                      4), silent = TRUE))
                    flush.console()
                  }
                  rm(j)
                }
                break
            }
            if (it==iter) {
                conv <- TRUE
                if (control$verbose) {
                  cat("\n\n Algorithm converged.\n")
                  cat("\n\niter:", it, "\n")
                  cat("log-likelihood:", sprintf("%.7f", lgLik[it]), 
                    "\n")
                  cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                    lgLik[it - 1]), "\n")
                  cat("fixed effects:", round(ybetas, 4), "\n")
                  cat("R_i:\n")
                  print(round(as.matrix(R_i), 4))
                  cat("\n")
                  print(round(cov2cor(as.matrix(R_i)), 4))
                  cat("\n")
                  for (j in 1:nyear) {
                    cat("\ngamma_teach_year", j, "\n")
                    print(round(as.matrix(gam_t[[j]]), 4))
                    cat("\n")
                    print(try(round(cov2cor(as.matrix(gam_t[[j]])), 
                      4), silent = TRUE))
                    flush.console()
                  }
                  rm(j)
                }
                break
            }
        }
        if ((control$verbose) & (it > 1)) {
            cat("\n\niter:", it, "\n")
            cat("log-likelihood:", sprintf("%.7f", lgLik[it]), 
                "\n")
            cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                lgLik[it - 1]), "\n")
            cat("fixed effects:", round(ybetas, 4), "\n")
            cat("R_i:\n")
            print(round(as.matrix(R_i), 4))
            cat("\n")
            print(round(cov2cor(as.matrix(R_i)), 4))
            cat("\n")
            for (j in 1:nyear) {
                cat("\ngamma_teach_year", j, "\n")
                print(round(as.matrix(gam_t[[j]]), 4))
                cat("\n")
                print(try(round(cov2cor(as.matrix(gam_t[[j]])), 
                  4), silent = TRUE))
                flush.console()
            }
            rm(j)
        }
        gam_t <- list()
        index1 <- 0
        for (j in 1:nyear) {
            gam_t[[j]] <- Matrix(0, Kg[j], Kg[j],doDiag=FALSE)
            temp_mat_j <- temp_mat[(index1 + 1):(index1 + nteacher[j] * 
                Kg[j]), (index1 + 1):(index1 + nteacher[j] * 
                Kg[j])]
            index2 <- c(1)
            for (k in 1:nteacher[j]) {
                gam_t[[j]] <- gam_t[[j]] + temp_mat_j[(index2):(index2 + 
                  Kg[j] - 1), (index2):(index2 + Kg[j] - 1)]
                index2 <- index2 + Kg[j]
            }
            index1 <- index1 + nteacher[j] * Kg[j]
            gam_t[[j]] <- suppressMessages(as(suppressMessages(symmpart(gam_t[[j]]/nteacher[j])), "sparseMatrix"))
        }
        rm(j, k, index2, index1)
        ybetasn <- numeric(n_ybeta)
        Gn <- Matrix(0, 0, 0,doDiag=FALSE)
        for (j in 1:nyear) {
            Gn <- bdiag(Gn, suppressMessages(kronecker(suppressMessages(.symDiagonal(nteacher[j])), 
                gam_t[[j]])))
        }
        rm(j)
        ybetasn <- update.ybeta(X, Y, Z, R_inv, eta.hat)
                   pattern.sum <- list()
for (p in unique(Z_mat$pat)) {
pattern.sum[[p]]<-R_mstep2(invsqrtW_=as.matrix(rep(1,Ny)),JYp_=as.matrix(Y.p[[p]]),loopsize_=pattern.count[[p]]/pattern.length[[p]],
  patternlength_=pattern.length[[p]],rownumber_=as.matrix(Y.p.rownumber[[p]]),ybetas_=as.matrix(ybetas),
  etahat_=as.matrix(eta.hat),tempmatR_=as.matrix(temp_mat),
  JXpi_=as.matrix(X.p[[p]]@i),JXpp_=as.matrix(X.p[[p]]@p),JXpx_=as.matrix(X.p[[p]]@x),JXpdim_=as.matrix(X.p[[p]]@Dim),
  JZpi_=as.matrix(Z.p[[p]]@i),JZpp_=as.matrix(Z.p[[p]]@p),JZpx_=as.matrix(Z.p[[p]]@x),JZpdim_=as.matrix(Z.p[[p]]@Dim))
}  
        R_i.parm <- ltriangle(as.matrix(R_i))
        R.cc <- 1
        hes.count <- 1
        R_i.parm.old <- R_i.parm
        s.prev <- numeric(length(R_i.parm))
        while (R.cc > 1e-04) {
            s <- pattern.f.score(R_i.parm = R_i.parm, nyear = nyear, 
                pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, 
                pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, 
                pattern.diag = pattern.diag, pattern.key = pattern.key, 
                pattern.sum = pattern.sum)
            j <- jacobian(pattern.f.score, c(R_i.parm), method = "simple", 
                nyear = nyear, pattern.parmlist2 = pattern.parmlist2, 
                pattern.count = pattern.count, pattern.length = pattern.length, 
                pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, 
                pattern.key = pattern.key, pattern.sum = pattern.sum)
if(it==1& (hes.count < 10)){     
           hesprod <- solve(j + max(c(diag(j), 5)) * diag(length(R_i.parm)), 
                  s)
                R.cc <- s %*% s
if(hes.count==9) R.cc<-0

          }else if ((it ==2) & (hes.count < 30)) {
                hesprod <- solve(j + max(c(diag(j), 5)*((1-hes.count/31)^2)) * diag(length(R_i.parm)), 
                  s)
                R.cc <- s %*% s
if(hes.count==29) R.cc<-0
            } else {
                hesprod <- solve(j, s)
                R.cc <- s %*% s
            }
            R_i.parm <- R_i.parm - hesprod
            hes.count <- hes.count + 1
            s.prev <- s
        }
        R_i <- ltriangle(R_i.parm)
        rm(R_i.parm)
        dimnames(R_i) <- list(NULL, NULL)
        if (length(mis.list) > 0) {
            R <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), 
                R_i))[-mis.list, -mis.list])
            R_inv <- suppressMessages(solve(R))
        }else {
            R <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), 
                R_i)))
            R_inv <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), 
                solve(R_i))))
        }
        ybetas <- ybetasn
        G <- Gn
    if (control$verbose)    cat("Iteration Time: ", proc.time()[3] - ptm, " seconds\n")
        flush.console()
    }
    names(ybetas) <- colnames(X_mat)
    thetas <- c(ltriangle(as.matrix(R_i)), reduce.G(G = G, nyear = nyear, 
        nteacher = nteacher, Kg = Kg))
    lgLik.hist <- lgLik
    lgLik <- lgLik[it]
    Hessian <- NA
     std_errors <- c(rep(NA, length(thetas)))
    if (control$hessian == TRUE) {
     if (control$verbose)   cat("Calculating Hessian of the variance components...")
        flush.console()
                 pattern.sum <- list()
for (p in unique(Z_mat$pat)) {
pattern.sum[[p]]<-R_mstep2(invsqrtW_=as.matrix(rep(1,Ny)),JYp_=as.matrix(Y.p[[p]]),loopsize_=pattern.count[[p]]/pattern.length[[p]],
  patternlength_=pattern.length[[p]],rownumber_=as.matrix(Y.p.rownumber[[p]]),ybetas_=as.matrix(ybetas),
  etahat_=as.matrix(eta.hat),tempmatR_=as.matrix(temp_mat),
  JXpi_=as.matrix(X.p[[p]]@i),JXpp_=as.matrix(X.p[[p]]@p),JXpx_=as.matrix(X.p[[p]]@x),JXpdim_=as.matrix(X.p[[p]]@Dim),
  JZpi_=as.matrix(Z.p[[p]]@i),JZpp_=as.matrix(Z.p[[p]]@p),JZpx_=as.matrix(Z.p[[p]]@x),JZpdim_=as.matrix(Z.p[[p]]@Dim))
}  
        if (control$hes.method == "richardson") {
            Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, 
                eta = eta.hat, ybetas = ybetas, X = X, Y = Y, 
                Z = Z, pattern.sum = pattern.sum, con = control, 
                year.count = year.count, n_ybeta = n_ybeta, nyear = nyear, 
                n_eta = n_eta, nstudent = nstudent, nteacher = nteacher, 
                Kg = Kg, cons.logLik = cons.logLik, mis.list = mis.list, 
                pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, 
                pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, 
                pattern.diag = pattern.diag, pattern.key = pattern.key, 
                Ny = Ny)))
        }
        else {
            Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, 
                method = "simple", eta = eta.hat, ybetas = ybetas, 
                X = X, Y = Y, Z = Z, pattern.sum = pattern.sum, 
                con = control, year.count = year.count, n_ybeta = n_ybeta, 
                nyear = nyear, n_eta = n_eta, nstudent = nstudent, 
                nteacher = nteacher, Kg = Kg, cons.logLik = cons.logLik, 
                mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, 
                pattern.count = pattern.count, pattern.length = pattern.length, 
                pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, 
                pattern.key = pattern.key, Ny = Ny)))
        }
        std_errors <- try(c(sqrt(diag(solve(Hessian)))), silent = TRUE)
        hes.warn <- FALSE
        Hessian <- round(Hessian, 4)
        if (any(eigen(Hessian)$values <= 0)) {
       if (control$verbose)     cat("Warning: Hessian not PD", "\n")
            std_errors <- c(rep(NA, length(thetas)))
            hes.warn <- TRUE
        }
    }
    c.temp <- crossprod(X, R_inv) %*% Z
    c.1 <- rbind(crossprod(X, R_inv) %*% X, t(c.temp))
    G.inv <- chol2inv(chol(G))
    c.2 <- rbind(c.temp, H.eta(G.inv, Z, R_inv))
    C_inv <- cbind(c.1, c.2)
    C <- solve(C_inv)
    eblup_stderror <- sqrt(diag(C)[-c(1:n_ybeta)])
    ybetas_stderror <- sqrt(diag(C)[1:n_ybeta])
    rm(C,C_inv,c.2,c.1,c.temp)
    eblup <- cbind(eta.hat, eblup_stderror)
    eblup <- round(eblup, 4)
    eblup <- as.data.frame(eblup)
    eblup.tracker <- as.data.frame(eblup.tracker)
    eblup <- as.data.frame(cbind(eblup.tracker, eblup))
    colnames(eblup) <- c("teacher_year", "teacher", "effect_year", 
        "EBLUP", "std_error")
    eblup$teacher <- as.character(eblup$teacher)
    t_lab <- as.vector(NULL)
    r_lab <- as.vector(NULL)
    for (j in 1:nyear) {
        ne <- (Kg[j] * (Kg[j] + 1))/2
        y <- c(NULL)
        x <- c(NULL)
        for (k in 1:Kg[j]) {
            x <- c(x, k:Kg[j])
            y <- c(y, rep(k, (Kg[j] - k + 1)))  
        }
        t_lab <- c(t_lab, paste("teacher effect from year", rep(j, 
            ne), ":[", x, ",", y, "]", sep = ""))
    }
    ne <- nyear * (nyear + 1)/2
    y <- c(NULL)
    x <- c(NULL)
    for (k in 1:nyear) {
        x <- c(x, k:nyear)
        y <- c(y, rep(k, (nyear - k + 1)))
    }
    r_lab <- paste("error covariance", ":[", x, ",", y, "]", 
        sep = "")
    rm(j, ne)
    effect_la <- c(names(ybetas), r_lab, t_lab)
    if (control$hessian == TRUE) {
        parameters <- round(cbind(c(ybetas, thetas), c(ybetas_stderror, 
            std_errors)), 4)
        colnames(parameters) <- c("Estimate", "Standard Error")
        rownames(parameters) <- as.character(effect_la)
    }
    if (control$hessian == FALSE) {
        parameters <- round(cbind(c(ybetas, thetas), c(ybetas_stderror, 
            rep(NA, length(thetas)))), 4)
        colnames(parameters) <- c("Estimate", "Standard Error")
        rownames(parameters) <- as.character(effect_la)
    }
    R_i <- round(R_i, 4)
 if (control$verbose)   cat("done.\n")
    mresid <- as.numeric(Y - X %*% ybetas)
    cresid <- as.numeric(mresid - Z %*% eta.hat)
    yhat <- as.numeric(X %*% ybetas + Z %*% eta.hat)
    yhat.m <- as.numeric(X %*% ybetas)
    Hessian <- round(Hessian, 5)
    R_i <- round(R_i, 4)
    for (i in 1:(nyear-1)) {
        gam_t[[i]] <- round(as.matrix(gam_t[[i]]), 4)
        colnames(gam_t[[i]]) <- c("current","future")
        rownames(gam_t[[i]]) <- colnames(gam_t[[i]])
    }
    i<-nyear
            gam_t[[i]] <- round(as.matrix(gam_t[[i]]), 4)
        colnames(gam_t[[i]]) <- c("current")
        rownames(gam_t[[i]]) <- colnames(gam_t[[i]])
    rchol <- try(chol(R_inv))
    yhat.s <- try(as.vector(rchol %*% (yhat)))
    sresid <- try(as.vector(rchol %*% Y - yhat.s))
    teach.cov <- lapply(gam_t, function(x) round(x, 4))
    list(loglik = lgLik, teach.effects = eblup, parameters = parameters, 
        Hessian = Hessian, R_i = as.matrix(R_i), teach.cov = gam_t, 
        mresid = mresid, cresid = cresid, y = Y, yhat = yhat, 
        stu.cov = NA, num.obs = Ny, num.student = nstudent, num.year = nyear, 
        num.teach = nteacher, yhat.m = yhat.m, sresid = sresid, 
        yhat.s = yhat.s,iter=it,persistence=control$persistence)
}
