countContingency <- function(patientData, sickSyms = unique(patientData$SICK_SYM), cohortPids = unique(patientData$PERSON_ID), continue = FALSE) {
    # sickSyms: should be a list of codes for which you want to count contingency tables for.
    # It defaults to all the codes found in patientData, but do not use it that way unless you know what you are doing, because it will take very long.
    # They need not be valid KCD codes for this to work.
    # So you can add custom risk factor/outcome codes and it will still work.
    # cohortPids: is a vector of the PERSON_IDs in the cohort.
    # This is necessary because in some study configurations, the cohort may contain patients with zero risk factors and outcomes, and thereby be absent in patientData.
    # It is also important because sometimes patientData may be contaminated with data for patients you do not want to count.
    # It defaults to all the PERSON_IDs found in patientData.
    # continue: won't work. Don't bother.
    #
    # Note to programmer: do not try to add age/sex filtering functionality in this function.
    # Instead write a separate wrapper function.
    # Let's keep the engine and the configuration separate.
    # Although the columns riskFactor.kcd and disease.kcd in the output of this function seem redundant,
    # their presence makes it easier to import into other programs including igraph.
    if (dim(patientData)[2] == 0) stop('Zero patient data.')
    if (!continue) {
        patientData <- subset(patientData, (SICK_SYM %in% sickSyms) & (PERSON_ID %in% cohortPids))
        return((function() {
             l <- length(sickSyms)
             a <- rep(0, l)
             names(a) <- sickSyms
             E <- matrix(rep(0, l^2), nrow = l, dimnames = list(sickSyms, sickSyms))
             i <- 0
             r <- (Reduce(function(l, m) {
                 nini <- setdiff(sickSyms, m$ini)
                 l$a[nini] <- l$a[nini] + 1
                 l$E[nini, m$ini] <- l$E[nini, m$ini] + 1
                 l$ED[m$fin, m$ini] <- l$ED[m$fin, m$ini] + 1
                 l$EcD[m$fin, nini] <- l$EcD[m$fin, nini] + 1
                 return(list(a = l$a, E = l$E, ED = l$ED, EcD = l$EcD))
             }, by(patientData, factor(patientData$PERSON_ID, levels = cohortPids), function(d) {
                 # <by> ignores data where the factor INDICES[] == NA, and will create a NULL object for levels that do not occur in the factor.
                 if (is.null(d)) return(list(ini = NULL, fin = NULL))
                 ini <- intersect(sickSyms, subset(d, INITIAL == TRUE)$SICK_SYM)
                 fin <- intersect(sickSyms, setdiff(subset(d, INITIAL == FALSE)$SICK_SYM, ini))
                 return(list(ini = ini, fin = fin))
             }, simplify = TRUE), list(a = a, E = E, ED = E, EcD = E)))
             riskFactor.kcd <- rep(sickSyms, rep(l, l))
             disease.kcd <- rep(sickSyms, l)
             return(data.frame(riskFactor.kcd, disease.kcd, C = rep(r$a, l), E = as.vector(r$E), ED = as.vector(r$ED), EcD = as.vector(r$EcD), row.names = paste0(rep(sickSyms, rep(l, l)), "|", sickSyms), stringsAsFactors = FALSE, check.names = FALSE, fix.empty.names = FALSE))
         }) ())
    }
# Below is the older version.
# The old version is designed such that work can be continued when terminated prematurely.
# It can be accessed by using continue = TRUE.
# It probably won't work because of variable name changes.
kRem <- unique(patientData$PERSON_ID[patientData$cohort == 1])
repeat {
    if (!file.exists("pid")) {
        kPre <- integer()
    } else kPre <- as.integer(scan("pid"))
    kRem <- setdiff(kRem, kPre)
    cat(length(kRem), "left.\n")
    if (length(kRem) == 0) {
        break
    } else if (length(kRem) < 10000 || !continue) {
        kPer <- kRem
    } else {
        kPer <- sample(kRem, 10000)
    }
    kWor <- subset(patientData, PERSON_ID %in% kPer)
    kIni <- subset(kWor, cohort == 1)
    kFol <- subset(kWor, cohort == 0)
    kLen <- length(sickSyms)
    if (length(kPre) == 0) {
        kCohort <- rep(0, kLen)
        kE <- rep(0, kLen ^ 2)
        kEcD <- kED <- kE
    } else {
        kRr <- read.csv("rr.csv")
        kCohort <- kRr[["C"]][1:kLen]
        kE <- kRr[["E"]]
        kED <- kRr[["ED"]]
        kEcD <- kRr[["EcD"]]
    }
    dim(kE) <- dim(kED) <- dim(kEcD) <- c(kLen, kLen)
    kIniMas <- rep(FALSE, kLen)
    names(kIniMas) <- sickSyms
    kFinMas <- kIniMas
    library(plyr)
    
    # Tally matrix algorithm. Because bootstrapping takes too long.
    for (i in 1:length(kPer)) {
        kFinMas[] <- FALSE
        kPerson <- kPer[i]
        if(i %in% 2^(1:20)) cat("Currently on", i, "th person\n")
        # Applying to cohort count
        kDel <- kIni$PERSON_ID == kPerson
        kIniMas <- sickSyms %in% kIni$SICK_SYM[kDel]
        kIni <- subset(kIni, !kDel)
        kCohort[!kIniMas] <- kCohort[!kIniMas] + 1
        # Applying to exposure count
        kE[!kIniMas, kIniMas] <- kE[!kIniMas, kIniMas] + 1
        # Applying to exposure disease count
        kDel <- kFol$PERSON_ID == kPerson
        kFinMas <- sickSyms %in% kFol$SICK_SYM[kDel]
        kFol <- subset(kFol, !kDel)
        kFinInc <- kFinMas & !kIniMas
        kED[kFinInc, kIniMas] <- kED[kFinInc, kIniMas] + 1
        # Applying to non-exposed disease count
        kEcD[kFinInc, !kIniMas] <- kEcD[kFinInc, !kIniMas] + 1
    }

        r = data.frame(C = rep(kCohort, kLen), E = as.vector(kE), ED = as.vector(kED), EcD = as.vector(kEcD))

    write(kPer, file = "pid", append = TRUE)

    }
    r
    
}
contingencyFromIncidences <- function(incidences, studyYear = 2004, numPatients = NULL, numSyms = NULL, syms = NULL, pids = NULL) {
    if (!is.null(numPatients) && !(is.integer(numPatients) && length(numPatients) == 1 && numPatients > 0)) {
        stop('integer(1) > 0 expected for numPatients.')
    }
    if (!is.null(numSyms) && !(is.integer(numSyms) && length(numSyms) == 1 && numSyms > 0)) {
        stop('integer(1) > 0 expected for numSyms.')
    }
    initial <- incidences$year < studyYear
    if (sum(initial) == 0) stop('Zero exposure.')
    pids2 <- unique(with(incidences, PERSON_ID[initial]))
    if (!is.null(pids)) pids2 <- intersect(pids2, pids)
    pids <- pids2; rm(pids2)
    if (!is.null(numPatients) && numPatients < length(pids)) pids <- pids[1:numPatients]
    syms2 <- names(sort(table(with(incidences, SICK_SYM[initial & PERSON_ID %in% pids])), decreasing = TRUE))
    if (!is.null(syms)) syms2 <- intersect(syms2, syms)
    syms <- syms2; rm(syms2)
    if (!is.null(numSyms) && numSyms < length(syms)) syms <- syms[1:numSyms]
    incidences <- subset.incidences(incidences, syms = syms, pids = pids)
    inciData <- with(incidences, data.frame0(PERSON_ID = PERSON_ID, SICK_SYM = SICK_SYM, INITIAL = year < studyYear))
    rm(incidences)
    time <- system.time(
    contingency <- countContingency(inciData, syms, pids)
    )
    attr(contingency, "time") <- time
    attr(contingency, "studyYear") <- studyYear
    return(contingency)
}
graphFromIncidences <- function(...) {
    gc()
    contingencies <- tryCatch(
        contingencyFromIncidences(...),
        error = function(e) {
            cat(conditionMessage(e))
            warning(conditionMessage(e))
            TRUE
        },
        finally = {
        }
    )
    if (class(contingencies) == 'logical') return(make_empty_graph())
    contingencyGraph <- graph_from_data_frame(contingencies)
    # Unfortunately, <simplify> strips edge attributes.
    contingencyGraph <- delete_edges(contingencyGraph, E(contingencyGraph)[which_loop(contingencyGraph)])
    contingencyGraph$time1 <- attr(contingencies, "time")
    contingencyGraph$studyYear <- attr(contingencies, "studyYear")
    return(contingencyGraph)
}
