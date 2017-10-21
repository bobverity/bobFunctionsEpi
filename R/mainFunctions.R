
# -----------------------------------
#' safeRead
#'
#' Reads in text file (tab-delimited by default) and replaces bad characters with replacements. Useful when readin in model output from the Imperial transmission model.
#' \cr\cr
#' WARNING - can overwrite original file with new values if \code{overwrite==TRUE}, although this is not the default.
#'
#' @export

safeRead <- function(fileName, delim="\t", useHeader=TRUE, report=TRUE, badCharacters=c("-1.#IND","1.#INF","-999"), replacements=c(0,0,0), overwrite=FALSE) {
    
    # read in raw data. All values are read in as characters at this stage, and headers are ignored
    data <- read.delim(fileName, sep=delim, header=F, colClasses="character")
    
    # replace bad characters
    badCount <- 0
    for (j in 1:length(badCharacters)) {
        badCount <- badCount + sum(data==badCharacters[j])
        data[data==badCharacters[j]] <- replacements[j]
    }
    
    # report to console
    if (report & badCount>0) {
        cat(paste(fileName,": ",badCount," bad characters replaced\n",sep=""))
    }
    
    # now convert variables to sensible types
    if (useHeader) {
        df <- data.frame(2:nrow(data))
        for (i in 1:ncol(data)) {
            df <- cbind(df, type.convert(data[-1,i]))
        }
        df <- df[,-1]
        names(df) <- as.character(data[1,])
    } else {
        df <- data.frame(1:nrow(data))
        for (i in 1:ncol(data)) {
            df <- cbind(df, type.convert(data[,i]))
        }
        df <- df[,-1]
        names(df) <- paste("X",1:ncol(data),sep="")
    }
    
    # return or write to file
    if (overwrite) {
        write.table(df, fileName, sep=delim, row.names=FALSE, col.names=useHeader, quote=FALSE)
    } else {
        return(df)
    }
}

# -----------------------------------
#' loadShapeFile
#'
#' Loads a shape file (i.e. a SpatialPolygonsDataFrame) with a given name from the inst/extdata folder.
#'
#' @export

loadShapeFile <- function(name="DRC_admin1") {
    filePath <- system.file("extdata", name, package="bobFunctionsEpi")
    if (filePath=="") {
        stop("Unable to find shape file in inst/extdata folder")
    }
    shp <- readOGR(filePath)
    return(shp)
}

# -----------------------------------
#' merge.SpatialPolygonsDataFrame
#'
#' An ordinary merge operation sometimes changes the order of rows, which when applied to SpatialPolygonsDataFrame data can mean data rows no longer correspond to the correct polygon. This function merges data while avoiding this issue.
#'
#' @export

merge.SpatialPolygonsDataFrame <- function(shp, df) {
    
    # extract data from shapefile and add key
    shp_data <- shp@data
    shp_data$mergeKey <- 1:nrow(shp_data)
    
    # merge with df and sort based on key
    m <- merge(shp_data,df,all.x=T)
    m <- m[order(m$mergeKey),]
    m <- subset(m,select=-mergeKey)
    
    # fix row names
    row.names(m) <- row.names(shp)
    
    # make final SpatialPolygonsDataFrame object
    s <- SpatialPolygonsDataFrame(geometry(shp), m)
    return(s)
}

# -----------------------------------
#' getPolyArea
#'
#' Reads in a shape file and extracts area of every polygon.
#'
#' @export

getPolyArea <- function(shp) {
    ret <- lapply(slot(shp,"polygons"), function(x){ slot(x,"area") })
    return(unlist(ret))
}

# -----------------------------------
#' rateRatio
#'
#' Computes point estimate and upper and lower confidence intervals on a ratio of rates. Default method is to enter raw counts and time periods, but if entering rates simply set time1=1 and time2=1.
#'
#' @export

rateRatio <- function(count1, time1, count2, time2, alpha=0.05) {
    
    point <- time2/time1*count1/count2
    if (count1==0 | time1==0 | count2==0 | time2==0) {
        if (count1==0) {
            return(list(point=0,LL=NaN,UL=NaN))
        } else {
            return(list(point=NaN,LL=NaN,UL=NaN))
        }
    }
    LL <- time2/time1*count1/(count2+1)*qf(0.025,2*(count2+1),2*count1)
    UL <- time2/time1*(count1+1)/count2*qf(1-0.025,2*(count1+1),2*count2)
    
    return(list(point=point,LL=LL,UL=UL))
}

