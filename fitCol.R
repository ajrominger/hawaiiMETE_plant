## function to color code deviations from theory, also useful for state variable plotting

fitCol <- function(x) {
    m <- 1/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    
    newx <- m * (x - min(x, na.rm = TRUE))
    
    colval <- colorRamp(hsv(c(0.17, 0.6), c(0.5, 1), c(1, 0.7)))(newx)
    colval[is.na(colval)] <- 0
    
    col <- rgb(colval, maxColorValue = 255)
    col[col == '#000000'] <- 'transparent'
    
    return(col)
}