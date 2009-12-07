library(ggplot2)

plot1 = function(d) {
    return(ggplot(d) + aes(Lineages) + geom_histogram(binwidth=1) + facet_grid(Region~.))
}

plot2 = function(d) {
    ggplot(d) + geom_point() + aes(Region, Lineages)
}
