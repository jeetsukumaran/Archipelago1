library(ggplot2)

plot1 = function(d) {
    return ggplot(d) + aes(Lineages) + geom_histogram(binwidth=1) + facet_grid(Region~.)
}

plot2 = function(d) {
    x = data.frame(row.names=c("Region", "Lineages", "Endemics"))
    for (r in unique(d$Region)) {
        x = rbind(x, data.frame(Region=r,
                                Lineages=mean(subset(d, Region==r, Lineages)),
                                Endemics=mean(subset(d, Region==r, Endemics))))
    }
    return ggplot(d) + aes(Regions) + geom_histogram(binwidth=1)
}
