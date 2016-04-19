# learn the data structures
# discrete genotypes vs. continous expression traits
library(qtl)
library(qpgraph)

set.seed(123567)
# simulate brassica cross
# this is for version 1.2 NOT 1.5 of the map
n.ind <- 124
genlen <- c(87, 110, 137, 61, 85, 91, 102, 75, 144, 102)
brass_markers <- c(136, 100, 212, 89, 72, 134, 139, 71, 177, 143)

BrMap <- sim.map(len = genlen, n.mar = brass_markers, eq.spacing = FALSE,
                 include.x = FALSE, anchor.tel = TRUE)
BrMap
plot(BrMap)

args(eQTLcrossParam)
eQTLcrossParam
args(reQTLcross)
eqtl <- reQTLcross(eQTLcrossParam(map = BrMap, genes = 50, type = "bc", cis = 0.5,
                     trans =rep(5,5)), a = 2, rho = 0.5)
eqtl
str(eqtl)
eqtl@model
plot(eqtl)
set.seed(23343)
cross <- sim.cross(BrMap, eqtl)
str(cross)
class(cross)

str(eqtl)
allcis <- ciseQTL(eqtl)
allcis

eQTLnetworkEstimate
showMethods("eQTLnetworkEstimate")

set.seed(1234)
g <- rgraphBAM(dRegularMarkedGraphParam(pI = 2, pY = 10, d = 3))
plot(g)
?convert2riself
?scanone

