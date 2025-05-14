
source("source.R")

####################### Input Network ##################

data("florentine_m")
G=florentine_m
plot(G)

######### Global setting ###########

M = 200 
alpha = 0.05

g.kernel= CalculateWLKernel
level = 3

########## Tesiting the fit of an ER model ###########

test_G = G
n = gorder(test_G)
C = rep(1,n)
Q_h = MLE.est(test_G, C)$estimates

p = matrix(Q_h, length(C), length(C))
diag(p) =0

GOF_IRG(test_G, C , p, M , g.kernel, level , alpha)

