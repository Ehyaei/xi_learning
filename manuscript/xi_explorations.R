############################################################
#                                                          #
#                 The Element of Algorithm                 #
#                                                          #
############################################################


n = 10000
u1 = rnorm(n)
u2 = rnorm(n)
u3 = rnorm(n)

x1 = u1
x2 = (x1 + u2)/sqrt(2)
x3 = (x1 + x2 + u3)/sqrt(3)

r12 = xicorln(x1, x2) %>% round(3)
r21 = xicorln(x2, x1) %>% round(3)
r13 = xicorln(x1, x3) %>% round(3)
r31 = xicorln(x3, x1) %>% round(3)
r23 = xicorln(x2, x3) %>% round(3)
r32 = xicorln(x3, x2) %>% round(3)

cat(sprintf("
X1 -> X2: %s 
X2 -> X1: %s 
X1 -> X3: %s 
X3 -> X1: %s 
X2 -> X3: %s 
X3 -> X2: %s", r12,r21, r13,r31, r23, r32))

