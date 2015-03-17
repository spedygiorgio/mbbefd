library(mbbefd)

	
		
x <- seq(0, 1, length=101)

#negative density?
cbind(x, dmbbefd(x, 1.5, 1.75), dmbbefdR(x, 1.5, 1.75))

#negative probability?
cbind(x, pmbbefd(x, 1.5, 1.75), pmbbefdR(x, 1.5, 1.75))

#should be dirac at x=1
cbind(x, pmbbefd(x, 1/2, 1), pmbbefdR(x, 1/2, 1))

#some deficiency at x=1
cbind(x, pmbbefd(x, 1/2, 1/2), pmbbefdR(x, 1/2, 1/2))

#some deficiency at x=1
cbind(x, pmbbefd(x, -1/3, 3/2), pmbbefdR(x, -1/3, 3/2))

#still the exposure curve is correctly defined
cbind(x, mbbefdExposure(x, -1/3, 1/2), gmbbefdR(x, -1/3, 1/2))


#some deficiency at x=1
cbind(x, pmbbefd(x, -.99, 3/2), pmbbefdR(x, -.99, 3/2))

