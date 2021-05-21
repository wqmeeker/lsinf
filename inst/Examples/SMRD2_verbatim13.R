
#chapter 13  verbatim


 
#-----------------------------------------------------------

library(StatInt)
library(lsinf)

### Example 13.9
psev((log(50)-4.81)/0.50)
scaledCovarianceMatrix("sev", -1.796)
(1.96)^{2}*(6.281)/(log(2.0)^2)

(1.96)^{2}*(scaledCovarianceMatrix("sev",
             (log(50)-4.81)/0.50)["v22"])/(log(2.0)^2)



### Example 13.10
pnorm((log(120)-6.21)/1.4)
scaledCovarianceMatrix("normal", (log(120)-6.21)/1.4)
(1.96)^2*(4.685)/(log(2.0)^2)

(1.96)^2*(scaledCovarianceMatrix("normal",
        (log(120)-6.21)/1.4)["v22"])/(log(2.0)^2)



### Example 13.11
(log(50)-4.81)/0.50 
pweibull(50, shape=2, scale=exp(4.81))
varianceFactorQuantile(0.10, psev(-1.796), "sev")
(1.96)^2*(0.50)^2*(7.608)/(log(1.5)^2)

(1.96)^2*(0.50)^2*(varianceFactorQuantile(0.10,
        psev((log(50)-4.81)/0.50), "sev"))/(log(1.5)^2)



### Example 13.12
pnorm((log(120)-6.21)/1.4)
varianceFactorQuantile(0.10, 0.1548, "normal")
(1.96)^2*(1.4)^2*(2.2)/(log(2.0)^2)

(1.96)^2*(1.4)^2*(varianceFactorQuantile(0.10,
         pnorm((log(120)-6.21)/1.4), "normal"))/(log(2.0)^2)
