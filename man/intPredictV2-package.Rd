\name{intPredictV2-package}
\alias{intPredictV2-package}
\alias{intPredictV2}
\docType{package}
\title{
\packageTitle{intPredictV2}
}
\description{
\packageDescription{intPredictV2}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{intPredictV2}
\packageIndices{intPredictV2}
A pipeline for different gene selection and class prediction methods.
}
\author{
\packageAuthor{intPredictV2}

Maintainer: \packageMaintainer{intPredictV2}
}
\keyword{ package }
\examples{
#require(devtools)
#install_github("syspremed/intPredictV2")
data(edata) 
data(labs)
seqp=c(seq(2,44,2),47)
fac=2/3
intpred(edata,labs,seqp,fac,seed=11)
}
