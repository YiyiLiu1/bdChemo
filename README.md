# bdChemo: Estimating dose-specific cell division and apoptosis rates from chemo-sensitivity experiments

This package implements the method for estimating dose-specific cell birth and death rates from chemo-sensitivity experiments as described in Liu and Crawford (2017).

To install the package, run:

~~~
  install.packages("devtools")
  library("devtools")
  install_github("YiyiLiu1/bdChemo")
~~~
To get started, run:
~~~
  load("bdChemo")
  vignette("example", package="bdChemo")
~~~
If you use this package, please cite the article below.
~~~
Yiyi Liu and Forrest Crawford (2017). "Estimating dose-specific cell division and apoptosis rates from chemo-sensitivity experiments", submitted.
~~~