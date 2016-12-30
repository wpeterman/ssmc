ssmc
============

### An R package that provides functions to conduct Monte Carlo simulations of demographic connectivity models. These functions can assess the contributions and importance of existing populations and can also assess the potential contributions of populations following habitat creation or restoration.  

To install this package, execute the following commands in R:

```
# Install 'devtools' package, if needed
if(!("devtools" %in% list.files(.libPaths()))) {
    install.packages("devtools", repo = "http://cran.rstudio.com", dep = TRUE) 
} 

library(devtools) # Loads devtools

install_github("wpeterman/ssmc") # Download package
library(ssmc) # Installs package and the other required packages needed
```
Once the package is installed, you can further explore the functions by opening the HTML 'Vignette' using the code below, or you can view the vignette directly [**here**](https://dl.dropboxusercontent.com/u/23513016/ssmc.pdf?raw=1 "Vignette")
```
vignette('ssmc')  # Opens tutorial in web browser
```
