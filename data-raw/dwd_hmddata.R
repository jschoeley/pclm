# Code used for downloading and creating of the test data objects in the packages

# devtools::install_github("mpascariu/MortalityLaws")
library(MortalityLaws)

country = 'SWE'
usr = "your@email.com"
passw = "your_password"

HMD_Dx <- ReadHMD(what = 'Dx', 
                  countries = country, 
                  username = usr,
                  password = passw,
                  save = FALSE)$data

HMD_Ex <- ReadHMD(what = 'Ex', 
                  countries = country, 
                  username = usr,
                  password = passw,
                  save = FALSE)$data

hmdDx <- HMD_Dx
hmdEx <- HMD_Ex

devtools::use_data(hmdDx, overwrite = TRUE)
devtools::use_data(hmdEx, overwrite = TRUE)
