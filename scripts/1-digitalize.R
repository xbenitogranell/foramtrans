#Script for digitalizing tables

library(rJava)
library(tabulizer)
library(pdftools)

location <- '~/R/foramtrans/papers/Culver_Horton_2005_Infaunal marsh foraminifera from the outer banks, North Carolina, USAThe Journal of Foraminiferal Research.pdf'
location

#extract the table
out <- extract_tables(location, pages = 8)
extract_areas(location, pages = 8)




