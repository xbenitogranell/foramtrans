library(tidyverse)

## Elegant way to read multiple excel sheets per Excel file
# load names of excel files 
files <- list.files(path = "datasets/", full.names = TRUE, pattern = ".xls")

# create function (see t that transpose the dataframe for later gather)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  sapply(sheets, function(f) t(as.data.frame(readxl::read_excel(filename, sheet = f, col_names = TRUE))), 
         simplify = FALSE)
}

# execute function for all excel files in "files"
all_data <- lapply(files, read_excel_allsheets)

# WA1
wa1 <- all_data[[1]][["WA-1"]] #it includes different studies per excel sheet
colnames(wa1) <- wa1[1,]
wa1 <- wa1[-1,]
wa1 <- wa1[,-1]
wa1[,1] <- rownames(wa1)
wa1[1:13] <- c("Alve&Murray1999")
wa1[15:30] <- c("Lutze1968")
wa1[31:nrow(wa1)] <- c("Phleger1970")

# extract site coordinates and environmental variables
site_coord_env <- wa1 %>%
  as.data.frame() %>%
  select(Percent, Area, Coordinates, "Sample no.", Salinity, c("Water depth (m)", "Temperature (°C)", "Sediment <63µm (%)", "TOC (%)")) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(Area, .direction = "down") %>%
  unite(site, "Percent", "Area", "Sample no.", sep = "_", remove = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>% #fill in downward values that are missing 
  mutate(habitat="marsh")

colnames(wa1) #check where spp are in columns
wa1df <- wa1 %>%
  as.data.frame() %>%
  select(seq(6,20)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and envionmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

# WA2
wa2 <- all_data[[1]][["WA-2"]] #it includes different studies per excel sheet
colnames(wa2) <- wa2[1,]
wa2 <- wa2[-1,]
wa2 <- wa2[,-1]
wa2[,1] <- rownames(wa2)
colnames(wa2)[2] <- "Area"

colnames(wa2)
site_coord_env <- wa2 %>%
  as.data.frame() %>%
  select(Percent, Area, Coordinates, "Altitude (m OD)") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(Area, .direction = "down") %>% #fill in downward values that are missing 
  unite(site, "Percent", "Area", sep = "_", remove = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa2)
wa2df <- wa2 %>%
  as.data.frame() %>%
  select(seq(5,30)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


# WA3
wa3 <- all_data[[1]][["WA-3"]] #it includes different studies per excel sheet
colnames(wa3) <- wa3[1,]
wa3 <- wa3[-1,]
wa3 <- wa3[,-1]
wa3[,1] <- rownames(wa3)
colnames(wa3)[1] <- "Site"
colnames(wa3)[4] <- "sample"
wa3[,2] <- c("Swallow2000")

colnames(wa3)
site_coord_env <- wa3 %>%
  as.data.frame() %>%
  select(Site, c("Site 2"), Coordinates) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  unite(site, "Site 2", "Site", sep = "_", remove = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa3)
wa3df <- wa3 %>%
  as.data.frame() %>%
  select(seq(6,13)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


# WA4
wa4 <- all_data[[1]][["WA-4"]] #it includes different studies per excel sheet
colnames(wa4) <- wa4[1,]
wa4 <- wa4[-1,]
wa4 <- wa4[,-1]
wa4[,1] <- rownames(wa4)
colnames(wa4)[1] <- "Site"
wa4[,4] <- "LeCampion1970"
colnames(wa4)[4] <- "study"

site_coord_env <- wa4 %>%
  as.data.frame() %>%
  select(Site, Coordinates, study) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  #fill("Living percent", .direction = "down") %>%
  unite(site, "study", "Site", sep = "_", remove = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa4)
wa4df <- wa4 %>%
  as.data.frame() %>%
  select(seq(5,65)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

# WA5
wa5 <- all_data[[1]][["WA-5"]] #it includes different studies per excel sheet
colnames(wa5) <- wa5[1,]
wa5 <- wa5[-1,]
wa5 <- wa5[,-1]
wa5[,1] <- rownames(wa5)
wa5[,3] <- "Cearretaetal2002"
colnames(wa5)[1] <- "site"
colnames(wa5)[5] <- "sample"

site_coord_env <- wa5 %>%
  as.data.frame() %>%
  select(site, Coordinates, Percent) %>%
  unite(site, "Percent", "site") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa5)
wa5df <- wa5 %>%
  as.data.frame() %>%
  select(seq(6,38)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

# WA6
wa6 <- all_data[[1]][["WA-6"]]
colnames(wa6) <- wa6[1,]
wa6 <- wa6[-1,]
wa6 <- wa6[,-1]
wa6[,1] <- rownames(wa6)
wa6 <- wa6[-13,]
colnames(wa6)[1] <- "site"
wa6[,1] <- c("Scottetal1979")
wa6[13,1] <- c("Petruccietal1983")
wa6[14,1] <- c("Petruccietal1983")


site_coord_env <- wa6 %>%
  as.data.frame() %>%
  select(site, Coordinates, Sample) %>%
  unite(site, "site", "Sample") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(lat, .direction = "up") %>%
  fill(long, .direction = "down") %>%
  fill(long, .direction = "up") %>%
  mutate(habitat="marsh")

colnames(wa6)
wa6df <- wa6 %>%
  as.data.frame() %>%
  select(seq(6,20)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

## Bind dataframes in wide format and spread
df <- bind_rows(wa1df, wa2df, wa3df, wa4df, wa5df, wa6df) %>%
  spread(variable, value_num)

colnames(df)

# first bind ind rows, then spread, and create a list of datasets

## Leaflet map (quickly visualize distribution of samples)
library(leaflet)

df %>% 
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") %>%
  addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
  addLayersControl(baseGroups = c("Toner Lite", "World Imagery")) %>%
  addCircleMarkers(~long, ~lat, 
                   clusterOptions = markerClusterOptions(),
                   popup = ~paste0("Habitat: ", as.character(habitat), "<br>",
                            "Site: ", as.character(df$site), "<br>")) 



## WA7
wa7 <- all_data[[2]][["WA-7"]] 

colnames(wa7) <- wa7[1,]
wa7 <- wa7[-1,]
wa7 <- wa7[,-1]
wa7[,1] <- rownames(wa7)
colnames(wa7)[1] <- "site"
wa7[,1] <- c("Scott&Martini1982")

site_coord_env <- wa7 %>%
  as.data.frame() %>%
  select(site, Coordinates, Sample, "Elevation above MSL cm") %>%
  unite(site, "site", "Sample") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa7)
wa7df <- wa7 %>%
  as.data.frame() %>%
  select("Balticammina pseudomacrescens") %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


## WA8
wa8 <- all_data[[2]][["WA-8"]]
colnames(wa8) <- wa8[1,]
wa8 <- wa8[-1,]
wa8 <- wa8[,-1]
wa8[,1] <- rownames(wa8)
colnames(wa8)[1] <- "site"
colnames(wa8)[3] <- "transect"
wa8[,"Station"] <- seq(length(wa8[,"Station"]))
wa8[,1] <- c("Scottetal1981")

colnames(wa8)

site_coord_env <- wa8 %>%
  as.data.frame() %>%
  select(site, transect, Coordinates, Station, "Height above mean sea level cm") %>%
  fill(transect, .direction = "down") %>% #fill in downward values that are missing 
  unite(site, "site", "transect", "Station") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa8)
wa8df <- wa8 %>%
  as.data.frame() %>%
  select(seq(7,18)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

## WA9
wa9 <- all_data[[2]][["WA-9"]]
colnames(wa9) <- wa9[1,]
wa9 <- wa9[-1,]
wa9 <- wa9[,-1]
wa9[,1] <- rownames(wa9)
colnames(wa9)[1] <- "site"
wa9[,1] <- c("Scott&Medioli1980")

colnames(wa9)

site_coord_env <- wa9 %>%
  as.data.frame() %>%
  select(site, Coordinates, Salinity, Station, Location) %>%
  fill(Location, .direction = "down") %>%
  unite(site, site, Station, Location) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa9)
wa9df <- wa9 %>%
  as.data.frame() %>%
  select(seq(6,36)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

## WA10
wa10 <- all_data[[2]][["WA-10"]]
colnames(wa10) <- wa10[1,]
wa10 <- wa10[-1,]
wa10 <- wa10[,-1]
wa10[,1] <- rownames(wa10)
colnames(wa10)[1] <- "site"
wa10[,1] <- c("Scott&Medioli1980")

colnames(wa10)

site_coord_env <- wa10 %>%
  as.data.frame() %>%
  select(site, Coordinates, "Elevation in cm above msl", Station, Location) %>%
  fill(Location, .direction = "down") %>%
  unite(site, site, Station, Location) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa10)
wa10df <- wa10 %>%
  as.data.frame() %>%
  select(seq(6,36)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

## WA11
wa11 <- all_data[[2]][["WA-11"]]
colnames(wa11) <- wa11[1,]
wa11 <- wa11[-1,]
wa11 <- wa11[,-1]
wa11[,1] <- rownames(wa11)
colnames(wa11)[1] <- "site"
colnames(wa11)[3] <- "area"
wa11[,1] <- c("Scott&Medioli1980")

colnames(wa11)

site_coord_env <- wa11 %>%
  as.data.frame() %>%
  select(site, Coordinates, Salinity, Location, Station, area) %>%
  fill(Location, .direction = "down") %>%
  fill(area, .direction = "down") %>%
  unite(site, site, Location, area, Station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa11)
wa11df <- wa11 %>%
  as.data.frame() %>%
  select(seq(7,37)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


## WA12
wa12 <- all_data[[2]][["WA-12"]]
colnames(wa12) <- wa12[1,]
wa12 <- wa12[-1,]
wa12 <- wa12[,-1]
wa12[,1] <- rownames(wa12)
colnames(wa12)[1] <- "site"
wa12[,1] <- c("Smithetal1994")

colnames(wa12)

site_coord_env <- wa12 %>%
  as.data.frame() %>%
  select(site, Coordinates, Transect, Station) %>%
  fill(Transect, .direction = "down") %>%
  unite(site, site, Transect, Station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa12)
wa12df <- wa12 %>%
  as.data.frame() %>%
  select(seq(5,12)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


## WA13
wa13 <- all_data[[2]][["WA-13"]]
colnames(wa13) <- wa13[1,]
wa13 <- wa13[-1,]
wa13 <- wa13[,-1]
wa13[,1] <- rownames(wa13)
colnames(wa13)[1] <- "site"
wa13[,1] <- c("Pattersonetal2004")

colnames(wa13)

site_coord_env <- wa13 %>%
  as.data.frame() %>%
  select(site, Coordinates, Sample, "Elevation m") %>%
  unite(site, site, Sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa13)
wa13df <- wa13 %>%
  as.data.frame() %>%
  select(seq(7,14)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

## WA14
wa14 <- all_data[[2]][["WA-14"]]
colnames(wa14) <- wa14[1,]
wa14 <- wa14[-1,]
wa14 <- wa14[,-1]
wa14[,1] <- rownames(wa14)
colnames(wa14)[1] <- "site"
colnames(wa14)[4] <- "sample"
wa14[,4] <- seq(nrow(wa14))
wa14[,1] <- c("Gehrels1994")

site_coord_env <- wa14 %>%
  as.data.frame() %>%
  select(site, Coordinates, Percent, sample, "elevation (MHW)") %>%
  fill(Percent, .direction = "down") %>%
  unite(site, site, Percent, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa14)
wa14df <- wa14 %>%
  as.data.frame() %>%
  select(seq(8,18)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)



## WA15
wa15 <- all_data[[2]][["WA-15"]]
colnames(wa15) <- wa15[1,]
wa15 <- wa15[-1,]
wa15 <- wa15[,-1]
wa15[,1] <- rownames(wa15)
colnames(wa15)[1] <- "site"
colnames(wa15)[4] <- "sample"
wa15[7,4] <- "Low marsh"
wa15[,1] <- c("Saffert&Thomas1998")
wa15[,"Percent"] <- seq(length(wa15[,"Percent"]))

site_coord_env <- wa15 %>%
  as.data.frame() %>%
  select(site, Coordinates, "Core site", Percent, sample) %>%
  #fill("Core site", .direction="down") %>%
  fill(sample, .direction = "down") %>%
  unite(site, site, Percent, sample, "Core site") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa15)
wa15df <- wa15 %>%
  as.data.frame() %>%
  select(seq(9,16)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


## WA16
wa16 <- all_data[[2]][["WA-16"]]
colnames(wa16) <- wa16[1,]
wa16 <- wa16[-1,]
wa16 <- wa16[,-1]
wa16[,1] <- rownames(wa16)
colnames(wa16)[1] <- "site"
colnames(wa16)[5] <- "sample"
wa16[,1] <- c("Steinceck&Bergstein1979")

site_coord_env <- wa16 %>%
  as.data.frame() %>%
  select(site, Coordinates, sample) %>%
  unite(site, site, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa16)
wa16df <- wa16 %>%
  as.data.frame() %>%
  select(seq(6,16)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


## WA17
wa17 <- all_data[[2]][["WA-17"]]
colnames(wa17) <- wa17[1,]
wa17 <- wa17[-1,]
wa17 <- wa17[,-1]
wa17[,1] <- rownames(wa17)
colnames(wa17)[1] <- "site"
colnames(wa17)[5] <- "sample"
wa17[,1] <- c("Hippensteeletal2000")


site_coord_env <- wa17 %>%
  as.data.frame() %>%
  slice(c(1,11,31)) %>% #filter high, intermediate and low marsh
  select(site, Coordinates, sample, Percent) %>%
  fill(Percent, .direction="down") %>%
  unite(site, site, Percent, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa17)
wa17df <- wa17 %>%
  as.data.frame() %>%
  slice(c(1,11,31)) %>% #filter high, intermediate and low marsh
  select(seq(7,11)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


## WA18
wa18 <- all_data[[2]][["WA-18"]]
colnames(wa18) <- wa18[1,]
wa18 <- wa18[-1,]
wa18[,1] <- rownames(wa18)
colnames(wa18)[1] <- "site"
wa18[,1] <- c("Collinsetal1995")


site_coord_env <- wa18 %>%
  as.data.frame() %>%
  select(site, Coordinates, "Transect number/location", "Station number", "Elevation above MSL (cm)") %>%
  fill("Transect number/location", .direction="down") %>%
  unite(site, site, "Transect number/location", "Station number") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa18)
wa18df <- wa18 %>%
  as.data.frame() %>%
  select(seq(8,28)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

## WA19
wa19 <- all_data[[2]][["WA-19"]]
colnames(wa19) <- wa19[1,]
wa19 <- wa19[-1,]
wa19[,1] <- rownames(wa19)
colnames(wa19)[1] <- "site"
wa19[,1] <- c("Goldstein&Harben1993")
wa19 <- wa19[,-c(2,6)]
wa19 <- wa19[c(1,6,10,15),]


site_coord_env <- wa19 %>%
  as.data.frame() %>%
  select(site, Coordinates, Percent) %>%
  fill(Percent, .direction="down") %>%
  unite(site, site, Percent) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa19)
wa19df <- wa19 %>%
  as.data.frame() %>%
  select(seq(6,23)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


## WA20
wa20 <- all_data[[2]][["WA-20"]]
colnames(wa20) <- wa20[1,]
wa20 <- wa20[-1,]
wa20[,1] <- rownames(wa20)
colnames(wa20)[1] <- "site"
wa20[,1] <- c("Goldstein&Watkins1998")
colnames(wa20)[2] <- "sample"
wa20[,2] <- seq(length(wa20[,2]))

site_coord_env <- wa20 %>%
  as.data.frame() %>%
  select(site, Coordinates, sample) %>%
  unite(site, site, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa20)
wa20df <- wa20 %>%
  as.data.frame() %>%
  select(seq(8,29)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


## Bind dataframes in wide format and spread
df1 <- bind_rows(wa1df, wa2df, wa3df, wa4df, wa5df, wa6df,
  wa7df, wa8df, wa9df, wa10df, wa11df, wa12df, wa13df, wa14df, wa15df,
                 wa16df, wa17df, wa18df, wa19df, wa20df) %>% spread(variable, value_num)
  

df1 %>% 
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") %>%
  addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
  addLayersControl(baseGroups = c("Toner Lite", "World Imagery")) %>%
  addCircleMarkers(~long, ~lat, 
                   clusterOptions = markerClusterOptions(),
                   popup = ~paste0("Habitat: ", as.character(habitat), "<br>",
                                   "Site: ", as.character(df1$site), "<br>")) 


 


