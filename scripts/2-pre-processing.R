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


## Bind dataframes in long format and spread
df1 <- bind_rows(wa1df, wa2df, wa3df, wa4df, wa5df, wa6df,
  wa7df, wa8df, wa9df, wa10df, wa11df, wa12df, wa13df, wa14df, wa15df,
                 wa16df, wa17df, wa18df, wa19df, wa20df) %>% spread(variable, value_num)
  
# Plot leaflet map
df1 %>% 
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") %>%
  addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
  addLayersControl(baseGroups = c("Toner Lite", "World Imagery")) %>%
  addCircleMarkers(~long, ~lat, 
                   clusterOptions = markerClusterOptions(),
                   popup = ~paste0("Habitat: ", as.character(habitat), "<br>",
                                   "Site: ", as.character(df1$site), "<br>")) 


####
## WA21
wa21 <- all_data[[3]][["WA-21"]]
colnames(wa21) <- wa21[1,]
wa21 <- wa21[-1,]
wa21[,1] <- rownames(wa21)
colnames(wa21)[1] <- "site"
wa21[,1] <- c("Scottetal1991")
colnames(wa21)[6] <- "sample"

site_coord_env <- wa21 %>%
  as.data.frame() %>%
  select(site, Coordinates, sample) %>%
  unite(site, site, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa21)
wa21df <- wa21 %>%
  as.data.frame() %>%
  select(seq(7,23)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

## WA22
wa22 <- all_data[[3]][["WA-22"]]
colnames(wa22) <- wa22[1,]
wa22 <- wa22[-1,]
wa22[,1] <- rownames(wa22)
colnames(wa22)[1] <- "site"
wa22[,1] <- c("Phleger1965")
colnames(wa22)[6] <- "sample"
wa22[,6] <- seq(length(wa22[,6]))


site_coord_env <- wa22 %>%
  as.data.frame() %>%
  select(site, Coordinates, sample, Location) %>%
  fill(Location, .direction = "up") %>% #fill in downward values that are missing 
  fill(Location, .direction = "down") %>% #fill in downward values that are missing 
  unite(site, site, sample, Location) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa22)
wa22df <- wa22 %>%
  as.data.frame() %>%
  select(seq(8,58)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA23
wa23 <- all_data[[3]][["WA-23"]]
colnames(wa23) <- wa23[1,]
wa23 <- wa23[-1,]
wa23[,1] <- rownames(wa23)
colnames(wa23)[1] <- "site"
wa23[,1] <- c("Phleger1965")


site_coord_env <- wa23 %>%
  as.data.frame() %>%
  select(site, Coordinates, Station) %>%
  unite(site, site, Station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa23)
wa23df <- wa23 %>%
  as.data.frame() %>%
  select(seq(6,39)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA24
wa24 <- all_data[[3]][["WA-24"]]
colnames(wa24) <- wa24[1,]
wa24 <- wa24[-1,]
wa24[,1] <- rownames(wa24)
colnames(wa24)[1] <- "site"
wa24[,1] <- c("Debeneyetal2002")
colnames(wa24)[5] <- "sample" 


site_coord_env <- wa24 %>%
  as.data.frame() %>%
  select(site, Coordinates, sample) %>%
  unite(site, site, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa24)

wa24df <- wa24 %>%
  as.data.frame() %>%
  select(seq(6,29)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

##WA25
wa25 <- all_data[[3]][["WA-25"]]
colnames(wa25) <- wa25[1,]
wa25 <- wa25[-1,]
wa25[,1] <- rownames(wa25)
colnames(wa25)[1] <- "site"
wa25[,1] <- c("Scott1989")
colnames(wa25)[5] <- "station" 


site_coord_env <- wa25 %>%
  as.data.frame() %>%
  select(site, Coordinates, Location, station) %>%
  fill(Location, .direction="down") %>%
  unite(site, site, Location, station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa25)

wa25 <- wa25 %>%
  as.data.frame() %>%
  select(seq(6,22)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA26
wa26 <- all_data[[3]][["WA-26"]]
colnames(wa26) <- wa26[1,]
wa26 <- wa26[-1,]
wa26[,1] <- rownames(wa26)
colnames(wa26)[1] <- "site"
wa26[,1] <- c("Barbosaetal2005")


site_coord_env <- wa26 %>%
  as.data.frame() %>%
  select(site, Coordinates, Sample, "Elevation cm") %>%
  unite(site, site, Sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa26)

wa26 <- wa26 %>%
  as.data.frame() %>%
  select(seq(7,26)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA27
wa27 <- all_data[[3]][["WA-27"]]
colnames(wa27) <- wa27[1,]
wa27 <- wa27[-1,]
wa27[,1] <- rownames(wa27)
colnames(wa27)[1] <- "site"
wa27[,1] <- c("Phleger1967")


site_coord_env <- wa27 %>%
  as.data.frame() %>%
  select(site, Coordinates, Location, "Marsh area", Station) %>%
  fill(Location, .direction = "down") %>%
  fill("Marsh area", .direction = "down") %>%
  unite(site, site, Location, "Marsh area", Station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa27)

wa27 <- wa27 %>%
  as.data.frame() %>%
  select(seq(8,42)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA28
wa28 <- all_data[[3]][["WA-28"]]
colnames(wa28) <- wa28[1,]
wa28 <- wa28[-1,]
wa28[,1] <- rownames(wa28)
colnames(wa28)[1] <- "site"
wa28[,1] <- c("Ozarkoetal1997")
wa28 <- wa28[c(1:5),]
wa28 <- wa28[,-c(6:15)]
wa28 <- wa28[,-4]


site_coord_env <- wa28 %>%
  as.data.frame() %>%
  select(site, Coordinates, "Sample station") %>%
  unite(site, site,"Sample station") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa28)

wa28 <- wa28 %>%
  as.data.frame() %>%
  select(seq(9,13)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA29
wa29 <- all_data[[3]][["WA-29"]]
colnames(wa29) <- wa29[1,]
wa29 <- wa29[-1,]
wa29[,1] <- rownames(wa29)
colnames(wa29)[5] <- "sample"
wa29[,1] <- c("Scottetal1976")
colnames(wa29)[1] <- "site"


site_coord_env <- wa29 %>%
  as.data.frame() %>%
  select(site, Coordinates, Percent, sample) %>%
  unite(site, site, Percent, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa29)

wa29 <- wa29 %>%
  as.data.frame() %>%
  select(seq(6,37)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA30
wa30 <- all_data[[3]][["WA-30"]]
colnames(wa30) <- wa30[1,]
wa30 <- wa30[-1,]
wa30[,1] <- rownames(wa30)
wa30[,1] <- c("Phleger1965")
colnames(wa30)[1] <- "site"


site_coord_env <- wa30 %>%
  as.data.frame() %>%
  select(site, Traverse, "Habitat type", Coordinates, Station) %>%
  fill(Traverse, .direction = "down") %>% #fill in downward values that are missing 
  fill("Habitat type", .direction = "down") %>% #fill in downward values that are missing 
  unite(site, site, Traverse, "Habitat type", Station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa30)

wa30 <- wa30 %>%
  as.data.frame() %>%
  select(seq(7,41)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA31
wa31 <- all_data[[3]][["WA-31"]]
colnames(wa31) <- wa31[1,]
wa31 <- wa31[-1,]
wa31[,1] <- rownames(wa31)
wa31[,1] <- c("Phleger&Ewing1962")
colnames(wa31)[1] <- "site"
colnames(wa31)[4] <- "Coordinates"


site_coord_env <- wa31 %>%
  as.data.frame() %>%
  select(site, Location, Coordinates, Station) %>%
  fill(Location, .direction = "down") %>% #fill in downward values that are missing 
  unite(site, site, Location, Station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa31)

wa31 <- wa31 %>%
  as.data.frame() %>%
  select(seq(6,56)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA32
wa32 <- all_data[[3]][["WA-32"]]
colnames(wa32) <- wa32[1,]
wa32 <- wa32[-1,]
wa32[,1] <- rownames(wa32)
wa32[,1] <- c("Scottetal1995")
colnames(wa32)[1] <- "site"
colnames(wa32)[3] <- "Coordinates"


site_coord_env <- wa32 %>%
  as.data.frame() %>%
  select(site, Coordinates, "Station number", "Elevation in cm above mean sea level") %>%
  unite(site, site, "Station number") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa32)

wa32 <- wa32 %>%
  as.data.frame() %>%
  select(seq(7,16)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)



##WA33
wa33 <- all_data[[3]][["WA-33"]]
colnames(wa33) <- wa33[1,]
wa33 <- wa33[-1,]
wa33[,1] <- rownames(wa33)
wa33[,1] <- c("Phleger1970")
colnames(wa33)[1] <- "site"
colnames(wa33)[3] <- "station"
colnames(wa33)[5] <- "sample"


site_coord_env <- wa33 %>%
  as.data.frame() %>%
  select(site, Coordinates, station, sample) %>%
  fill(station, .direction = "down") %>% #fill in downward values that are missing 
  unite(site, site, station, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa33)

wa33 <- wa33 %>%
  as.data.frame() %>%
  select(seq(6,15)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA34
wa34 <- all_data[[3]][["WA-34"]]
colnames(wa34) <- wa34[1,]
wa34 <- wa34[-1,]
wa34[,1] <- rownames(wa34)
wa34[,1] <- c("Haywardetal1999")
wa34[,2] <- seq(length(wa34[,2]))
colnames(wa34)[1] <- "site"
colnames(wa34)[2] <- "sample"
wa34[,3] <- "Kaipara"
colnames(wa34)[3] <- "station"


site_coord_env <- wa34 %>%
  as.data.frame() %>%
  select(site, Coordinates, station, sample) %>%
  unite(site, site, station, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa34)

wa34 <- wa34 %>%
  as.data.frame() %>%
  select(seq(11,24)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA35
wa35 <- all_data[[3]][["WA-35"]]
colnames(wa35) <- wa35[1,]
wa35 <- wa35[-1,]
wa35[,1] <- rownames(wa35)
wa35[,1] <- c("Hortonetal2005")
colnames(wa35)[1] <- "site"
colnames(wa35)[5] <- "sample"


site_coord_env <- wa35 %>%
  as.data.frame() %>%
  select(site, Coordinates, Percent, sample) %>%
  fill(Percent, .direction = "down") %>% #fill in downward values that are missing 
  unite(site, site, Percent, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="marsh")

colnames(wa35)

wa35 <- wa35 %>%
  as.data.frame() %>%
  select(seq(6,45)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)

##WA36
wa36 <- all_data[[4]][["WA-36"]]
colnames(wa36) <- wa36[1,]
wa36 <- wa36[-1,]
wa36[,1] <- rownames(wa36)
wa36[,1] <- c("Korsun1999")
colnames(wa36)[1] <- "site"
colnames(wa36)[3] <- "Coordinates"
colnames(wa36)[7] <- "station"


site_coord_env <- wa36 %>%
  as.data.frame() %>%
  select(site, Coordinates,station, "Water depth m") %>%
  unite(site, site, station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa36)

wa36 <- wa36 %>%
  as.data.frame() %>%
  select(seq(8,33)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA37
wa37 <- all_data[[4]][["WA-37"]]
colnames(wa37) <- wa37[1,]
wa37 <- wa37[-1,]
wa37[,1] <- rownames(wa37)
wa37[,1] <- c("Alve&urray1999")
colnames(wa37)[1] <- "site"


site_coord_env <- wa37 %>%
  as.data.frame() %>%
  select(site, Coordinates, Area, "Sample no.", "Water depth (m)") %>%
  fill(Area, .direction = "down") %>%
  unite(site, site, Area, "Sample no.") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa37)

wa37 <- wa37 %>%
  as.data.frame() %>%
  select(seq(8,34)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA38
wa38 <- all_data[[4]][["WA-38"]]
colnames(wa38) <- wa38[1,]
wa38 <- wa38[-1,]
wa38[,1] <- rownames(wa38)
wa38[,1] <- c("VanVoorthuysen1960")
colnames(wa38)[1] <- "site"


site_coord_env <- wa38 %>%
  as.data.frame() %>%
  select(site, Coordinates, Sample) %>%
  unite(site, site, Sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa38)

wa38 <- wa38 %>%
  as.data.frame() %>%
  select(seq(5,9)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)



##WA39
wa39 <- all_data[[4]][["WA-39"]]
colnames(wa39) <- wa39[1,]
wa39 <- wa39[-1,]
wa39[,1] <- rownames(wa39)
wa39[,1] <- c("Wang1983")
colnames(wa39)[1] <- "site"
wa39[,6] <- seq(length(wa39[,6]))
colnames(wa39)[6] <- "sample"

site_coord_env <- wa39 %>%
  as.data.frame() %>%
  select(site, Coordinates, sample) %>%
  unite(site, site, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa39)

wa39 <- wa39 %>%
  as.data.frame() %>%
  select(seq(7,16)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)



##WA40
wa40 <- all_data[[4]][["WA-40"]]
colnames(wa40) <- wa40[1,]
wa40 <- wa40[-1,]
wa40[,1] <- rownames(wa40)
wa40[,1] <- c("Murray1965")
colnames(wa40)[1] <- "site"
wa40 <- wa40[-7,]

site_coord_env <- wa40 %>%
  as.data.frame() %>%
  select(site, Coordinates, Station) %>%
  unite(site, site, Station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa40)

wa40 <- wa40 %>%
  as.data.frame() %>%
  select(seq(8,77)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA41
wa41 <- all_data[[4]][["WA-41"]]
colnames(wa41) <- wa41[1,]
wa41 <- wa41[-1,]
wa41[,1] <- rownames(wa41)
wa41[,1] <- c("Murray1983")
colnames(wa41)[1] <- "site"
colnames(wa41)[5] <- "station"

site_coord_env <- wa41 %>%
  as.data.frame() %>%
  select(site, Coordinates, station) %>%
  unite(site, site, station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa41)

wa41 <- wa41 %>%
  as.data.frame() %>%
  select(seq(6,46)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)


##WA42
wa42 <- all_data[[4]][["WA-42"]]
colnames(wa42) <- wa42[1,]
wa42 <- wa42[-1,]
wa42[,1] <- rownames(wa42)
wa42[,1] <- c("Murray1968")
colnames(wa42)[1] <- "site"
colnames(wa42)[4] <- "station"

site_coord_env <- wa42 %>%
  as.data.frame() %>%
  select(site, Coordinates, station, Sample) %>%
  fill(station, .direction = "down") %>%
  unite(site, site, station, Sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa42)

wa42 <- wa42 %>%
  as.data.frame() %>%
  select(seq(6,40)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value)



##WA43
wa43 <- all_data[[4]][["WA-43"]]
colnames(wa43) <- wa43[1,]
wa43 <- wa43[-1,]
wa43[,1] <- rownames(wa43)
wa43[,1] <- c("Alve&Murray1994")
colnames(wa43)[1] <- "site"
colnames(wa43)[4] <- "station"
colnames(wa43)[5] <- "sample"

site_coord_env <- wa43 %>%
  as.data.frame() %>%
  select(site, Coordinates, station, sample) %>%
  fill(station, .direction = "down") %>%
  unite(site, site, station, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa43)

wa43 <- wa43 %>%
  as.data.frame() %>%
  select(seq(6,29)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value) %>%
  replace(is.na(.), 0) #replace NA with 0


##WA44
wa44 <- all_data[[4]][["WA-44"]]
colnames(wa44) <- wa44[1,]
wa44 <- wa44[-1,]
wa44[,1] <- rownames(wa44)
wa44[,1] <- c("Murray&Alve2000")
colnames(wa44)[1] <- "site"
colnames(wa44)[4] <- "station"
colnames(wa44)[5] <- "sample"

site_coord_env <- wa44 %>%
  as.data.frame() %>%
  select(site, Coordinates, station, sample) %>%
  fill(station, .direction = "down") %>%
  unite(site, site, station, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa44)

wa44 <- wa44 %>%
  as.data.frame() %>%
  select(seq(6,31)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value) %>%
  replace(is.na(.), 0) #replace NA with 0


##WA45 --> Warsash sta.2 from Alve and Murray 2000

## WA46
wa46 <- all_data[[4]][["WA-46"]]
colnames(wa46) <- wa46[1,]
wa46 <- wa46[-1,]
wa46[,1] <- rownames(wa46)
wa46[,1] <- c("Elison1984")
colnames(wa46)[1] <- "site"
colnames(wa46)[5] <- "station"
colnames(wa46)[6] <- "sample"

site_coord_env <- wa46 %>%
  as.data.frame() %>%
  select(site, Coordinates, station, sample) %>%
  unite(site, site, station, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa46)

wa46 <- wa46 %>%
  as.data.frame() %>%
  select(seq(7,15)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value) 


## WA47
wa47 <- all_data[[4]][["WA-47"]]
colnames(wa47) <- wa47[1,]
wa47 <- wa47[-1,]
wa47[,1] <- rownames(wa47)
wa47[,1] <- c("LeCampion1970")
wa47 <- wa47[-c(2,7,11,16,21,24,25,26,31,36),]
colnames(wa47)[1] <- "site"
colnames(wa47)[4] <- "habitat"
colnames(wa47)[5] <- "station"

site_coord_env <- wa47 %>%
  as.data.frame() %>%
  select(site, Coordinates, habitat, station) %>%
  fill(habitat, .direction = "down") %>%
  unite(site, site, habitat, station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa47)

wa47 <- wa47 %>%
  as.data.frame() %>%
  select(seq(6,66)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value) 


## WA48
wa48 <- all_data[[4]][["WA-48"]]
colnames(wa48) <- wa48[1,]
wa48 <- wa48[-1,]
wa48[,1] <- rownames(wa48)
wa48[,1] <- c("Cearreta1988a")
colnames(wa48)[1] <- "site"
colnames(wa48)[4] <- "year"
colnames(wa48)[5] <- "habitat"

site_coord_env <- wa48 %>%
  as.data.frame() %>%
  select(site, Coordinates, site, year, habitat, Station) %>%
  fill(year, .direction = "down") %>%
  fill(habitat, .direction = "down") %>%
  unite(site, site, year, habitat, Station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa48)

wa48 <- wa48 %>%
  as.data.frame() %>%
  select(seq(7,53)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value) %>%
  replace(is.na(.), 0) #replace NA with 0


## WA49
wa49 <- all_data[[4]][["WA-49"]]
colnames(wa49) <- wa49[1,]
wa49 <- wa49[-1,]
wa49[,1] <- rownames(wa49)
wa49[,1] <- c("Cearreta1988b")
colnames(wa49)[1] <- "site"
colnames(wa49)[4] <- "year"
colnames(wa49)[5] <- "habitat"

site_coord_env <- wa49 %>%
  as.data.frame() %>%
  select(site, Coordinates, site, year, habitat) %>%
  fill(year, .direction = "down") %>%
  fill(habitat, .direction = "down") %>%
  unite(site, site, year, habitat) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa49)

wa49 <- wa49 %>%
  as.data.frame() %>%
  select(seq(6,58)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value) %>%
  replace(is.na(.), 0) #replace NA with 0



## WA50
wa50 <- all_data[[4]][["WA-50"]]
colnames(wa50) <- wa50[1,]
wa50 <- wa50[-1,]
wa50[,1] <- rownames(wa50)
wa50[,1] <- c("Donnicietal1997")
colnames(wa50)[1] <- "site"
colnames(wa50)[7] <- "sample"

site_coord_env <- wa50 %>%
  as.data.frame() %>%
  select(site, Coordinates, year, sample, "sample number") %>%
  unite(site, site, year, sample, "sample number") %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa50)

wa50 <- wa50 %>%
  as.data.frame() %>%
  select(seq(9,47)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value) %>%
  replace(is.na(.), 0) #replace NA with 0

##WA51 --> discarded because I couldn't find the original paper showing the geographical coordinates

## WA52
wa52 <- all_data[[5]][["WA-52"]]
colnames(wa52) <- wa52[1,]
wa50 <- wa52[-1,]
wa52[,1] <- rownames(wa52)
wa52[,1] <- c("Scottetal1977")
wa52 <- wa52[-1,]
colnames(wa52)[1] <- "site"
colnames(wa52)[5] <- "sample"


site_coord_env <- wa52 %>%
  as.data.frame() %>%
  select(site, Coordinates, sample) %>%
  unite(site, site, sample) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

colnames(wa52)

wa52 <- wa52 %>%
  as.data.frame() %>%
  select(seq(6,21)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value) %>%
  replace(is.na(.), 0) #replace NA with 0

##WA53
wa53 <- all_data[[5]][["WA-53"]]
colnames(wa53) <- wa53[1,]
wa53 <- wa53[-1,]
wa53[,1] <- rownames(wa53)
wa53[,1] <- c("Scottetal1980")
colnames(wa53)[1] <- "site"


site_coord_env <- wa53 %>%
  as.data.frame() %>%
  select(site, Coordinates, Station) %>%
  unite(site, site, Station) %>%
  separate(Coordinates, into = c("lat", "long"), sep = ";", convert = TRUE) %>%
  fill(lat, .direction = "down") %>% #fill in downward values that are missing 
  fill(long, .direction = "down") %>%
  mutate(habitat="lagoon")

site_coord_env <- site_coord_env[-nrow(site_coord_env),] #delete last row (shelf sample)
colnames(wa53)

wa53 <- wa53[-nrow(wa53),] %>%
  as.data.frame() %>%
  select(seq(6,27)) %>% # species columns
  bind_cols(site_coord_env) %>% #bind site, coordinates and environmental columns
  gather(key=variable,value=value, -site, -habitat) %>% #gather
  mutate(value_num=as.numeric(value)) %>%
  select(-value) %>%
  replace(is.na(.), 0) #replace NA with 0
