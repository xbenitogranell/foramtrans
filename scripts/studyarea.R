library(leaflet)

#https://www.icgc.cat/en/Public-Administration-and-Enterprises/Services/Online-services-Geoservices/WMS-Geoindex/WMS-Ebre-Projecte-Life-EBRO 
wms <- 'https://geoserveis.icgc.cat/icgc_lifeebro/wms/service'

#pal <- colorNumeric("RdYlBu", values(cladium_pot))

leaflet()  %>% 
  setView(lng = 0.704648, lat = 40.709175, zoom = 9) %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") %>%
  addWMSTiles(
    wms,
    layers = c("delta_svi_sx","delta_1580"),
    options = WMSTileOptions(format = "image/png", transparent = T)) %>% 
  # addRasterImage(cladium_pot, colors = pal, opacity=0.4) %>%
  # addLegend(pal = pal, values = values(cladium_pot), title = "Cladium spp. potential distribution") %>%
  # addPolygons(data=cladium_sf, weight=2,col = 'red') %>%
  # addPolygons(data=lagoons_sf, weight=1,col = 'blue') %>%
  # addPolygons(data=reedbeds_sf, weight=1,col = 'blue') %>%
  # addPolygons(data=salicornia_sf, weight=1,col = 'blue') %>%
  addMiniMap(zoomLevelOffset = -6)
