#TRY 
rm(list = ls())
graphics.off()
#install.packages("rtry")
library(rtry)
library(tidyverse)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
#https://cran.r-project.org/web/packages/rtry/vignettes/rtry-introduction.html

#CARGAR DATOS
TRY <- rtry_import("~/Documents/MS Ecología/TFM/DATOS/36144_03102024161056 /36144.txt")
FA <- read.csv(file='/Users/vanepalmav/Documents/MS Ecología/TFM/DATOS/taxon_temp_vanesa.csv',
                 header=TRUE, sep=';') 
FA_taxon <- read.csv(file='/Users/vanepalmav/Documents/MS Ecología/TFM/DATOS/assemblages_temp_vanesa.csv',
                  header=TRUE, sep=';') 

#measurement
#rescaled_measurement
head(TRY)
unique(TRY$TraitID)
num_coincidencias <- sum(TRY$AccSpeciesName %in% FA$taxon_clean)
colnames(FA)

#Analizar por trait (ejemplo)
# Filtrar el dataframe para TraitID igual a 28
filtered_df <- TRY %>%
  filter(TraitID == 131) 
#609 no sirve

TRY_cuali <- TRYNF

# Mostrar el dataframe filtrado
print(filtered_df)

TRY_filt <- TRY[TRY$AccSpeciesName %in% FA$taxon_clean, ]
sum(TRY_filt$ErrorRisk > 4, na.rm = TRUE)

#Seleccionar solo ciertas columnas de FA solo cuando coinciden especies de TRY y FA
TRYNF <- TRY[TRY$AccSpeciesName %in% FA$taxon_clean & !is.na(TRY$TraitID),
             c("AccSpeciesName","TraitID", "DataName", "StdValue","ErrorRisk")]

#Filtrar valores NA
summary(TRYNF$StdValue) #NAs are traits that have cathegorical values (e "low")
TRYNF <- TRYNF[!is.na(TRYNF$StdValue),]

#Filtrar errores <4
TRYNF<-TRYNF[TRYNF$ErrorRisk < 5, ]

plot(TRYNF$StdValue)
summary(TRYNF$ErrorRisk)
plot(TRYNF$StdValue[TRYNF$ErrorRisk<4])

#Calculo valor medio, desv de StdValue y mean de ErrorRisk
TRYNF <- TRYNF %>%
  group_by(AccSpeciesName, TraitID) %>%  
  summarise(
    n_registros = sum(!is.na(StdValue)),                   # Número de registros no NA en StdValue
    mean_StdValue = mean(StdValue, na.rm = TRUE),          # Promedio de StdValue
    sd_StdValue = sd(StdValue, na.rm = TRUE),              # Desviación estándar de StdValue
    mean_ErrorRisk = mean(ErrorRisk, na.rm = TRUE),        # Promedio de ErrorRisk
    sd_ErrorRisk = sd(ErrorRisk, na.rm = TRUE)             # Desviación estándar de ErrorRisk
  ) %>%
  ungroup()

#Si sd es NA significa que hay solo 1 observacion

#Summary e histograma de los valores de c/trait
# Resumen por TraitID
summary_by_trait <- TRYNF %>%
  group_by(TraitID) %>%
  summarise(
    mean_of_means = mean(mean_StdValue, na.rm = TRUE),    # Media de las medias
    sd_of_means = sd(mean_StdValue, na.rm = TRUE),        # Desviación estándar de las medias
    min_mean = min(mean_StdValue, na.rm = TRUE),          # Mínimo de las medias
    max_mean = max(mean_StdValue, na.rm = TRUE),          # Máximo de las medias
    n_species = n()                                       # Número de especies
  )

# Mostrar resumen
print(summary_by_trait)

# Histograma de mean_StdValue por TraitID
ggplot(TRYNF, aes(x = mean_StdValue)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black") +  # Ajusta el binwidth según tus datos
  facet_wrap(~ TraitID, scales = "free") +   # Un histograma separado por cada TraitID
  labs(title = "Histograma de mean_StdValue por TraitID", x = "mean_StdValue", y = "Frecuencia") +
  theme_minimal()

#Tabla de especies versus frecuencias de traitname
trait_frequencies <- TRYNF_mean %>%
  # Contar combinaciones de AccSpeciesName y TraitName
  count(AccSpeciesName, TraitName, name = "Frequency") %>%
  # Convertir a formato ancho
  pivot_wider(
    names_from = TraitName,  # Las columnas serán los valores de TraitName
    values_from = Frequency, # Las celdas serán las frecuencias
    values_fill = 0          # Rellenar con 0 para combinaciones sin ocurrencias
  )

################################

#MERGE FA Y FA2
FA_merge<- merge(FA[, c("id_comm", "taxon_clean","measurement")], FA_taxon, by = "id_comm", all.x = TRUE)

# Mantener solo las columnas id_study, id_comm, taxon_clean, measurement y age
FA_clean <- FA_merge %>%
  select(id_study,id_comm, taxon_clean, age, measurement)

#ANALISIS COLONIZACION (x aparicion)

########TEST CON 1 ESTUDIO (luego automatizar)
#Filtrar por estudio
filter_id <- FA_clean %>%
  filter(id_study == "CU_Aravena et al. 2002_1 1A Saplings random plots or quadrats") 

#Filtrar por sp
filter_sp <- filter_id  %>%
  filter(taxon_clean == "Weinmannia trichosperma") 

#Ordenar age de menor a mayor
is.numeric(filter_sp$age)
filter_sp <- filter_sp %>%
  arrange(age)

# Generar todas las combinaciones de edades para cada id_study y taxon_clean
age_mix <- filter_sp %>%
  group_by(id_study, taxon_clean) %>%  
  arrange(age) %>%
  summarise(data = list(combn(age, 2, simplify = FALSE)), .groups = "drop") %>%  
  unnest(cols = c(data)) %>%  
  rename(age_mix = data) %>% 
  mutate(age_n = map_dbl(age_mix, 1),  
         age_m = map_dbl(age_mix, 2))  

# Unir con las mediciones de cada combinación
FA_comparisons <- age_mix %>%
  left_join(filter_sp, by = c("id_study", "taxon_clean", "age_n" = "age")) %>%
  rename(measurement_n = measurement) %>%
  left_join(filter_sp, by = c("id_study", "taxon_clean", "age_m" = "age")) %>%
  rename(measurement_m = measurement) %>%
  mutate(
    comparacion = case_when(
      age_n == age_m ~ NA_character_,  # Si las edades son iguales, comparacion es NA
      measurement_n != 0 & measurement_m != 0 ~ "permanece",       # Permanece
      measurement_n == 0 & measurement_m == 0 ~ "no aparece",     # No aparece
      measurement_n == 0 & measurement_m != 0 ~ "aparece",        # Aparece
      measurement_n != 0 & measurement_m == 0 ~ "desaparece",     # Desaparece
      TRUE ~ NA_character_                                       # Caso por defecto
    ),
    age_comparison = paste(age_n, age_m, sep = "_")               # Combinación de edades
  )

#Contabilizar NC y C, calcular %
#CODIGO PARA TODOS LOS ESTUDIOS Y SP

# Función para análisis y comparación entre combinaciones de edades
analisis_study_sp <- function(df) {
  # Ordenar por edad dentro de cada grupo
  df <- df %>%
    arrange(age) %>%
    # Generar todas las combinaciones posibles de edad dentro de cada grupo
    summarise(data = list(combn(age, 2, simplify = FALSE)), .groups = "drop") %>%
    unnest(cols = c(data)) %>%  # Desanidar las combinaciones
    rename(age_mix = data) %>%  # Renombrar la columna para claridad
    mutate(age_n = map_dbl(age_mix, 1),  # Extraer la primera edad de la combinación
           age_m = map_dbl(age_mix, 2))  # Extraer la segunda edad de la combinación
  

# Agrupar por id_study y taxon_clean, luego aplicar la función para cada combinación
FA_result <-
  FA_clean %>%
  group_by(id_study, taxon_clean) %>%  # Agrupar por id_study y taxon_clean
  nest() %>%  # Anidar los datos por grupo
  ungroup() %>%  # Desagrupar después de hacer el `nest()`
  ###############ERROR SOS

str(FA_result)

###################
#UNIR FA_result CON TRYNF
TRY_FA_result <- FA_result %>%
  inner_join(
    TRYNF_mean %>%
      filter(AccSpeciesName %in% FA_result$taxon_clean & !is.na(TraitName)) %>%
      select(AccSpeciesName, TraitName, mean_StdValue, sd_StdValue, mean_ErrorRisk, sd_ErrorRisk),
    by = c("taxon_clean" = "AccSpeciesName"),
    relationship = "many-to-many"
  )

#Calcular CWM por C/NC
#QUE HACER SI NO TENGO ABUNDANCIA?
# CWM <- 

#RASGOS CUANTITATIVOS

#install.packages("devtools")
#devtools::install_github("bmaitner/RBIEN")
#https://github.com/bmaitner/RBIEN/blob/master/tutorials/RBIEN_tutorial.Rmd
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12861
#library(BIEN)
#library(ape)

#?BIEN
#BIEN<-BIEN_trait_list()

##############################
#Landcover
library(terra)
library(sf)

LC_class1 <- rast("/Users/vanepalmav/Documents/MS Ecología/TFM/DATOS/LANDCOVER/Consensus_reduced_class_1.tif")
LC_class2 <- rast("/Users/vanepalmav/Documents/MS Ecología/TFM/DATOS/LANDCOVER/Consensus_reduced_class_2.tif")
LC_class3 <- rast("/Users/vanepalmav/Documents/MS Ecología/TFM/DATOS/LANDCOVER/Consensus_reduced_class_3.tif")

#Cambiar crs de landcover a uno plano 
# Definir el sistema de referencia a mollweide EPSG:54009
mollweide_crs <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"

# Reproyectar cada raster
LC_class1mw <- project(LC_class1, mollweide_crs)
LC_class2mw <- project(LC_class2, mollweide_crs)
LC_class3mw <- project(LC_class3, mollweide_crs)

#Suma de raster para que todos esten en uno
LC <- LC_class1mw+LC_class2mw+LC_class3mw

plot(LC)
LC

#Extraer dataframe de coordenadas por id_comm
coords <- FA_taxon |>
  group_by(id_comm) |>
  summarise(lon = mean(exact_long),
            lat = mean(exact_lat)) |>
  ungroup()

#Conteo coordenadas distintas
coords$union <- paste(coords$lon, coords$lat, sep = "_")
obs <- unique(coords$union)
length(obs) #198 de 354 coordenadas

#Convertir a Spatialpoints
puntos <- vect(coords, geom = c("lon", "lat"), crs = "EPSG:4326")

#Reproyectar 
puntos_mw <- project(puntos, mollweide_crs)

#Interseccion de puntos y raster, para extraer poligonos en el buffer de 10km
# Convertir los puntos a un objeto sf 
puntos_sf <- st_as_sf(puntos_mw)

# Crear un buffer de 10 km alrededor de cada punto
buffer_10km <- st_buffer(puntos_sf, dist = 10000)  # Distancia en metros

# Convertir el raster a polígonos
lC_polygon <- as.polygons(LC)

# Interseccion del buffer con el raster convertido a polígono
polygon_intersect <- st_intersection(buffer_10km, st_as_sf(lC_polygon))

# Visualizar los polígonos resultantes
plot(st_geometry(polygon_intersect), main = "Polígonos Intersectados")
plot(st_geometry(buffer_10km), add = TRUE, border = 'red', lwd = 2)  # Opcional: mostrar el buffer

# Guardar el resultado en un archivo
st_write(polygon_intersect, "ruta/al/resultado_poligonos_intersectados.shp")

#Calcular metricas
#install.packages("landscapemetrics")
library(landscapemetrics)



