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
library(purrr)
library(rgbif)
#https://cran.r-project.org/web/packages/rtry/vignettes/rtry-introduction.html

#CARGAR DATOS
TRY <- rtry_import("~/Documents/MS Ecología/TFM/DATOS/TFM/36144.txt")
FA <- read.csv(file='/Users/vanepalmav/Documents/MS Ecología/TFM/DATOS/TFM/taxon_temp_vanesa.csv',
                 header=TRUE, sep=';') 
FA_taxon <- read.csv(file='/Users/vanepalmav/Documents/MS Ecología/TFM/DATOS/TFM/assemblages_temp_vanesa.csv',
                  header=TRUE, sep=';') 
TRY_disp <- read.csv('/Users/vanepalmav/Documents/MS Ecología/TFM/DATOS/TFM/TRY_dispersion.csv')

#Analizar por trait (ejemplo)
# Filtrar el dataframe para TraitID igual a 28
#unique(TRY$TraitID)
#filtered_df <- TRY %>%
  #filter(TraitID == 4083) 
#609 no sirve Plant propagation type
#4083 pocos datos
#3364 pocos datos
#28 Sindrome de dispersion
#357 Plant vegetation reproduction

# Eliminar filas donde TraitID sea 4083 (SLA),28,
TRY <- TRY %>%
  filter(TraitID != 4083,
         TraitID != 28,
         TraitID != 609,
         TraitID != 357)

#Seleccionar solo ciertas columnas de FA solo cuando coinciden especies de TRY y FA
TRYNF <- TRY[TRY$AccSpeciesName %in% FA$taxon_clean & !is.na(TRY$TraitID),
             c("AccSpeciesName","TraitID","TraitName","DataName", "StdValue","ErrorRisk")]

###LIMPIEZA DE NOMBRE DE ACCSPECIESNAME
#Para TRY_disp en caso de que NA en gbifspecies repetir taxon_original
TRY_disp <- TRY_disp %>%
  mutate(gbif_species = ifelse(is.na(gbif_species) | gbif_species == "NA NA", 
                               taxon_original, gbif_species))

#Solo especies no subespecies
TRYNF <- TRYNF %>%
  mutate(AccSpeciesName = sapply(strsplit(AccSpeciesName, " "), 
                                 function(x) paste(x[1:2], collapse = " ")))

TRY_disp <- TRY_disp %>%
  mutate(gbif_species = sapply(strsplit(gbif_species, " "),
                               function(x) paste(x[1:2], collapse = " ")))

#Solo 1a letra mayuscula
# Definir la función 
mayus_first <- function(text) {
  tolower(text) %>% 
    sub("^(.)", "\\U\\1", ., perl = TRUE)
}

# Aplicarlo a la columna AccSpeciesName
TRYNF <- TRYNF %>% 
  mutate(AccSpeciesName = mayus_first(AccSpeciesName))

#Filtrar por error < 4
TRY_filt<-TRYNF[TRYNF$ErrorRisk < 4, ]

#Filtrar valores NA
summary(TRYNF$StdValue) #NAs son valores categoricos que pueden haber quedado
TRY_filt <- TRY_filt[!is.na(TRY_filt$StdValue),]

####MODA DE DATOS DE TRY_DISP
###Datos repetidos para TRY_disp
# Función para calcular la moda
get_mode <- function(x) {
  ux <- unique(x)  
  tab <- tabulate(match(x, ux))  
  ux[tab == max(tab)] 
}

# Aplicar la moda por especie (gbif_species)
TRY_disp_mode <- TRY_disp %>%
  group_by(gbif_species) %>%  
  summarise(
    mode_trait_value = list(get_mode(trait_value)),  
    .groups = "drop"
  )

###MEDIDAS DE TENDENCIA Y DISPERSION
#Calculo valor medio, desv de StdValue y mean de ErrorRisk
TRY_mean <- TRY_filt %>%
  group_by(AccSpeciesName, TraitName) %>%  
  summarise(
    n_registros = sum(!is.na(StdValue)),                   # Número de registros no NA en StdValue
    mean_StdValue = mean(StdValue, na.rm = TRUE),          # Promedio de StdValue
    sd_StdValue = sd(StdValue, na.rm = TRUE),              # Desviación estándar de StdValue
    mean_ErrorRisk = mean(ErrorRisk, na.rm = TRUE),        # Promedio de ErrorRisk
    sd_ErrorRisk = sd(ErrorRisk, na.rm = TRUE)             # Desviación estándar de ErrorRisk
  ) %>%
  ungroup()

length(unique(TRY_mean$AccSpeciesName))
#Si sd es NA significa que hay solo 1 observacion

#Media por familias y generos
TRY_mean <- TRY_mean %>%
  mutate(
    gbif_data = map(AccSpeciesName, ~ tryCatch(
      name_backbone(name = .x), error = function(e) NULL
    )),
    familia = map_chr(gbif_data, ~ .x$family %||% NA_character_),
    genero = map_chr(gbif_data, ~ .x$genus %||% NA_character_)
  ) %>%
  select(-gbif_data) 

# Cambiar nombres
unique(TRYNF_mean$TraitName)
TRY_mean <- TRY_mean %>%
  mutate(TraitName = case_when(
    TraitName == "Stem specific density (SSD, stem dry mass per stem fresh volume) or wood density" ~ "Wood density",
    TraitName == "Seed germination rate (germination efficiency)" ~ "Seed germination rate",
    TraitName == "Plant biomass and allometry: Seed number per plant" ~ "Seed number per plant",
    TraitName == "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA) of total leaf area" ~ "Specific leaf area",
    TraitName == "Photosynthesis A/Ci curve: photosynthetic rate per leaf area" ~ "Photosynthetic rate per leaf area",
    TRUE ~ TraitName
  ))

# Aplicar el logaritmo y logit
TRY_mean <- TRY_mean %>%
  mutate(mean_StdValue = case_when(
    TraitName == "Seed dry mass" ~ log(mean_StdValue), 
    TraitName == "Seed number per plant" ~ log(mean_StdValue),
    TraitName == "Seed germination rate" ~ qlogis(mean_StdValue),
    TRUE ~ mean_StdValue                                
  ))

#data frame de especies
sp <- unique(TRY_mean$AccSpeciesName) %>%
  as.data.frame() %>%
  rename ("SpeciesName"=".")

#barplot de frecuencias de rasgos
trait_counts2 <- table(TRY_mean$TraitName)

# Crear el barplot sin etiquetas en el eje x
bp <- barplot(trait_counts2, 
              xlab = "TraitName", 
              ylab = "Cantidad de Datos", 
              main = "Distribución de Datos por TraitName", 
              las = 2,    
              xaxt = "n") 

# Agregar etiquetas en diagonal
text(x = bp, 
     y = par("usr")[3] - 1,          
     labels = names(trait_counts2), 
     srt = 45,                      
     adj = 1,                       
     xpd = TRUE,                    
     cex = 0.5)                     

#Summary e histograma de los valores de c/trait
# Resumen por TraitID
summary_by_trait <- TRY_mean %>%
  group_by(TraitName) %>%
  summarise(
    mean_of_means = mean(mean_StdValue, na.rm = TRUE),    
    sd_of_means = sd(mean_StdValue, na.rm = TRUE),        
    min_mean = min(mean_StdValue, na.rm = TRUE),          
    max_mean = max(mean_StdValue, na.rm = TRUE),          
    n_species = n()                                       
  )
# Mostrar resumen
print(summary_by_trait)

# Histograma de mean_StdValue por TraitName
ggplot(TRY_mean, aes(x = mean_StdValue)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black") +  
  facet_wrap(~ TraitName, scales = "free") +   # Un histograma separado por cada TraitName
  labs(title = "Histograma de mean_StdValue por TraitID", x = "mean_StdValue", y = "Frecuencia") +
  theme_minimal()

#Tabla de especies versus frecuencias de traitname
#filtrar datos NA
TRY_mean <- TRY_mean %>%
  filter(!is.na(TraitName) & TraitName != "")

trait_frequencies <- TRY_mean %>%
  count(AccSpeciesName, TraitName, name = "Frequency") %>%
  pivot_wider(
    names_from = TraitName,  
    values_from = Frequency, 
    values_fill = 0          
  )

# Crear el data frame con los porcentajes promedio
trait_perc <- trait_frequencies %>%
  select(-AccSpeciesName) %>% 
  summarise(across(everything(), ~ sum(.x) / nrow(trait_frequencies) * 100))  

################################

#MERGE FA Y FA_taxon
#Filtrar FA, eliminar no rank
table(FA$rank)
FA <- FA %>% 
  filter(rank != "no rank")

#Para subespecies quedarme solo con 2 primeras palabras

FA_merge<- merge(FA[, c("id_comm", "taxon_clean","measurement")], FA_taxon, by = "id_comm", all.x = TRUE)

# Mantener solo las columnas id_study, id_comm, taxon_clean, measurement y age
FA_clean <- FA_merge %>%
  select(id_study,id_comm, taxon_clean, age, measurement)

length(unique(FA_clean$taxon_clean))
FA_clean <- FA_clean %>% 
  mutate(taxon_clean = mayus_first(taxon_clean))


#ANALISIS COLONIZACION (x aparicion)

########TEST CON 1 ESTUDIO (luego automatizar)
#Filtrar por estudio
unique(FA_clean$id_study)
filter_id <- FA_clean %>%
  filter(id_study == "PR_Baeten et al. 2010_Baeten Belgium vegetation recovery") 

#Filtrar por sp
unique(filter_id$taxon_clean)
filter_sp <- filter_id  %>%
  filter(taxon_clean == "Viola arvensis") 

# Generar todas las combinaciones de edades para cada id_study y taxon_clean
age_mix <- filter_sp %>%
  group_by(id_study, taxon_clean) %>%  
  arrange(age) %>%
  summarise(data = list(combn(age, 2, simplify = FALSE)), .groups = "drop") %>%  
  unnest(cols = c(data)) %>%  
  rename(age_mix = data) %>% 
  mutate(age_n = map_dbl(age_mix, 1),  
         age_m = map_dbl(age_mix, 2),
         age_pair = paste(age_n, age_m, sep = "_"),
         age_diff = abs(age_n - age_m))  

#dataframe de combinaciones de años que sirven
age_unique <- unique(filter_sp$age) %>%
  as.data.frame() %>%
  arrange (.) %>%
  rename("age" = ".")

#vector de pares de años ordenados
años<- paste(age_unique$age[-length(age_unique$age)], age_unique$age[-1], sep = "_")

#filtrar dataframe age_mix por coincidencias
age_mix_filt <- age_mix %>%
  filter(age_pair %in% años)

# Unir con las mediciones de cada combinación
FA_comparisons <- age_mix_filt %>%
  left_join(filter_sp, by = c("id_study", "taxon_clean", "age_n" = "age"),relationship = "many-to-many") %>%
  rename(measurement_n = measurement) %>%
  left_join(filter_sp, by = c("id_study", "taxon_clean", "age_m" = "age"),relationship = "many-to-many") %>%
  rename(measurement_m = measurement) %>%
  mutate(
    comparacion = case_when(
      age_n == age_m ~ NA_character_,  # Si las edades son iguales, comparacion es NA
      measurement_n != 0 & measurement_m != 0 ~ "permanece",       # Permanece
      measurement_n == 0 & measurement_m == 0 ~ "no aparece",     # No aparece
      measurement_n == 0 & measurement_m != 0 ~ "aparece",        # Aparece
      measurement_n != 0 & measurement_m == 0 ~ "desaparece",     # Desaparece
      TRUE ~ NA_character_                                       
    ),
  )

#CODIGO PARA TODOS LOS ESTUDIOS Y SP
    #FA_result<-

# Agrupar por id_study y taxon_clean, luego aplicar la función para cada combinación
FA_result <-
  FA_clean %>%
  group_by(id_study, taxon_clean) %>%  # Agrupar por id_study y taxon_clean
  nest() %>%  # Anidar los datos por grupo
  ungroup() %>%  # Desagrupar después de hacer el `nest()`
  ###############ERROR SOS

str(FA_result)

###################
#UNIR FA_coloniza CON TRYNF
TRY_FA_result <- FA_clean %>%
  inner_join(
    TRY_mean %>%
      filter(AccSpeciesName %in% FA_clean$taxon_clean & !is.na(TraitName)) %>%
      select(AccSpeciesName, TraitName, mean_StdValue, sd_StdValue, mean_ErrorRisk, sd_ErrorRisk),
    by = c("taxon_clean" = "AccSpeciesName"),
    relationship = "many-to-many"
  )

trait_frequencies2 <- TRY_FA_result %>%
  count(taxon_clean, TraitName, name = "Frequency") %>%
  pivot_wider(
    names_from = TraitName,  
    values_from = Frequency, 
    values_fill = 0          
  )

trait_accum <- trait_frequencies2 %>%
  select(-taxon_clean) %>% 
  summarise(across(everything(), ~ sum(.x != 0)))

#UNIR FA_result CON TRY_cuali
TRY_FA_disp <- FA_result %>%
  inner_join(
    TRY_disp %>%
      filter(AccSpeciesName %in% FA_result$taxon_clean & !is.na(TraitName)) %>%
      select(AccSpeciesName, TraitName, mean_StdValue, sd_StdValue, mean_ErrorRisk, sd_ErrorRisk),
    by = c("taxon_clean" = "AccSpeciesName"),
    relationship = "many-to-many"
  )

#test 1 estudio
#study1 <- filter
#boxplot
study1 <- TRY_FA_result %>%
  filter(id_study == "Nombre_Estudio")

study1 <- study1 %>%
  mutate(TraitName = gsub(" ", "\n", TraitName))  # Reemplaza espacios por saltos de línea

ggplot(study1, aes(x = comparacion, y = mean_StdValue)) +
  geom_boxplot() +
  facet_grid(TraitName ~ age, scales = "free") +  # Facetas por 'TraitName' y 'age'
  labs(x = "Comparación", 
       y = "Mean StdValue", 
       title = "Boxplot de Comparación vs. Mean StdValue, separado por age y TraitName") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Textos verticales
    strip.text.y = element_text(size = 10)  # Tamaño del texto en las facetas (filas)
  )


#Calcular CWM por C/NC
#QUE HACER SI NO TENGO ABUNDANCIA?

#RASGOS CUANTITATIVOS

# Crear el gráfico jitter
ggplot(TRYNF, aes(x = comparación, y = DataName)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.6) +  # Ajusta la dispersión y transparencia
  labs(
    x = "Comparación",
    y = "DataName",
    title = "Relación entre Comparación y DataName"
  ) +
  theme_minimal()

##############################ANALISIS ESPACIAL
library(terra)
library(sf)

#Landcover
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

#Filtrar valores, conservar solo >60% cobertura
LC_filtrado <- classify(LC, rcl = matrix(c(-Inf, 60, NA), ncol = 3, byrow = TRUE))
#LC binario (1=bosque)
LC_mayores_60 <- classify(LC, rcl = matrix(c(-Inf, 60, NA, 60, Inf, 1), ncol = 3, byrow = TRUE))

#writeRaster(LC_mayores_60, "LC_bin.tif", overwrite = TRUE)

#Extraer dataframe de coordenadas por id_comm
coords <- FA_taxon |>
  group_by(id_comm) |>
  summarise(lon = mean(exact_long),
            lat = mean(exact_lat)) |>
  ungroup()

#Conteo coordenadas distintas
coords$union <- paste(coords$lon, coords$lat, sep = "_")
length(unique(coords$union))
#198 de 354 coordenadas

#Convertir a Spatialpoints
puntos <- vect(coords, geom = c("lon", "lat"), crs = "EPSG:4326")
#Reproyectar 
puntos_mw <- project(puntos, mollweide_crs)
#Interseccion de puntos y raster para extraer poligonos en el buffer de 10km
# Convertir los puntos a sf 
puntos_sf <- st_as_sf(puntos_mw)
# Crear un buffer de 10 km alrededor de cada punto
buffer_10km <- st_buffer(puntos_sf, dist = 10000)  # Distancia en metros

######CALCULO METRICAS
library(landscapemetrics)
#https://r-spatialecology.github.io/landscapemetrics/

# Lista para almacenar resultados de las métricas de paisaje
resultados_metrica <- list()

# Iterar sobre cada buffer
for (i in 1:nrow(buffer_10km)) {
  # Extraer el buffer actual
  buffer_actual <- buffer_10km[i, ]
  # Recortar y enmascarar el raster según el buffer actual
  LC_buffer <- crop(LC_mayores_60, vect(buffer_actual))
  LC_buffer <- mask(LC_buffer, vect(buffer_actual))
  # Verificar si hay suficiente cobertura forestal para calcular métricas
  if (all(is.na(values(LC_buffer)))) {
    next  # Saltar al siguiente buffer si el raster está vacío
  }
  # Convertir el raster
  LC_raster <- as(LC_buffer, "SpatRaster")
  # Calcular las métricas de paisaje a nivel de clase
  metrica_area <- lsm_c_area_mn(landscape = LC_raster) %>%
    rename(area_mn = value)
  metrica_enn <- lsm_c_enn_mn(landscape = LC_raster) %>%
    rename(enn_mn = value)
  metrica_shape <- lsm_c_shape_mn(landscape = LC_raster) %>%
    rename(shape_mn = value)
  # Combinar todas las métricas en un único data frame
  metrica_total <- full_join(metrica_cpland, metrica_area, by = "class") %>%
    full_join(metrica_enn, by = "class") %>%
    full_join(metrica_shape, by = "class")
  # Agregar id del buffer y columna `id_comm`
  metrica_total <- metrica_total %>%
    mutate(buffer_id = i,
           id_comm = buffer_actual$id_comm)
  # Almacenar en list
  resultados_metrica[[i]] <- metrica_total
}

# Combinar resultados en un único data frame
resultados_df <- bind_rows(resultados_metrica)
print(areas_clases_por_buffer)

#% DE COBERTURA POR BUFFER
resultados_cpland <- list()

for (i in 1:nrow(buffer_10km)) {
  buffer_actual <- buffer_10km[i, ]
  LC_filtrado_buffer <- crop(LC_filtrado, vect(buffer_actual))
  LC_filtrado_buffer <- mask(LC_filtrado_buffer, vect(buffer_actual))
  if (all(is.na(values(LC_filtrado_buffer)))) {
    next  
  }
  # Calcular el área total del buffer en m2
  area_buffer <- as.numeric(st_area(buffer_actual))
  # Calcular CPLAND
  celdas_valores <- values(LC_filtrado_buffer)
  celdas_area <- prod(res(LC_filtrado))  # Area  c/pixel en m2
  promedio_ponderado <- mean(celdas_valores, na.rm = TRUE)  # Promedio 
  # Convertir a % respecto al buffer completo
  cpland_global <- (promedio_ponderado / 100) * 100  
  # Guardar CPLAND y el id del buffer
  resultados_cpland[[i]] <- data.frame(
    buffer_id = i,
    id_comm = buffer_actual$id_comm,
    cpland_global = cpland_global
  )
}

# Combinar los resultados de CPLAND
cpland_df <- bind_rows(resultados_cpland)




