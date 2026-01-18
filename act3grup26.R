############################################################
### Análisis de un conjunto de datos de origen biológico mediante técnicas de machine learning supervisadas y no supervisadas
### 1) Reducción dimensionalidad: PCA, MDS, Isomap,LLE,LE,MVU,UMAP
### 2) Clusterización: K-means, Dendogramas, Heatmap
### 3) Supervisado: SVM, KNN, SVM-Kernel,Árboles de decisión,Random Forest + métricas (CM, Precision, Sensitivity, Specificity, F1)
############################################################

library(tidyverse)
library(caret)
library(factoextra)
library(randomForest)
library(ggplot2)
library(stats)  
library(Rtsne)

# ###############################################################################
# Generación de un único dataframe
# ###############################################################################

f_c <- list.files(pattern="classes.csv$", ignore.case=TRUE)
df_c<-read.csv(f_c, header = FALSE,sep = ";",col.names = c("ID","Clase"))


colnames_g <- readLines("column_names.txt")
f_g <- list.files(pattern="gene_expression.csv$", ignore.case=TRUE)

df <- df_c %>% # Se unifica el dataframe que contiene las clases con el que contiene la expresión génica
  mutate(read.csv(
    f_g,
    header = FALSE,
    # Se asigna el título (nombre de gen) de cada columna del dataframe que contiene los datos de expresión génica
    col.names = colnames_g,
    sep = ";"
  ))

# ###############################################################################
# Imputación de NAs
# ###############################################################################

df_NAs <- df %>%
  select(where( ~ any(is.na(.)))) # La imputación no procede puesto que el dataset no contiene valores faltantes

# ###############################################################################
# Exploración de los datos
# ###############################################################################

table (df$Clase) # Evaluar el balance las clases

genes <- df %>% select(DUOXA1:TTC31) # dataframe con los datos de expresión génica de cada gen
range(genes, na.rm = TRUE) # Consultar el intervalo de valores de expresión génica
res_gen<-t(summary(genes)) # Resumen de estadistícos descriptivos para cada gen


# Gráfico de Densidad para evaluar la distribución global de los valores de expresión génica

g_lar <- genes %>% #Convertir a formato largo para ggplot
  pivot_longer(everything(), names_to = "Gen", values_to = "Expresion")

ggplot(g_lar, aes(x = Expresion)) +
  geom_density(fill = "orange", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribución de la Expresión Génica", 
       x = "Nivel de Expresión", y = "Densidad") 

# Histograma de los ceros
# 
porcentaje_ceros <- colMeans(genes == 0) * 100 # Cálculolo del % de ceros para cada gen

# Histograma para identificar genes con altos porcentaje de ceros
hist(porcentaje_ceros, 
     main = "Histograma de Sparsity (Ceros)",
     xlab = "% de ceros por gen", 
     col = "lightgreen")


# Se genera un nuevo dataframe donde se excluyen los genes con % de ceros superior a 80 (por considerarse menos informativos)

df_filtrado <- bind_cols(df %>% select(ID, Clase), genes[, porcentaje_ceros <=80])

# Escalado de datos 

X_scaled <- df_filtrado %>%
  select(where(is.numeric)) %>% scale() #Aplicación de Z-score (Media=0, SD=1) a la expresión génica


# ###############################################################################
# Reducción dimensionalidad....................................................
# ###############################################################################

# PCA--------------------------------------------------------------------------

# Se aplica PCA a la matriz de datos previamente escalada
pca <- prcomp(X_scaled, center=FALSE)

# Se extraen coordenadas de las muestras en los dos primeros componentes
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Clase = df_filtrado$Clase #Agrupar según clase
)

# Visualización de la estructura global de los datos con PCA
ggplot(pca_df, aes(PC1, PC2, color = Clase)) +
  geom_point(size = 2, alpha = 0.85) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA (PC1 vs PC2)", x = "PC1", y = "PC2")



# ###############################################################################
# Clusterización
# ###############################################################################


# ###############################################################################
# Métodos de aprendizaje supervisado
# ###############################################################################
