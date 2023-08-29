
#############################################################################
# Importing data: 

library(readxl)

# Define the file paths
# DG_path <- "/Users/lukenikolic/Library/Mobile Documents/com~apple~CloudDocs/Documents/MDS DISS/DiGenova_database_wt%_RAI(mol%)_K#(mol%).xlsx"
# Popa_path <- "/Users/lukenikolic/Library/Mobile Documents/com~apple~CloudDocs/Documents/MDS DISS/popa_et_al_2021_wt%_RAI(mol%)_K#(mol%).xlsx"

# Get DG sheets and assign to variables

FMQ_DG_df <- data.frame(read_excel(DG_path,
                                   sheet = "FMQ"))


# Get Popa sheets and assign to variables




FMQ_popa_df <- data.frame(read_excel(Popa_path,
                                     sheet = "FMQ"))




# Add eruptive_type binary column and move vairables to global env. 

# DG 


DG_names <- c("FMQ_DG_df")


for (name in DG_names) {
  df <- get(name)
  df$eruptive_type <- c(rep(0, 571), rep(1, nrow(df) - 571))
  assign(name, df, envir = .GlobalEnv)
}


columns_to_remove <- c(1, 2, 3, 4, 10, 13, 16, 17, 18, 19)

for (name in DG_names) {
  df <- get(name)
  df <- df[,-columns_to_remove]
  assign(name, df, envir = .GlobalEnv)
}



# Same for Popa

Popa_names <- c("FMQ_popa_df")


for (name in Popa_names) {
  df <- get(name)
  df$eruptive_type <- ifelse(df$'Eruptive.style' == "effusive", 0, 1)
  assign(name, df, envir = .GlobalEnv)
}


columns_popa_remove <- c(1, 2, 3, 4, 5, 6, 16) # shift by 1 for 1-based indexing

for (name in Popa_names) {
  df <- get(name)
  df <- df[,-columns_popa_remove]
  assign(name, df, envir = .GlobalEnv)
}


# Now Popa and DG can be concatenated


combined_FMQ_df <- rbind(FMQ_DG_df, FMQ_popa_df)




dataframes_buffers <- list(combined_FMQ_df)



names_buffers <- c("combined_FMQ_clean_range_df")

for (i in seq_along(dataframes_buffers)) {
  df <- dataframes_buffers[[i]]
  
  
  # Select all columns except the last one and check if any of them are equal to 0
  zero_in_columns <- apply(df[, -ncol(df)], 1, function(x) any(x == 0))
  
  # Keep only the rows where none of the selected columns have a value of 0
  clean_df <- df[!zero_in_columns,]
  
  # Reset row names to make them sequential
  rownames(clean_df) <- NULL
  
  # Assign the cleaned dataframe to a new variable in the global environment
  assign(names_buffers[i], clean_df, envir = .GlobalEnv)
}


# Now buffer data sets in right form for clr transformation

#############################################################################

### CLR transformation

## Reduced data set

library(compositions)



dataframes_buffers_ToBeTransformed <- list(combined_FMQ_clean_range_df)

names_buffers_ToBeTransformed <- c("combined_FMQ_clr_transformed_df_range")

for (i in seq_along(dataframes_buffers_ToBeTransformed)) {
  df <- dataframes_buffers_ToBeTransformed[[i]]
  
  # Extract the columns "RAI.mol..", "K..mol.." and "eruptive_type"
  RAI_column <- df["RAI.mol.."]
  k_column <- df["K..mol.."]
  eruptive_type_column <- df["eruptive_type"]
  
  # Extract the compositional data (first 9 columns)
  comp_data <- df[, 1:9]
  
  # Perform CLR transformation on the compositional data
  clr_transformed_data <- clr(comp_data)
  
  # Convert the transformed data back to a data frame
  clr_transformed_df <- data.frame(clr_transformed_data)
  
  # Name the columns as in the original data frame
  colnames(clr_transformed_df) <- colnames(comp_data)
  
  # Add back the columns "RAI.mol..", "K..mol.." and "eruptive_type" to the data frame
  clr_transformed_df["RAI.mol.."] <- RAI_column
  clr_transformed_df["K..mol.."] <- k_column
  clr_transformed_df["eruptive_type"] <- eruptive_type_column
  
  # Assign the transformed data frame to the global environment
  assign(names_buffers_ToBeTransformed[i], clr_transformed_df, envir = .GlobalEnv)
}


##################################################################################



## FMQ: 

# Set seed for reproducibiilty of PCA

set.seed(69)

### PCA for FMQ data set: 

pr.out_FMQ <- prcomp(combined_FMQ_clr_transformed_df_range[,-12] , scale =TRUE)
names(pr.out_FMQ)
summary(pr.out_FMQ)

# eigenvalues
pr.out_FMQ$sdev^2

# PC loadings
pr.out_FMQ$rotation

# PC score
pr.out_FMQ$x

# Generate scree plot

library(ggplot2)
library(factoextra)
fviz_screeplot(pr.out_FMQ, addlabels = TRUE)+
  ggtitle("Proportion of Variance Explained by Each Principal Component")+
  xlab("Dimensions (number of principal components)")+
  ylab("Percentage of Explained Variance in the Dataset")

# Scree plot and summary of pr.out_reduced suggests first 4 PCs keeps 
# 89% of data, and first 3 explains 82% of data therefore will choose 4. 



var_FMQ <- get_pca_var(pr.out_FMQ)

var_FMQ$contrib

# Plot contribution of variables for each PC (1-4)

fviz_contrib(pr.out_FMQ, choice = "var", axes = 1, top = 10)
fviz_contrib(pr.out_FMQ, choice = "var", axes = 2, top = 10)
fviz_contrib(pr.out_FMQ, choice = "var", axes = 3, top = 10)
fviz_contrib(pr.out_FMQ, choice = "var", axes = 4, top = 10)

library(ggplot2)

fviz_contrib(pr.out_FMQ, choice = "var", axes = 1, top = 10) + 
  ggtitle("Contribution of Original clr-transformed Variables to Principal Component 1")+
  xlab("Predictor Variables")+
  ylab("Contribution to Variance (%)")

fviz_contrib(pr.out_FMQ, choice = "var", axes = 2, top = 10) + 
  ggtitle("Contribution of Original clr-transformed Variables to Principal Component 2")+
  xlab("Predictor Variables")+
  ylab("Contribution to Variance (%)")

fviz_contrib(pr.out_FMQ, choice = "var", axes = 3, top = 10) + 
  ggtitle("Contribution of Original clr-transformed Variables to Principal Component 3")+
  xlab("Predictor Variables")+
  ylab("Contribution to Variance (%)")

fviz_contrib(pr.out_FMQ, choice = "var", axes = 4, top = 10) + 
  ggtitle("Contribution of Original clr-transformed Variables to Principal Component 4")+
  xlab("Predictor Variables")+
  ylab("Contribution to Variance (%)")



# Biplot for first two PCs >> way of representing data in 2-D space

library(factoextra)
library(ggplot2)


# Create a color palette to map "effusive" to red and "explosive" to blue

# Creating copy of pr.out to change RAI and K# names 

combined_FMQ_clr_transformed_df_range_copy <- combined_FMQ_clr_transformed_df_range
names(combined_FMQ_clr_transformed_df_range_copy)[names(combined_FMQ_clr_transformed_df_range_copy) == "RAI.mol.."] <- "RAI"
names(combined_FMQ_clr_transformed_df_range_copy)[names(combined_FMQ_clr_transformed_df_range_copy) == "K..mol.."] <- "K#"


pr.out_FMQ_copy <- prcomp(combined_FMQ_clr_transformed_df_range_copy[,-12] , scale =TRUE)

color_palette <- c("effusive" = "skyblue", "explosive" = "lightcoral")

# Create a factor variable for eruptive_type
eruptive_type_factor_FMQ <- factor(combined_FMQ_clr_transformed_df_range$eruptive_type,
                                   levels = c(0, 1),
                                   labels = c("effusive", "explosive"))



library(factoextra)
library(ggplot2)

eruptive_type_factor_FMQ <- factor(combined_FMQ_clr_transformed_df_range$eruptive_type,
                                   levels = c("Class0", "Class1"),
                                   labels = c("effusive", "explosive"))



library(ggplot2)

biplot_plot <- fviz_pca_biplot(pr.out_FMQ_copy,
                               col.ind = eruptive_type_factor_FMQ,
                               palette = color_palette,
                               addEllipses = FALSE,
                               invisible = "var",
                               label = "none",
                               col.var = "black", repel = FALSE,
                               legend.title = "Type of Eruption") +
  xlab("PC 1 (53.4% variance)") +
  ylab("PC 2 (18.8% variance)") +
  guides(fill = guide_legend(title = "Type of Eruption")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

biplot_plot

variable_biplot_plot <- fviz_pca_biplot(pr.out_FMQ_copy,
                                        geom = "arrow", # To represent variables as arrows
                                        invisible = "ind", # Hide individuals
                                        col.var = "black", repel = TRUE,
                                        show.legend = FALSE) + # No need for the legend
  xlab("PC 1 (53.4% variance)") +
  ylab("PC 2 (18.8% variance)") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 14)
    
    
  )


variable_biplot_plot


# Save 











# Now need to see how to analyse biplot << come back to this




set.seed(69)

library(caret)
library(MLmetrics)
library(pROC)
library(caTools)
library(Metrics)


combined_FMQ_clr_transformed_df_range$eruptive_type <- as.factor(combined_FMQ_clr_transformed_df_range$eruptive_type)
levels(combined_FMQ_clr_transformed_df_range$eruptive_type) <- c("Class0", "Class1")


folds_FMQ <- createFolds(combined_FMQ_clr_transformed_df_range$eruptive_type, k = 10)

# Initialize variables to hold metrics for each fold
precision_class0_FMQ <- recall_class0_FMQ <- f1_class0_FMQ <- 0
precision_class1_FMQ <- recall_class1_FMQ <- f1_class1_FMQ <- 0
accuracy_total_FMQ <- 0

# Perform cross-validation
for (fold in folds_FMQ) {
  test_indexes <- fold
  train_indexes <- setdiff(1:nrow(combined_FMQ_clr_transformed_df_range), test_indexes)
  
  # Divide into training and testing sets
  train <- combined_FMQ_clr_transformed_df_range[train_indexes, ]
  test <- combined_FMQ_clr_transformed_df_range[test_indexes, ]
  
  # Apply PCA to training data
  pr.out_FMQ <- prcomp(train[,-12], scale = TRUE)
  
  # Transform training data using the first 4 principal components
  train_transformed <- data.frame(pr.out_FMQ$x[, 1:4])
  train_transformed$eruptive_type <- train$eruptive_type
  
  # Transform testing data
  test_transformed <- as.data.frame(predict(pr.out_FMQ, newdata = test[,-12]))
  test_transformed <- test_transformed[, 1:4]
  test_transformed$eruptive_type <- test$eruptive_type
  
  # Fit logistic regression model
  model <- glm(eruptive_type ~ PC1 + PC2 + PC3 + PC4, data = train_transformed, family = "binomial")
  
  # Predict probabilities on the test set
  probs <- predict(model, test_transformed, type = "response")
  
  # Convert probabilities to class labels
  predictions <- ifelse(probs > 0.5, "Class1", "Class0")
  predictions <- factor(predictions, levels = c("Class0", "Class1"))
  test_labels <- factor(test_transformed$eruptive_type, levels = c("Class0", "Class1"))
  
  # Compute confusion matrix
  cm <- confusionMatrix(predictions, test_labels)
  
  # Compute precision, recall, and F1 score for each class
  precision_class0_FMQ <- precision_class0_FMQ + cm$byClass[5]
  recall_class0_FMQ <- recall_class0_FMQ + cm$byClass[6]
  f1_class0_FMQ <- f1_class0_FMQ + cm$byClass[7]
  
  precision_class1_FMQ <- precision_class1_FMQ + cm$byClass[1]
  recall_class1_FMQ <- recall_class1_FMQ + cm$byClass[2]
  f1_class1_FMQ <- f1_class1_FMQ + cm$byClass[3]
  
  # Compute accuracy
  accuracy_total_FMQ <- accuracy_total_FMQ + cm$overall[1]
}

# Calculate average precision, recall, F1 score, and accuracy across 10 folds
average_precision_class0_FMQ <- precision_class0_FMQ / 10
average_recall_class0_FMQ <- recall_class0_FMQ / 10
average_f1_class0_FMQ <- f1_class0_FMQ / 10

average_precision_class1_FMQ <- precision_class1_FMQ / 10
average_recall_class1_FMQ <- recall_class1_FMQ / 10
average_f1_class1_FMQ <- f1_class1_FMQ / 10

average_accuracy_FMQ <- accuracy_total_FMQ / 10

# Print the results
cat("Average metrics for Class0:\n")
cat("Precision:", average_precision_class0_FMQ, "\nRecall:", average_recall_class0_FMQ, "\nF1 Score:", average_f1_class0_FMQ, "\n")
cat("\nAverage metrics for Class1:\n")
cat("Precision:", average_precision_class1_FMQ, "\nRecall:", average_recall_class1_FMQ, "\nF1 Score:", average_f1_class1_FMQ, "\n")
cat("\nAverage Accuracy:", average_accuracy_FMQ)

# Found average accuracy is roughly 80% 


library(caret)
library(pROC)
library(ggplot2)

folds_roc_FMQ <- createFolds(combined_FMQ_clr_transformed_df_range$eruptive_type, k = 10)

# Initialize variables to hold the ROC curve and AUC for each fold
fpr_list_FMQ <- vector("list", length(folds_roc_FMQ))
tpr_list_FMQ <- vector("list", length(folds_roc_FMQ))
auc_list_FMQ <- numeric(length(folds_roc_FMQ))

# Perform cross-validation
for (i in seq_along(folds_roc_FMQ)) {
  test_indexes <- folds_roc_FMQ[[i]]
  train_indexes <- setdiff(1:nrow(combined_FMQ_clr_transformed_df_range), test_indexes)
  
  # Divide into training and testing sets
  train <- combined_FMQ_clr_transformed_df_range[train_indexes, ]
  test <- combined_FMQ_clr_transformed_df_range[test_indexes, ]
  
  # Apply PCA to training data
  pr.out_FMQ <- prcomp(train[,-12], scale = TRUE)
  
  # Transform training data using the first 4 principal components
  train_transformed <- data.frame(pr.out_FMQ$x[, 1:4])
  train_transformed$eruptive_type <- train$eruptive_type
  
  # Transform testing data
  test_transformed <- as.data.frame(predict(pr.out_FMQ, newdata = test[,-12]))
  test_transformed <- test_transformed[, 1:4]
  test_transformed$eruptive_type <- test$eruptive_type
  
  # Fit logistic regression model
  model <- glm(eruptive_type ~ PC1 + PC2 + PC3 + PC4, data = train_transformed, family = "binomial")
  
  # Predict probabilities on the test set
  probs <- predict(model, test_transformed, type = "response")
  
  # Compute ROC curve for this fold
  roc_obj <- roc(test_transformed$eruptive_type, probs, levels = c("Class1", "Class0"))
  fpr_list_FMQ[[i]] <- coords(roc_obj, "all", ret = "1-specificity")[, 1]
  tpr_list_FMQ[[i]] <- coords(roc_obj, "all", ret = "sensitivity")[, 1]
  
  # Compute AUC for this fold
  auc_list_FMQ[i] <- roc_obj$auc
  
}

# Compute the average AUC
avg_auc_FMQ <- mean(auc_list_FMQ)

# Compute the average ROC curve
common_fpr_FMQ <- seq(0, 1, by = 0.01)
avg_tpr_FMQ <- sapply(common_fpr_FMQ, function(x) mean(sapply(1:length(fpr_list_FMQ), function(i) {
  fpr_FMQ <- fpr_list_FMQ[[i]]
  tpr_FMQ <- tpr_list_FMQ[[i]]
  approx(fpr_FMQ, tpr_FMQ, xout = x)$y
})))

# Print the average AUC
cat("Average AUC across 10 folds:", avg_auc_FMQ, "\n")

avg_roc_label <- paste0("Average ROC Curve (AUC = ", round(avg_auc_FMQ, 2), ")")

p <- ggplot() + 
  geom_line(aes(x = common_fpr_FMQ, y = avg_tpr_FMQ, color = avg_roc_label), linetype = "solid", size = 1) +
  scale_color_manual(
    values = c(avg_roc_label = "black"), 
    name = ""
  ) +
  theme_minimal()
p

p <- p + 
  geom_line(aes(x = c(0, 1), y = c(0, 1), color = "Random Classifier (AUC = 0.5)"), linetype = "dashed", size = 1) +
  scale_color_manual(
    values = c(avg_roc_label = "black", "Random Classifier (AUC = 0.5)" = "red"), 
    name = ""
  )
p

p <- p + 
  geom_line(aes(x = c(0, 0, 1), y = c(0, 1, 1), color = "Perfect Classifier (AUC = 1.0)"), linetype = "dashed", size = 1) +
  scale_color_manual(
    values = c(avg_roc_label = "black", 
               "Random Classifier (AUC = 0.5)" = "red", 
               "Perfect Classifier (AUC = 1.0)" = "green"),
    name = ""
  )
p

p <- p +
  scale_linetype_manual(
    values = c(avg_roc_label = "solid", 
               "Random Classifier (AUC = 0.5)" = "dashed", 
               "Perfect Classifier (AUC = 1.0)" = "dashed"),
    name = ""
  ) +
  labs(x = "False Positive Rate", y = "True Positive Rate", title = "Average ROC Curve from 10-fold Cross-Validation (Logistic Regression, FMQ Data Set)") +
  theme(legend.title = element_blank())
p

# Going to export data from R into python: 

# Combining vectors into a dataframe
df_export <- data.frame(
  FPR = common_fpr_FMQ,
  TPR = avg_tpr_FMQ
)

# Writing the dataframe to a CSV file











#############################################################################


##############################################################################

# EDA - combined FMQ data 



library(tidyverse)
library(hrbrthemes)
library(viridis)

library(tidyverse)
library(hrbrthemes)
library(viridis)


## Plot for class imbalance - data set used is one just before clr trans

library(ggplot2)


response_var_FMQ <- combined_FMQ_clean_range_df[[ncol(combined_FMQ_clean_range_df)]]

# Create a data frame for plotting
class_imbalance_df_FMQ <- data.frame(
  class = factor(response_var_FMQ),
  count = 1
)


plot <- ggplot(class_imbalance_df_FMQ, aes(x = class, fill = class)) +
  geom_bar(aes(y = ..count..), alpha = 0.7) +
  geom_text(stat='count', aes(label=..count.., y=..count..), vjust=-0.5) +
  scale_fill_manual(values = c("skyblue", "lightcoral"), labels = c("0 (effusive)", "1 (explosive)"),
                    name = "Eruption Type") +
  labs(
    title = "Class Distribution for Type of Eruption",
    x = "Type of Eruption",
    y = "Number of Instances",
    fill = "Class"
  ) +
  theme_classic()

# Get the maximum count to set y limits
max_count <- max(summary(as.factor(class_imbalance_df_FMQ$class)))

# Adjust the y limits
plot + coord_cartesian(ylim = c(0, max_count * 1.1))


# Class imbalance graph shows there to be ratio of 
# roughly 3:2 << motivation for straitifed k-fold

predictor_variables_clr <- combined_FMQ_clr_transformed_df_range[, -ncol(combined_FMQ_clr_transformed_df_range)]

# Convert the data to 'long' format
long_data_clr <- predictor_variables_clr %>%
  gather(key = "variable", value = "value")

# Create a ggplot showing the density distribution of each predictor variable
ggplot(long_data_clr, aes(x = value)) +
  geom_density(fill = "blue", alpha = 0.5) +
  facet_wrap(~ variable, scales = "free") + # Create a grid of plots
  labs(
    x = "Value",
    y = "Density",
    title = "Distributions of Predictor Variables"
  ) +
  theme_minimal()


######## Redo density plots with RAI and K# given proper names

# Create a copy of the dataset
predictor_variables_clr_copy <- predictor_variables_clr

# Change the names in the copy
colnames(predictor_variables_clr_copy)[colnames(predictor_variables_clr_copy) == "K..mol.."] <- "K#"
colnames(predictor_variables_clr_copy)[colnames(predictor_variables_clr_copy) == "RAI.mol.."] <- "RAI"

# Convert the data to 'long' format using the modified dataset
long_data_clr <- predictor_variables_clr_copy %>%
  gather(key = "variable", value = "value")

# Create a ggplot showing the density distribution of each predictor variable
ggplot(long_data_clr, aes(x = value)) +
  geom_density(fill = "blue", alpha = 0.5) +
  facet_wrap(~ variable, scales = "free") + # Create a grid of plots
  labs(
    x = "CLR-Transformed Value",
    y = "Density",
    title = "Distributions of Predictor Variables"
  ) +
  theme_minimal()

ggplot(long_data_clr, aes(x = value)) +
  geom_density(fill = "blue", alpha = 0.5) +
  facet_wrap(~ variable) + # Create a grid of plots
  labs(
    x = "CLR-Transformed Value",
    y = "Density",
    title = "Distributions of Predictor Variables"
  ) +
  theme_minimal()


################ Re-do denisity plots with original combined data not clr

predictor_variables_org <- combined_FMQ_clean_range_df[, -ncol(combined_FMQ_clean_range_df)]
# Change the names in the copy
colnames(predictor_variables_org)[colnames(predictor_variables_org) == "K..mol.."] <- "K#"
colnames(predictor_variables_org)[colnames(predictor_variables_org) == "RAI.mol.."] <- "RAI"

# Convert the data to 'long' format using the modified dataset
long_data_org <- predictor_variables_org %>%
  gather(key = "variable", value = "value")

# Create a ggplot showing the density distribution of each predictor variable
ggplot(long_data_org, aes(x = value)) +
  geom_density(fill = "blue", alpha = 0.5) +
  facet_wrap(~ variable, scales = "free") + # Create a grid of plots
  labs(
    x = "Value (wt%)",
    y = "Density",
    title = "Distributions of Predictor Variables"
  ) +
  theme_minimal()





## Correlation matrix and pairs plot for clr-transformed data


library(ggplot2)
library(reshape2)
library(tidyverse)
predictor_variables_names <- colnames(combined_FMQ_clr_transformed_df_range)[colnames(combined_FMQ_clr_transformed_df_range) != "eruptive_type"]

# Calculate Spearman's rank correlation for each predictor variable
spearman_correlations <- sapply(predictor_variables_names, function(var_name) {
  cor.test(combined_FMQ_clr_transformed_df_range[[var_name]], combined_FMQ_clr_transformed_df_range$eruptive_type, method = "spearman")$estimate
})

# Print the result
print(spearman_correlations)

# Interestingly, SiO2 has low correlation, but there might be reason
# for this. 

(average_value_SiO2 <- mean(combined_FMQ_clean_range_df$SiO2))

## Average value for SiO2 in my dataset is 74 wt% - therefore dealing with 
## primarily rhyolitic magmas? 



library(reshape2) 

variables_selected_spear <- c(predictor_variables_names, "eruptive_type")


correlation_matrix_spearman <- cor(combined_FMQ_clr_transformed_df_range[variables_selected_spear], method = "spearman")


melted_corr_matrix_spear <- melt(correlation_matrix_spearman)

# Create a ggplot
ggplot(data = melted_corr_matrix_spear, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value), color = "white") +
  geom_text(aes(label = round(value, 2)), vjust = 1) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("Spearman's Correlation Matrix for Predictor and Response Variables")



### Lattice plot for different box plots for predictor variables

library(ggplot2)
library(tidyr)
library(viridis)

# Pivot the data to a longer format
long_data <- combined_FMQ_clr_transformed_df_range %>%
  select(predictor_variables_names) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")




# Create the ggplot with boxplots
p <- ggplot(long_data, aes(x = Variable, y = Value)) +
  geom_boxplot() +
  facet_wrap(~ Variable, scales = "free") +
  theme(axis.text.x = element_blank())+
  ggtitle("Box Plots for Predictor Variables in CLR-Transformed Data")

# Print the plot
print(p)


## Can't add colour since too many pred variables. Could use viridis but 
## don't like how each colour differes minimally

### Pre-modelling EDA to be continued if needed

# Summary stats table for clr-transformed data
library(stargazer)
stargazer(combined_FMQ_clr_transformed_df_range[,-12])

# Summary table for FMQ_clean_range df - i.e. non-clr transformed 
stargazer(combined_FMQ_clean_range_df[,-12])


# Head of first few rows for non-clr data

stargazer(combined_FMQ_clean_range_df[1:10,], summary = FALSE, rownames = FALSE)










library(tidyverse)
library(ggplot2)

# Transform the dataframe into a long format excluding Volcano and eruptive_type
long_df_python <- python_dataframe %>%
  select(-Volcano, -eruptive_type) %>%
  gather(key = "variable", value = "value")


# Create a boxplot
p <- ggplot(long_df_python, aes(x = variable, y = value)) +
  geom_boxplot(fill = "skyblue") +
  labs(y = "Value", x = "Variable") +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10)
  )

# Display the plot
print(p)




















#################################################################################

## EDA post-modelling: 

# RAI, K#, Na2O, Fe2O3, CaO, SiO2 << important variables for models

# ^^ subject to change


# Investigating interesting variables: 

# Relationship between RAI and response variable found to not be the same as 
# suggested by DG but this is probs because of assumed O2 buffer - models 
# still show RAI to be important which is consistent w DG's findings

### RAI 

group_summary_RAI.mol.. <- combined_FMQ_clean_range_df_copy %>% 
  group_by(eruptive_type) %>%
  summarise(mean = mean(`RAI.mol..`), median = median(`RAI.mol..`))


# Boxplots 

ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `RAI.mol..`, fill = as.factor(eruptive_type))) +
  geom_boxplot() +
  labs(title = "RAI Values for Different Eruptive Types",
       x = "Eruptive Type",
       y = "RAI(mol%)") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type")


# Bplots show generally RAI lower for explosive eruptions than effusive. 
# This is consistent with DGs findings


ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `RAI.mol..`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  labs(title = "Jitter Plot of RAI vs Eruptive Type",
       x = "Eruptive Type",
       y = "RAI(mol%)") +
  theme_minimal()


# Therefore findings partially consistent with numerical threshold set 
# by DG. However not surpising it is diff, since o2 buffer (FMQ) 
# with constant temp has been assumed, which will change vals of FeO and 
# Fe2O3 immensely and therefoe values of RAI. 
# 




### K# 

# K# indicated to be strong predictor variable, will look at relationship 
# with this and binary response variable

### K#

library(ggplot2)




# Create a copy of the dataframe
combined_FMQ_clean_range_df_copy <- combined_FMQ_clean_range_df

# Convert the eruptive_type variable in the copy to a factor
combined_FMQ_clean_range_df_copy$eruptive_type <- as.factor(combined_FMQ_clean_range_df_copy$eruptive_type)



library(dplyr)

# Calculate group-wise average and median
group_summary_K <- combined_FMQ_clean_range_df_copy %>% 
  group_by(eruptive_type) %>%
  summarise(mean = mean(`K..mol..`), median = median(`K..mol..`))

# Define custom labels and colors
custom_labels_copy <- c("0" = "Effusive", "1" = "Explosive")
custom_colors_copy <- c("0" = "skyblue", "1" = "lightcoral")

# Plot
ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `K..mol..`, fill = as.factor(eruptive_type))) +
  geom_boxplot() +
  labs(title = "K#(mol%) Values for Different Eruptive Types",
       x = "Eruptive Type",
       y = "K#(mol%)") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type")



ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `K..mol..`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 0.55, linetype = "dashed", color = "black") +
  labs(title = "Jitterplot of K# vs Eruptive Type",
       x = "Eruptive Type",
       y = "K#(mol%)") +
  theme_minimal()



test_p <- ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `K..mol..`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 1.5) +
  labs(x = "Eruptive Type",
       y = "K#(mol%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 24),
    legend.position = "none",  # This line removes the legend
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 28)
  )


test_p











### Na2O


group_summary_Na2O <- combined_FMQ_clean_range_df_copy %>% 
  group_by(eruptive_type) %>%
  summarise(mean = mean(`Na2O`), median = median(`Na2O`))


# Boxplots 

ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `Na2O`, fill = as.factor(eruptive_type))) +
  geom_boxplot() +
  labs(title = "Na2O Values for Different Eruptive Types",
       x = "Eruptive Type",
       y = "Na2O(wt%)") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type")



ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `Na2O`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "black") +
  labs(title = "Jitterplot of Na2O vs Eruptive Type",
       x = "Eruptive Type",
       y = "Na2O (wt%)") +
  theme_minimal()



group_summary_SiO2 <- combined_FMQ_clean_range_df_copy %>% 
  group_by(eruptive_type) %>%
  summarise(mean = mean(`SiO2`), median = median(`SiO2`))


# Boxplots 

ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `SiO2`, fill = as.factor(eruptive_type))) +
  geom_boxplot() +
  labs(title = "SiO2 Values for Different Eruptive Types",
       x = "Eruptive Type",
       y = "SiO2(wt%)") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type")




ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `SiO2`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 78.5, linetype = "dashed", color = "black")+
  labs(title = "Jitterplot of SiO2 vs Eruptive Type",
       x = "Eruptive Type",
       y = "SiO2 (wt%)") +
  theme_minimal()





### Fe2O3


group_summary_Fe2O3 <- combined_FMQ_clean_range_df_copy %>% 
  group_by(eruptive_type) %>%
  summarise(mean = mean(`Fe2O3`), median = median(`Fe2O3`))


# Boxplots 

ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `Fe2O3`, fill = as.factor(eruptive_type))) +
  geom_boxplot() +
  labs(title = "Fe2O3 Values for Different Eruptive Types",
       x = "Eruptive Type",
       y = "Fe2O3(wt%)") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type")

## Showing Fe2O3 usually higher for effusive eruptions



ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `Fe2O3`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  labs(title = "Jitterplot of Fe2O3 vs Eruptive Type",
       x = "Eruptive Type",
       y = "Fe2O3 (wt%)") +
  theme_minimal()






### CaO


group_summary_CaO <- combined_FMQ_clean_range_df_copy %>% 
  group_by(eruptive_type) %>%
  summarise(mean = mean(`CaO`), median = median(`CaO`))


# Boxplots 

ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `CaO`, fill = as.factor(eruptive_type))) +
  geom_boxplot() +
  labs(title = "CaO Values for Different Eruptive Types",
       x = "Eruptive Type",
       y = "CaO(wt%)") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type")

## Not much difference shown between groups



ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `CaO`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  labs(title = "Jitterplot of CaO vs Eruptive Type",
       x = "Eruptive Type",
       y = "CaO (wt%)") +
  theme_minimal()



### Now done EDA for important variables found from feature selection 
### in models: RAI, K#, Na2O, Fe2O3, SiO2, CaO

########## Now some post modelling EDA for volcanoes that dominant the data

## Read in from csv file from python data as three separate dfs



## Now looking at most important variables but for each volcanic system



## RAI:

## Cordon Caulle 

# Create a copy of the dataframe
cordon_caulle_df_copy <- cordon_caulle_df 

# Convert the eruptive_type variable in the copy to a factor
cordon_caulle_df_copy$eruptive_type <- as.factor(cordon_caulle_df_copy$eruptive_type)


ggplot(cordon_caulle_df_copy, aes(x = eruptive_type, y = `RAI.mol..`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  labs(title = "Jitter Plot of RAI vs Eruptive Type (Cordon Caulle)",
       x = "Eruptive Type",
       y = "RAI(mol%)") +
  theme_minimal()

##  Yellowstone

yellowstone_df_copy <- yellowstone_df 

# Convert the eruptive_type variable in the copy to a factor
yellowstone_df_copy$eruptive_type <- as.factor(yellowstone_df_copy$eruptive_type)


ggplot(yellowstone_df_copy, aes(x = eruptive_type, y = `RAI.mol..`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  labs(title = "Jitter Plot of RAI vs Eruptive Type (Yellowstone)",
       x = "Eruptive Type",
       y = "RAI(mol%)") +
  theme_minimal()



ggplot(yellowstone_df_copy, aes(x = eruptive_type, y = `K..mol..`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 1.5)+
  labs(
    x = "Eruptive Type",
    y = "K#(mol%)") +
  theme_minimal()+
  theme(
    text = element_text(size = 24),
    legend.position = "none",  # This line removes the legend
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 28)
  )

# Non yellowstone for K# 


non_yellow_df_copy <- non_yellow 

# Convert the eruptive_type variable in the copy to a factor
non_yellow_df_copy$eruptive_type <- as.factor(non_yellow_df_copy$eruptive_type)

ggplot(non_yellow_df_copy, aes(x = eruptive_type, y = `K..mol..`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 1.5)+
  labs(
    x = "Eruptive Type",
    y = "K#(mol%)") +
  theme_minimal()+
  theme(
    text = element_text(size = 24),
    legend.position = "none",  # This line removes the legend
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 28)
  )




######### Modified Plots for other important variables

# RAI

ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `RAI.mol..`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1.5) +
  labs(x = "Eruptive Type",
       y = "RAI(mol%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 24),
    legend.position = "none",  # This line removes the legend
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 28)
  )


# Na_2O

ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `Na2O`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "black", size = 1.5) +
  labs(x = "Eruptive Type",
       y = "Na2O(wt%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 24),
    legend.position = "none",  # This line removes the legend
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 28)
  )


# SiO2


ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `SiO2`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 77.5, linetype = "dashed", color = "black", size = 1.5) +
  labs(x = "Eruptive Type",
       y = "SiO2(wt%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 24),
    legend.position = "none",  # This line removes the legend
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 28)
  )

# Fe2o3

ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `Fe2O3`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "black", size = 1.5) +
  labs(x = "Eruptive Type",
       y = "Fe2O3(wt%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 24),
    legend.position = "none",  # This line removes the legend
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 28)
  )

# CaO

ggplot(combined_FMQ_clean_range_df_copy, aes(x = eruptive_type, y = `CaO`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type") +
  labs(x = "Eruptive Type",
       y = "CaO(wt%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 24),
    legend.position = "none",  # This line removes the legend
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 28)
  )




### Bivariate plots of CaO with other predictor variables

## RAI 

ggplot(combined_FMQ_clean_range_df_copy, aes(x = CaO, y = `Fe2O3`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy, labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1.5) +  # Adjust or remove as needed
  labs(x = "CaO", y = "RAI(mol%)") +
  theme_minimal() +
  theme(
    text = element_text(size = 24),
    legend.position = "none",
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 28)
  )












### Lipari 

lipari_df_copy <- lipari_df 

# Convert the eruptive_type variable in the copy to a factor
lipari_df_copy$eruptive_type <- as.factor(lipari_df_copy$eruptive_type)


ggplot(lipari_df_copy, aes(x = eruptive_type, y = `RAI.mol..`)) +
  geom_jitter(aes(color = eruptive_type), width = 0.2, size = 2, alpha = 0.6) +
  scale_color_manual(values = custom_colors_copy,labels = custom_labels_copy, name = "Eruptive Type") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  labs(title = "Jitter Plot of RAI vs Eruptive Type (Lipari)",
       x = "Eruptive Type",
       y = "RAI(mol%)") +
  theme_minimal()
















###########################################################################










































