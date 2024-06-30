# This script is the second part of the base model prediction pipeline.
# Its goal is to gather the base layer predictions and merge them into
# a single data frame

# Load the packages
library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)


setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")


predictions_folder <- "Data/Model_output/Predictions/Base_predictions/indiv/"

# Create the full empty predictions data frame
predictions <- read.csv(paste0(predictions_folder, "ref_data.csv"),
  row.names = 1
) %>%
  select("act_index")

# Loop over the individual predictions csv files
for (file in list.files(predictions_folder)) {
  # Read the file
  preds <- read.csv(paste0(predictions_folder, file),
    row.names = 1
  )

  # Add the predictions to the full data frame
  predictions <- merge(
    predictions,
    preds,
    by = "act_index",
    all = TRUE
  )
}

predictions <- predictions %>%
  column_to_rownames("act_index") %>%
  rename("SampleID" = Row.names) %>%
  select("SampleID", everything())

# Save the predictions with the row names since they represent the actual index
write.csv(predictions, "Data/Model_output/Predictions/Base_predictions/Full_base_predictions.csv",
  row.names = TRUE
)

write.csv(predictions, "Data/Model_input/Base_predictions/Full_base_predictions.csv",
          row.names = TRUE
)
