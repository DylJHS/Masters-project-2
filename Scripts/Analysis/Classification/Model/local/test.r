library(dplyr)
library(tidyverse)

datetime <- Sys.time() %>%
str_replace_all(" ", "_") %>%
str_replace_all(":", "_")


print(datetime)