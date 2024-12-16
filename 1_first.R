---
title: "first"
format: html
params:
  downsample_size: 100
  
---

## 


```{r}
#| title: setup

library(Seurat)
library(tidyseurat)
library(tidyverse)
library(janitor)
library(here)


```



```{r}
#| title: function

library(stringr)

reformat_to_snake_case <- function(char_vector) {
  # Convert to lowercase
  char_vector <- tolower(char_vector)
  
  # Replace special characters with spaces
  char_vector <- str_replace_all(char_vector, "[^a-z0-9]+", " ")
  
  # Trim leading and trailing whitespace
  char_vector <- str_trim(char_vector)
  
  # Replace spaces with underscores to create snake case
  char_vector <- str_replace_all(char_vector, " ", "_")
  
  return(char_vector)
}


# Function to downsample a Seurat object
downsample_seurat <- function(seurat_obj, cells_per_ident = 100) {
  # Get the cell names grouped by Ident
  cell_groups <- split(Cells(seurat_obj), Idents(seurat_obj))
  
  # Downsample each group to the desired number of cells
  sampled_cells <- unlist(lapply(cell_groups, function(cells) {
    if (length(cells) > 0) { # Check if the group is not empty
      if (length(cells) > cells_per_ident) {
        sample(cells, cells_per_ident) # Randomly sample if more cells exist
      } else {
        cells # Retain all cells if fewer than required
      }
    } else {
      NULL # Return NULL for empty groups
    }
  }))
  
  # Ensure there are no nulls and subset the Seurat object
  sampled_cells <- sampled_cells[!is.null(sampled_cells)]
  subset(seurat_obj, cells = sampled_cells)
}

```



```{r}
#| title: Load external data and reformat 

dat <- readRDS(here("..", "..",
                    "raw_data",
                    "wendisch_et_al_final",
                    "oli_github",
                    "Monocytes.Rds"))
dat$condition_clean <- reformat_to_snake_case(dat$condition )
Idents(dat) <- "condition_clean"

dat <- UpdateSeuratObject(dat)

dat <- dat |> group_by(condition_clean) |> slice_sample(n=params$downsample_size)

if (!file.exists(here("intermediate_data"))) dir.create("intermediate_data") #needed with snakemake?

file_name <- paste("small_seurat", as.character(params$downsample_size), "obj.rds", sep = "_")
write_rds(dat, here("intermediate_data",file_name))
```
