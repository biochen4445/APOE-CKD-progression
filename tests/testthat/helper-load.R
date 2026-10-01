suppressPackageStartupMessages(library(data.table))
for (f in list.files(file.path("..", "..", "R", "functions"), pattern = "\\.R$", full.names = TRUE)) source(f)
