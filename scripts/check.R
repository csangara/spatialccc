test <- R.utils::commandArgs(trailingOnly = TRUE)

#print(test)
#print(test[1])

file_name <- stringr::str_split(basename(test[1]), "\\.")[[1]][1]
dir_name <- dirname(test[1])
print(dir_name)
ext <- stringr::str_split(basename(test[1]), "\\.")[[1]][2]

print(file_name)
print(ext)
