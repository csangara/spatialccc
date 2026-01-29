methods_cheng <- c("CellAgentChat",
             "GCNG", "NCEM", "spaCI", "CLARIFY", "DeepCOLOR", "DeepTALK", "CLARA", "PearlST",
             "CellNEST", "NichCompass", "GraphComm", "HoloNet", "CellGAT", "SpaGraphCCI",
             "SpaCCC", "DeepLinc", "GITIII", "OrgaCCC", "STCase", "AMICI", "STEAMBOAT",
             "SPROUT", "SpaCeNet", "LRnetST",
             "cell2cell", "DcjComm", "SpaCET", "SVCA", "Giotto", "BATCOM",
             "stLearn", "CytoSignal", "CellPhoneDB", "BulkSignalR", "IGAN",
             "Niche-DE", "Spacia", "CellChat", "SpatialDM", "NICHES",
             "scCellFie", "TWCOM", "Squidpy", "CSOmap", "MESSI",
             "CCPLS", "CellNeighborEX", "CPPLS-MLP")
methods_cheng %>% sort
methods_armingol <- c("STcomm", "COMMOT", "SpaTalk", "Renoir",
              "stMLnet", "MISTy", "SpaOTsc", "Neighbor-seq", "RNA-Magnet")

intersect(tolower(methods_cheng), tolower(methods_armingol))

all_methods <- c(methods_cheng, methods_armingol)
length(all_methods)
