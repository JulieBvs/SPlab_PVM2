###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
#


#### General

GLOBAL_DESCRIPTION = "PVM2 - All_brain_cortex_non_neural_cells"

SCIENTIFIC_GROUP = "SPLAB"
SCIENTIFIC_PROJECT_NAME = "SPlab_PVM2"
EXPERIMENT_PROJECT_NAME = "All_brain_cortex_non_neural_cells"

#SAMPLE_NAME = ""
#SAMPLE_ID = ""


#### List of genes to be monitored
#SIGNATURE_MACROPHAGE = c( "Mrc1", "Cd68", "Fcgr3", "Lyve1", "H2-Aa", "Ptprc", "Cx3cr1", "Spi1", "Fcrls", "Siglech", "Sall1")
#CORE_MACRPHAGE_SIGN <- c("Ptplad2", "Tlr4", "Fcgr1", "Fgd4", "Cd14", "Tbxas1", "Fcgr3")

MARKERS_MICROGLIA = c("Ptprc", "Cx3cr1", "Emr1", "Aif1", "C1qa", "Cd68", "Fcgr1", "Spi1", "Itgam", "Siglech", "Trem2", "P2ry12", "Hexb") #"Ptprc" and "Adgre1" low
MARKERS_PVM = c("Lyve1", "Mrc1", "Cd163", "Ptprc", "Cx3cr1", "Emr1", "Aif1", "C1qa", "Cd68", "Fcgr1", "Spi1", "Csf1r") #Spi1 low
MARKERS_MHCII = c("H2-Aa")
MARKERS_FOS_COMPLEX = c("Atf3", "Fos", "Jun")
MARKERS_PROLIFERATION = c("Mki67", "Top2a", "Birc5")

MARKERS_ALL = c(MARKERS_PVM, MARKERS_MICROGLIA[(!MARKERS_MICROGLIA %in% MARKERS_PVM)], MARKERS_MHCII, MARKERS_PROLIFERATION)


#### Input / Output

# Output folder name in data folder (for R session object, lists of cells/genes) 
PATH_PROJECT = file.path( "your_path", 
                          SCIENTIFIC_GROUP,
                          "BIOINFO",
                          SCIENTIFIC_PROJECT_NAME,
                          EXPERIMENT_PROJECT_NAME)

PATH_PROJECT_RAWDATA = file.path( PATH_PROJECT, "00_RawData")
PATH_PROJECT_OUTPUT = file.path( PATH_PROJECT, "05_Output")

# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_",     
                            EXPERIMENT_PROJECT_NAME, "_"
                            #startTimeFileName, "_",
                            #paramsHash, "_"
)


#### Color palets
col_comet = c("1"="#E6E7E8","2"="#3A97FF","3"="#8816A7","4"="black")
col_greyMagma = c("2"="grey", "4"="#FB8861FF", "5"="#B63679FF", "3"="#51127CFF", "1"="#000004FF")
col_rainbow = c("#E96767", "#EC9342", "#FFCA3A", "#8AC926", "#4DADE8", "#9579B9", "#E25CB4",
                "#DB2020", "#DA7316", "#F0AE00", "#6D9D1E", "#1882C0", "#71529A", "#D02494",
                "#9F1717", "#AE5B11", "#C48D00", "#517416", "#115C8A", "#584178", "#9D1C70",
                "#EF9292", "#F2B57D", "#FFDA77", "#B6E36A", "#7BC4EE", "#AD98C9", "#EA8AC9")
#barplot(rep(1,28), col=col_rainbow)


#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



