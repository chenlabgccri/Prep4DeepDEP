library(Prep4DeepDEP)

# use example data of 3 cell lines
path <- system.file("examples/",package = "Prep4DeepDEP")

# example of exp.data
exp.data <- read.delim( paste0(path,"exp.file.txt"),sep = "\t",
                        header = TRUE,stringsAsFactors = FALSE)
head(exp.data)
dim(exp.data)

# example of mut.data
mut.data <- read.delim( paste0(path,"mut.file.txt"),sep = "\t",
                        header = TRUE,stringsAsFactors = FALSE)
head(mut.data)
dim(mut.data)

# example meth.data
meth.data <- read.delim( paste0(path,"meth.file.txt"),sep = "\t",
                        header = TRUE,stringsAsFactors = FALSE)
head(meth.data)
dim(meth.data)

# example cna data
cna.data <- read.delim( paste0(path,"cna.file.txt"),sep = "\t",
                        header = TRUE,stringsAsFactors = FALSE)
head(cna.data)
dim(cna.data)

# example dep.data for prediction (with a list of DepOIs only) and
# training (DepOIs with CRISPR gene dependeny scores)
dep.data.p <- read.delim( paste0(path,"dep.file.prediction.txt"),sep = "\t",
                        header = TRUE,stringsAsFactors = FALSE)
head(dep.data.p)
dim(dep.data.p)

dep.data.t <- read.delim( paste0(path,"dep.file.training.txt"),sep = "\t",
                        header = TRUE,stringsAsFactors = FALSE)
head(dep.data.t)
dim(dep.data.t)

# run Prep4DeepDEP "training" mode with single omics
Prep4DeepDEP(exp.data = exp.data, dep.data = dep.data.t, mode = "training",
             filename.out = "./test1")

# run Prep4DeepDEP "training" mode with full four omics
Prep4DeepDEP(exp.data = exp.data,
             mut.data = mut.data,
             meth.data = meth.data,
             cna.data = cna.data,
             dep.data = dep.data.t,
             mode = "training",
             filename.out = "./test2")

# run Prep4DeepDEP "prediction" mode with user's 500 DepOIs
Prep4DeepDEP(exp.data = exp.data,
             mut.data = mut.data,
             meth.data = meth.data,
             cna.data = cna.data,
             dep.data = dep.data.p,
             mode = "prediction",
             filename.out = "./test3")

# run Prep4DeepDEP "prediction" mode with default 1298 DepOIs
Prep4DeepDEP(exp.data = exp.data,
             mut.data = mut.data,
             meth.data = meth.data,
             cna.data = cna.data,
             dep.data = NULL,
             mode = "prediction",
             filename.out = "./test4")

