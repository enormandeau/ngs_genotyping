# Information about the program
cat("\n")
cat("HIGGYPOP Version 1.6","\n")
cat("Haplotype Inference and GenotypinG Yn POPulation samples with individual barcodes","\n")

cat("HIGGYPOP depends on the following packages: gee, ape, sequinr","\n")

# Loads R packages gee et ape for sequence analysis
library(gee, warn.conflicts = FALSE)
library(ape, warn.conflicts = FALSE)
library(seqinr, warn.conflicts = FALSE)

divide <- function(d,Th,Ac)
# The "divide" function splits the dendrogram obtained from the distance matrix at its deepest 
# node and returns a list of two distance matrices corresponding to the two subtrees obtained.
# If a subtree contains a proportion of the total number of individual reads which is below the 
# threshold value "Th", a 0 is stored in the list instead of the distance matrix. If the deepest
# branch of the dendrogram to split is shorter than threshold "Ac", the original matrix is kept.
	{
	# "hcw" stores the result of hierarchical clustering based on Ward's minimum variance method
	hcw <- hclust(d,method="ward")
	# A cluster membership is assigned to each sequence
	memb <- cutree(hcw,k=2)
	L <- list(names(memb[memb==1]),names(memb[memb==2]))
	# The two distance matrices in the output are first defined as zero
	d1 <- 0
	d2 <- 0
	# If the dendrogram to be divided only consists of closely related sequences, then the
	# original distance matrix is kept if it contains a sufficent proportion of the reads.
	# This corresponds in the dendrogram to a branch that contains a collection of reads
	# from a single allele, plus eventual reads containing a few seqeuncing errors.
	if ((length(d)/N) >= Th & max(hcw$height) <= Ac)
		{d1 <- d}
	# If the dendrogram is divided into two sub-trees, the distance matrix of the first
	# sub-tree is stored in "d1" if it contains a sufficent proportion of related reads.
	if (((length(L[[1]]))/N) >= Th & max(hcw$height) > Ac)
		{d1 <- as.dist((Matrix_dist)[L[[1]],L[[1]]], diag = FALSE, upper = FALSE)}
	# Same thing for the distance matrix of the second sub-tree, stored in "d2"
	if (((length(L[[2]]))/N) >= Th & max(hcw$height) > Ac)
		{d2 <- as.dist((Matrix_dist)[L[[2]],L[[2]]], diag = FALSE, upper = FALSE)}
	# The outpout is a list of two distance matrices: either the original distance matrix
	# and the zero matrix, the original matrix without its deepest branch containing
	# errors and the zero matrix, or two non-null matrices obtained by spitting the input
	# distance matrix.
	list(d1,d2)
	}

Reduc <- function(A)
# The "Reduc" function takes a list containg non-zero distance matrices
# and zero distance matrices and returns a list from which zero distance
# matrices have been removed. This function cleans the output from the
# "divide" function in order to save computation time.
	{
	Clst = NULL
	for (j in 1:length(A))
		{
		if (length(A[[j]])!=1)
		Clst <- c(Clst,list(A[[j]]))
		}
	if (is.null(Clst)) {Clst=a}
	Clst
	}

ColLab <- function(DEND)
# The "ColLab" function colors the leaves of an individual dendrogram to
# show the leaves that were assigned to clusters and those that were
# rejected by the algorithm.
	{
	if(is.leaf(DEND))
		{
		At <- attributes(DEND)
		# Starts by setting the grey color to each leaf
		attr(DEND, "nodePar") <- c(At$nodePar, list(lab.col = "transparent", lab.cex=0.4, col="grey", pch=16, cex=0.5))
		# And then gives a different color to the leaves that belong to each cluster
		Paltette = rainbow(Nclt)
		for (p in 1:Nclt)
			{
			if (length(which(COLOR[[p]] == At$label)) == 1)
				{
				attr(DEND, "nodePar") <- c(At$nodePar, list(lab.col = "transparent", lab.cex=0.4, col=Paltette[p], pch=16, cex=0.5))
				}
			}
		}
	DEND
	}

# Asks for choosing between Phase 1 and phase 2
cat("\n")
cat("Enter 1 to perform phase 1: allele detection from individual reads\n")
cat("Enter 2 to perform phase 2: individual genotypes detection from phase 1 output)\n")
cat("(1/2)\n")
Phase <- readLines(n=1)
if(Phase==1)
	{
	# Prompts the user to specify the folder to upload and choose parameter values for the analysis
	cat("Phase 1: detecting alleles from individual reads","\n")
	# Prompts for the full path of the folder where individual fasta files have been stored
	cat("\n")
	cat("Enter the full path of the folder containing individual fasta files","\n")
	Path <- readLines(n=1)
	setwd(Path)

	# Asks the user to check if the number of individual fasta files detected is correct
	List_fasta <- list.files(Path)
	cat(paste("Number of individual fasta files = ",length(List_fasta),"\n"))
	cat("Is this correct?","\n")
	cat("(y/n)","\n")
	answer1 <- readLines(n=1)
	if(answer1=="n")cat("Check the folder path","\n")
	if(answer1!="y"&answer1!="n")cat("Enter y or n","\n")
	if(answer1=="y")
	{cat("Starting phase 1: Identification of individual sequence clusters","\n")

	# Prompts for the internal branch length parameter value. A value of 0.06 will be OK in most of
	# the cases to distinguish alleles from sequencing errors. A lower value will result in a higher
	# rate of exclusion of reads containing sequencing errors, and an increased power to distinguish
	# closely related alleles that differ from a single base mutation.
	cat("Choose internal branch length parameter for sequence cluster discrimination","\n")
	cat("(e.g. 0.06)","\n")
	Ac <- as.numeric(readLines(n=1))
	}

	# Asks the user to choose between an automated and a manual version of phase 1. In the automated
	# version, the entered Ac parameter value is used for each individual, the number of clusters is
	# automatically detected and each cluster's consensus sequence is exported for each individual.
	# In the manual version, an automatic detection of cluster is also performed in a first step,
	# but the user may decide to change the internal branch length parameter value and the minimal
	# amount of sequences in each cluster after looking at the results.
	cat("Check clustering for each individual before exporting consensus sequences?","\n")
	cat("Enter y if you want to have a full control for each individual","\n")
	cat("Enter n if you want to perform automatic clustering (large datasets)","\n")
	cat("(y/n)","\n")
	answer2 <- readLines(n=1)
	if(answer2=="n")
		{
		cat("Starting automatic clustering and exporting individual consensus sequences","\n","\n")
		# For each individual having a 'fasta' file in the folder, the number of sequence clusters is
		# determined and stored in NCLUST, and added to LCLUST which stores the number of clusters of
		# each indivudual in the order the are called.
		LCLUST <- NULL
		cat("Enter minimal proportion of sequences to define a cluster","\n")
		sequ_min <- as.numeric(readLines(n=1))
		# cat("Enter MAXIMAL proportion of sequences to define a cluster","\n")
		sequ_max <- sequ_min
		sequ = seq(sequ_min,sequ_max,0.01)
		for (l in 1:(length(List_fasta)))
			{
			# Loads the fasta alignment of individual "l"
			input.aln <- read.dna(file=List_fasta[l], format = "fasta")
			cat("Now processing file", List_fasta[l], "\n")
			# Computes all pairwise raw distances and places them into a distance matrix
			Matrix_dist <- (dist.dna(input.aln, model = "raw", as.matrix = TRUE, pairwise.deletion = FALSE))
			# Makes a lower triangular distance matrix of class 'dist' # 
			d <- as.dist(Matrix_dist, diag = FALSE, upper = FALSE)
			# Performs one step of hierarchical clustering based on Ward's minimum variance method
			HCW <- hclust(d,method="ward")
			DEND <- as.dendrogram(HCW)
			NCLUST <- NULL
			# Explores the number of retained clusters as a function of the amount of sequence data in each cluster
			for (Thr1 in sequ)
				{
				# Thr1 <- seq(0.01,0.5,0.01)[k]
				N <- length(colnames(Matrix_dist))
				W <- list(d,divide(d,Thr1,Ac))
				S <- length(unlist(W[[1]]))
				while (length(unlist(W[[2]])) < S)
					{
					a <- NULL
					for (i in 1:length(W[[length(W)]]))
						{
						if(length(W[[length(W)]][[i]])!=1)
							{w <- divide(W[[length(W)]][[i]],Thr1,Ac)}
						if(length(W[[length(W)]][[i]])==1)
							{w <- list(0,0)}
						a <- c(a,w)
						}
					S <- length(unlist(W[[2]]))
					a <- Reduc(a)
					W <- list(W,a)
					}
				NCLUST <- c(NCLUST,(length(Reduc(W[[2]]))))
				}
			LCLUST <- c(LCLUST,rle(NCLUST)[[2]][rle(NCLUST)[[1]]==max(rle(NCLUST)[[1]])][1])
			Thr2 <- max(sequ[NCLUST==LCLUST[length(LCLUST)]])
			W <- list(d,divide(d,Thr2,Ac))
			S <- length(unlist(W[[1]]))
			while (length(unlist(W[[2]])) < S)
				{
				a <- NULL
				for (i in 1:length(W[[length(W)]]))
					{
					if(length(W[[length(W)]][[i]])!=1)
						{w <- divide(W[[length(W)]][[i]],Thr2,Ac)}
					if(length(W[[length(W)]][[i]])==1)
						{w <- list(0,0)}
					a <- c(a,w)
					}
				S <- length(unlist(W[[2]]))
				a <- Reduc(a)
				W <- list(W,a)
				}
			input.aln2 <- read.alignment(file=List_fasta[l], format = "fasta")
			FASL <- NULL
			COLOR <- NULL
			Nclt <- LCLUST[length(LCLUST)]
			Prefix <- unique(gsub(pattern = "_.+",replacement = "", x=input.aln2$nam))
			Seqnames <- NULL
			for (n in 1:Nclt)
				{
				r <- rep(FALSE,dim(input.aln)[1])
				for (m in 1:length(names(as.matrix((W[[2]])[[n]])[1,])))
					{
					r <- r+(rownames(input.aln)==names(as.matrix((W[[2]])[[n]])[1,])[m])
					}
				FASL <-c(FASL,list(toupper(c2s(consensus((as.matrix(input.aln2))[r==1,])))))
				COLOR <- c(COLOR, list(as.factor(row.names((as.matrix(input.aln2))[r==1,]))))
				Seqnames <- c(Seqnames, paste(Prefix,"_",n,sep = ""))
				}
			dend_colored <- dendrapply(DEND,ColLab)
			plot(dend_colored, xlab=paste("Automatically generated dendrogram of individual",Prefix))
			cat(paste("  -> Detecting",Nclt,"clusters with more than", Thr2 * 100, "percents\n"))
			write.fasta(sequences=FASL,names=Seqnames, nbchar = 60, file.out=paste(Prefix,".fas"), open = "a")
			}
		}

	if(answer2=="y")
		{
		cat("Starting clustering with individual check option","\n","\n")
		LCLUST <- NULL
		for (l in 1:(length(List_fasta)))
			{
			input.aln <- read.dna(file=List_fasta[l], format = "fasta")	
			Matrix_dist <- dist.dna(input.aln, model = "raw", as.matrix = TRUE, pairwise.deletion = FALSE)
			d <- as.dist(Matrix_dist, diag = FALSE, upper = FALSE)
			HCW <- hclust(d,method="ward")
			DEND <- as.dendrogram(HCW)
			NCLUST <- NULL
			for (k in 7:13)
				{
				Thr1 <- seq(0.01,0.5,0.01)[k]
				N <- length(colnames(Matrix_dist))
				W <- list(d,divide(d,Thr1,Ac))
				S <- length(unlist(W[[1]]))
				while (length(unlist(W[[2]])) < S)
					{
					a <- NULL
					for (i in 1:length(W[[length(W)]]))
						{
						if(length(W[[length(W)]][[i]])!=1)
							{w <- divide(W[[length(W)]][[i]],Thr1,Ac)}
						if(length(W[[length(W)]][[i]])==1)
							{w <- list(0,0)}
						a <- c(a,w)
						}
					S <- length(unlist(W[[2]]))
					a <- Reduc(a)
					W <- list(W,a)
					}
				NCLUST <- c(NCLUST,(length(Reduc(W[[2]]))))
				}
			LCLUST <- c(LCLUST,rle(NCLUST)[[2]][rle(NCLUST)[[1]]==max(rle(NCLUST)[[1]])][1])
			Thr2 <- max(seq(0.07,0.13,0.01)[NCLUST==LCLUST[length(LCLUST)]])
			W <- list(d,divide(d,Thr2,Ac))
			S <- length(unlist(W[[1]]))
			while (length(unlist(W[[2]])) < S)
				{
				a <- NULL
				for (i in 1:length(W[[length(W)]]))
					{
					if(length(W[[length(W)]][[i]])!=1)
						{w <- divide(W[[length(W)]][[i]],Thr2,Ac)}
					if(length(W[[length(W)]][[i]])==1)
						{w <- list(0,0)}
					a <- c(a,w)
					}
				S <- length(unlist(W[[2]]))
				a <- Reduc(a)
				W <- list(W,a)
				}
			input.aln2 <- read.alignment(file=List_fasta[l], format = "fasta")
			FASL <- NULL
			COLOR <- NULL
			Nclt <- LCLUST[length(LCLUST)]
			Prefix <- unique(gsub(pattern = "_.+",replacement = "", x=input.aln2$nam))
			Seqnames <- NULL
			for (n in 1:Nclt)
				{
				r <- rep(FALSE,dim(input.aln)[1])
				for (m in 1:length(names(as.matrix((W[[2]])[[n]])[1,])))
					{
					r <- r+(rownames(input.aln)==names(as.matrix((W[[2]])[[n]])[1,])[m])
					}
				FASL <-c(FASL,list(toupper(c2s(consensus((as.matrix(input.aln2))[r==1,])))))
				COLOR <- c(COLOR, list(as.factor(row.names((as.matrix(input.aln2))[r==1,]))))
				Seqnames <- c(Seqnames, paste(Prefix,"_",n,sep = ""))
				}
			dend_colored <- dendrapply(DEND,ColLab)
			plot(dend_colored, xlab=paste("Automatically generated dendrogram of individual",Prefix))
			cat(paste("Detecting",Nclt,"clusters in individual",Prefix,"\n"))
			# Prompts to check the graphical output
			cat("Accept clustering?","\n")
			cat("(y/n)","\n")
			answer3 <- readLines(n=1)
			while (answer3=="n")
				{
				cat("Enter new internal branch length parameter value","\n")
				cat("(e.g. 0.06)","\n")
				Ac2 <- as.numeric(readLines(n=1))
				cat("Enter new minimal proportion of sequences per cluster","\n")
				cat("(e.g. 0.10)","\n")
				Thr3 <- as.numeric(readLines(n=1))
				W <- list(d,divide(d,Thr3,Ac2))
				S <- length(unlist(W[[1]]))
				while (length(unlist(W[[2]])) < S)
					{
					a <- NULL
					for (i in 1:length(W[[length(W)]]))
						{
						if(length(W[[length(W)]][[i]])!=1)
							{w <- divide(W[[length(W)]][[i]],Thr3,Ac2)}
						if(length(W[[length(W)]][[i]])==1)
							{w <- list(0,0)}
						a <- c(a,w)
						}
					S <- length(unlist(W[[2]]))
					a <- Reduc(a)
					W <- list(W,a)
					}
				Nclt <- length(Reduc(W[[2]]))
				input.aln2 <- read.alignment(file=List_fasta[l], format = "fasta")
				FASL <- NULL
				COLOR <- NULL
				Prefix <- unique(gsub(pattern = "_.+",replacement = "", x=input.aln2$nam))
				Seqnames <- NULL
				for (n in 1:Nclt)
					{
					r <- rep(FALSE,dim(input.aln)[1])
					for (m in 1:length(names(as.matrix((W[[2]])[[n]])[1,])))
						{
						r <- r+(rownames(input.aln)==names(as.matrix((W[[2]])[[n]])[1,])[m])
						}
					FASL <-c(FASL,list(toupper(c2s(consensus((as.matrix(input.aln2))[r==1,])))))
					COLOR <- c(COLOR, list(as.factor(row.names((as.matrix(input.aln2))[r==1,]))))
					Seqnames <- c(Seqnames, paste(Prefix,"_",n,sep = ""))
					}
				dend_colored <- dendrapply(DEND,ColLab)
				plot(dend_colored, xlab=paste("New dendrogram of individual",Prefix))
				cat(paste("Now detecting",Nclt,"clusters in individual",Prefix,"\n"))
				# Prompts to check the graphical output
				cat("Accept clustering?","\n")
				cat("(y/n)","\n")
				answer3 <- readLines(n=1)
				}
			cat(paste("Consensus sequences of individual",Prefix,"have been exported","\n","\n"))
			write.fasta(sequences=FASL,names=Seqnames, nbchar = 60, file.out=paste(Prefix,".fas"), open = "a")
			}
		}
	}

if(Phase==2)
	{
	# Prompts the user to specify the path of the file containing an alignment of the consensus sequences from phase 1
	cat("Phase 2: detecting individual genotypes","\n")
	# Prompts for the full path of the folder where the fasta alignment has been stored
	cat("\n")
	cat("Enter the full path of the folder containing the fasta alignment","\n")
	Path <- readLines(n=1)
	setwd(Path)
	cat("Enter the full name of the file to upload (e.g. toto.fas)","\n")
	Filename <- readLines(n=1)
	input.aln <- read.dna(file=Filename, format = "fasta")
	Matrix_dist <- dist.dna(input.aln, model = "raw", as.matrix = TRUE, pairwise.deletion = FALSE)
	d <- as.dist(Matrix_dist, diag = FALSE, upper = FALSE)
	NI <- unique(gsub(pattern = "_.+",replacement = "", x=colnames(Matrix_dist)))
	cat(paste("Detecting",length(NI),"individuals in the alignment:"),"\n")
	print(NI)
	# Performs one step of hierarchical clustering based on Ward's minimum variance method
	HCW <- hclust(d,method="ward")
	DEND <- as.dendrogram(HCW)
	# Prompts for internal branch length parameter and minimal allele counts
	cat("Enter the internal branch length parameter value","\n")
	cat("(e.g. 0.06)","\n")
	Ac <- as.numeric(readLines(n=1))
	cat("Enter the minimal number of individuals sharing a given allele for this allele to be scored","\n")
	cat("(e.g. 2)","\n")
	NMIN <- as.numeric(readLines(n=1))
	N <- length(colnames(Matrix_dist))
	Thr <- NMIN/N
	W <- list(d,divide(d,Thr,Ac))
	S <- length(unlist(W[[1]]))
	while (length(unlist(W[[2]])) < S)
		{
		a <- NULL
		for (i in 1:length(W[[length(W)]]))
			{
			if(length(W[[length(W)]][[i]])!=1)
				{w <- divide(W[[length(W)]][[i]],Thr,Ac)}
			if(length(W[[length(W)]][[i]])==1)
				{w <- list(0,0)}
			a <- c(a,w)
			}
		S <- length(unlist(W[[2]]))
		a <- Reduc(a)
		W <- list(W,a)
		}
	Nclt <- length(Reduc(W[[2]]))
	cat(paste("Detecting",Nclt,"alleles in",length(NI),"individuals"),"\n")
	input.aln2 <- read.alignment(file=Filename, format = "fasta")
		FASL <- NULL
		COLOR <- NULL
		Genotypes <- NULL
		Seqnames <- NULL
		for (n in 1:Nclt)
			{
			r <- rep(FALSE,dim(input.aln)[1])
			for (m in 1:length(names(as.matrix((W[[2]])[[n]])[1,])))
				{
				r <- r+(rownames(input.aln)==names(as.matrix((W[[2]])[[n]])[1,])[m])
				}
			FASL <-c(FASL,list(toupper(c2s(consensus((as.matrix(input.aln2))[r==1,])))))
			COLOR <- c(COLOR, list(as.factor(row.names((as.matrix(input.aln2))[r==1,]))))
			Genotypes <- c(Genotypes, paste(COLOR[[n]],"_A",n,sep = ""))
			Seqnames <- c(Seqnames, paste("A_",n,sep = ""))
			}
		dend_colored <- dendrapply(DEND,ColLab)
		plot(dend_colored, xlab=paste("dendrogram of",Nclt,"retained alleles, each found in at least",NMIN,"individuals"))
	# Prompts to check the graphical output
	cat("Accept clustering?","\n")
	cat("(y/n)","\n")
	answer4 <- readLines(n=1)
	while (answer4=="n")
		{
		cat("Enter new internal branch length parameter value","\n")
		cat("(e.g. 0.06)","\n")
		Ac2 <- as.numeric(readLines(n=1))
		cat("Enter the new minimal number of individuals sharing a given allele for this allele to be scored","\n")
		cat("(e.g. 2)","\n")
		NMIN2 <- as.numeric(readLines(n=1))
		Thr2 <- NMIN2/N
		W <- list(d,divide(d,Thr2,Ac2))
		S <- length(unlist(W[[1]]))
		while (length(unlist(W[[2]])) < S)
			{
			a <- NULL
			for (i in 1:length(W[[length(W)]]))
				{
				if(length(W[[length(W)]][[i]])!=1)
					{w <- divide(W[[length(W)]][[i]],Thr2,Ac2)}
				if(length(W[[length(W)]][[i]])==1)
					{w <- list(0,0)}
				a <- c(a,w)
				}
			S <- length(unlist(W[[2]]))
			a <- Reduc(a)
			W <- list(W,a)
			}
		Nclt <- length(Reduc(W[[2]]))
		FASL <- NULL
		COLOR <- NULL
		Genotypes <- NULL
		Seqnames <- NULL
		for (n in 1:Nclt)
			{
			r <- rep(FALSE,dim(input.aln)[1])
			for (m in 1:length(names(as.matrix((W[[2]])[[n]])[1,])))
				{
				r <- r+(rownames(input.aln)==names(as.matrix((W[[2]])[[n]])[1,])[m])
				}
			FASL <-c(FASL,list(toupper(c2s(consensus((as.matrix(input.aln2))[r==1,])))))
			COLOR <- c(COLOR, list(as.factor(row.names((as.matrix(input.aln2))[r==1,]))))
			Genotypes <- c(Genotypes, paste(COLOR[[n]],"_A",n,sep = ""))
			Seqnames <- c(Seqnames, paste("A_",n,sep = ""))
			}
		dend_colored <- dendrapply(DEND,ColLab)
		plot(dend_colored, xlab=paste("dendrogram of",Nclt,"retained alleles, each found in at least",NMIN2,"individuals"))
		cat(paste("Now detecting",Nclt,"alleles in",length(NI),"individuals"),"\n")
		# Prompts to check the graphical output
		cat("Accept clustering?","\n")
		cat("(y/n)","\n")
		answer4 <- readLines(n=1)
		}
	print(data.frame(Genotypes))
	output_file = paste("genotypes_B", as.character(Ac), "_N", as.character(NMIN), "_", format(Sys.time(), "%Y-%m-%d"), ".txt")
	output_file = gsub(" ", "", output_file)
	# write.table(data.frame(Genotypes),output_file,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="")
	write.fasta(sequences=FASL,names=Seqnames, nbchar = 60, file.out="allele_database.fasta", open = "a")
	}

