
### Install and Load Seurat

	# if you have not installed the Seurat package on your computer, run this command
	# (only needs to be run once)
	install.packages("Seurat")

	# load Seurat package
	library(Seurat)

	# raw data downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124952 
	# tutorial data available at https://www.singlecellworkshop.com/workshop-raw-data.zip


### Import count matrix

	# read in the data directories containing the cell/barcode matrices
	pfc2.data <- Read10X(data.dir = "raw-data/pfc-sample2")
	pfc3.data <- Read10X(data.dir = "raw-data/pfc-sample3")
	pfc5.data <- Read10X(data.dir = "raw-data/pfc-sample5")
	pfc7.data <- Read10X(data.dir = "raw-data/pfc-sample7")


### Create Seurat object

	# create a new Seurat object for each sample
	# min.cells = 3, only genes detected in at least 3 cells will be included
	# min.features = 200, only cells with at least 200 genes detected will be included
	pfc2 <- CreateSeuratObject(counts = pfc2.data, project = "pfc-demo", min.cells = 3, min.features = 200)
	pfc3 <- CreateSeuratObject(counts = pfc3.data, project = "pfc-demo", min.cells = 3, min.features = 200)
	pfc5 <- CreateSeuratObject(counts = pfc5.data, project = "pfc-demo", min.cells = 3, min.features = 200)
	pfc7 <- CreateSeuratObject(counts = pfc7.data, project = "pfc-demo", min.cells = 3, min.features = 200)


### Inspect the new objects

	pfc2
	pfc3
	pfc5
	pfc7

	# remove the raw data to save processing space
	rm(pfc2.data, pfc3.data, pfc5.data, pfc7.data)


### Saving Seurat objects

	# save the Seurat object in .RDS format for easy access later
	# this way we will not need to re-run the CreateSeuratObject command before every analysis
	# give this file an informative name so that we know what its contents are
	saveRDS(pfc2, "pfc2-unprocessed-200g-seurat-object.RDS")
	saveRDS(pfc3, "pfc3-unprocessed-200g-seurat-object.RDS")
	saveRDS(pfc5, "pfc5-unprocessed-200g-seurat-object.RDS")
	saveRDS(pfc7, "pfc7-unprocessed-200g-seurat-object.RDS")

	# to reload this object from a saved file, use the following command
	# since we already have the data loaded for now, we do not need to use this command
	pfc2 <- readRDS("pfc2-unprocessed-200g-seurat-object.RDS")
	pfc3 <- readRDS("pfc3-unprocessed-200g-seurat-object.RDS")
	pfc5 <- readRDS("pfc5-unprocessed-200g-seurat-object.RDS")
	pfc7 <- readRDS("pfc7-unprocessed-200g-seurat-object.RDS")


### Inspecting object metadata

	# inspect the metadata for one of our objects using the 'head' function
	head(pfc2@meta.data)


### Adding custom metadata

	# to add custom metadata use the following command structure
	# you can replace 'sample_number' with any desired metadata slot name
	# this is useful for comparing/combining multiple samples or treatments
	pfc2@meta.data$sample_number <- "2" 
	pfc3@meta.data$sample_number <- "3" 
	pfc5@meta.data$sample_number <- "5" 
	pfc7@meta.data$sample_number <- "7" 


### Merging Seurat objects

	# merge multiple Seurat objects
	pfc <- merge(x = pfc2, y = list(pfc3, pfc5, pfc7))

	# remove individual objects to save processing space
	rm(pfc2, pfc3, pfc5, pfc7)

	# inpect our new combined object
	pfc


### Mitochondrial gene ratio

	# an important metadata slot to add in every experiment is the ratio of mitochondrial genes
	# detected in each cell - this can be used as a proxy for cell quality in most preparations
	pfc[["percent.mt"]] <- PercentageFeatureSet(object = pfc, pattern = "^mt-")

	# inspect the metadata for our object using the 'head' function to see our new metadata column(s)
	head(pfc@meta.data)


### Quality Control Filtering

	# next we will plot some QC metrics including nFeature_RNA, nCount_RNA, and percent.mt
	# examine these plots to see if any cutoffs should be considered
	VlnPlot(pfc, features = c("nFeature_RNA"), pt.size=0)
	VlnPlot(pfc, features = c("nCount_RNA"), pt.size=0)
	VlnPlot(pfc, features = c("percent.mt"), pt.size=0)

	# I like to begin with percent.mt, since it is usually the most effective at 
	# removing low quality samples
	# percent.mt cutoffs typically range from 5-10% depending on the sample
	# we can zoom into the lower range of percent.mt to help select a cutoff point
	VlnPlot(pfc, features = c("percent.mt"), pt.size=0, y.max=15)

	# 5% seems to be a reasonable cutoff, so let's subset our data to include
	# only cells with less than 5% percent.mt 
	pfc <- subset(pfc, subset = percent.mt < 5)

	# inspect the object again to see how many cells were removed
	pfc

	# next we will examine nFeature_RNA (AKA nGene)
	VlnPlot(pfc, features = c("nFeature_RNA"), pt.size=0, y.max=2000)
	VlnPlot(pfc, features = c("nFeature_RNA"), pt.size=0, y.max=1000)

	# we can trim these cells now, or cluster them and trim later
	# for now let's only include cells with >600 nFeature_RNA since those are likely to be
	# the highest quality
	pfc <- subset(pfc, subset = nFeature_RNA > 600)

	# inspect the object again to see how many cells were removed
	pfc

	# inpect our QC metrics again
	VlnPlot(pfc, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size=0)



### SCTransform normalization & scaling

	# we have trimmed low quality cells and are ready to proceed with analysis
	# first, scale and normalize our data using the SCTransform function
	# this may take several minutes to execute, and progress will display in the console
	pfc <- SCTransform(pfc)

	# this is another good point to save your data for easy reloading later
	saveRDS(pfc, "pfc-sctransformed.RDS")
	pfc <- readRDS("pfc-sctransformed.RDS")

### Dimensionality reduction with PCA

	pfc <- RunPCA(pfc, npcs = 60)


### Clustering 

	pfc <- FindNeighbors(pfc, dims = 1:60)
	pfc <- FindClusters(pfc, resolution = 0.7)


### UMAP non-linear dimensional reduction 

	pfc <- RunUMAP(pfc, dims = 1:60)


### Plotting clustered cells

	DimPlot(pfc, label=T)

	# to remove legend (quicker plotting)
	DimPlot(pfc, label=T) + NoLegend()

	# to group data by a metadata slot
	DimPlot(pfc, group.by="sample_number")

	# to split data by a metadata slot
	DimPlot(pfc, split.by="sample_number")

	# display the number of cells in each cluster
	table(pfc@meta.data$seurat_clusters)

	# display the number of cells in each cluster split by batch
	table(pfc@meta.data$seurat_clusters, pfc@meta.data$sample_number)


### Identify Marker Genes

	# check our default active assay (should be "SCT")
	DefaultAssay(pfc)

	# switch back to raw RNA count data in the "RNA" assay
	DefaultAssay(pfc) <- "RNA"

	# quickly normalize the raw RNA count data before running DGE analysis
	pfc <- NormalizeData(pfc)

	# find marker genes for cluster 11
	# only positively enriched genes, must be expressed in 5% of cells in cluster 11
	c11_markers <- FindMarkers(pfc, ident.1 = 11, only.pos = T, min.pct = 0.05)

	# examine the top 20 significant marker genes in cluster 11
	head(c11_markers, 20)

	# export the complete list of marker genes (for easy viewing outside of R)
	write.csv(c11_markers, "c11_markergenes.csv")



	# find positively associated marker genes for all clusters
	markers <- FindAllMarkers(pfc, only.pos=T)
	# export the complete marker gene list to a .csv file for easy viewing
	write.csv(markers, "pfc-allmarkers.csv")



### Visualizing Gene Expression

	# visualize some of the marker genes for cluster 11
	FeaturePlot(pfc, features = c("Gad1", "Gad2"), order=T)
	DotPlot(pfc, features = c("Gad1","Gad2","Slc32a1","Prox1","Pnoc"))
	VlnPlot(pfc, features = c("Gad1"))



### Subsetting cells for further analysis

	# select a subset of clusters ("neurons") for further analysis
	neurons <- subset(pfc, idents = c("0","1","4","7"))
	# examine the new subset
	neurons

	# switch to the SCT-normalized data before re-running PCA/clustering
	DefaultAssay(neurons) <- "SCT"

	# repeat PCA and clustering on the new subset
	neurons <- RunPCA(neurons, npcs = 30)
	neurons <- FindNeighbors(neurons, dims = 1:30)
	neurons <- FindClusters(neurons, resolution = 0.7)
	neurons <- RunUMAP(neurons, dims = 1:30)

	# remember to switch back to the RNA assay for visualization/gene expression
	DefaultAssay(neurons) <- "RNA"

	# plot the reclustered subset
	DimPlot(neurons, label=T)






### Challenges:

		#Use R help to look up additional parameters for the DotPlot() function
		
		#Try changing the dot.scale parameter - how do your plots change?

		#What is the default `min.pct` cutoff for the FindMarkers() function?

		#What happens to the output if you omit the `only.pos` parameter when running FindMarkers() on cluster 11?

		#Try including more or fewer PCs when clustering - how does that affect the results?

		#Try increasing/decreasing the resolution of the FindClusters() function and examine how clustering changes. 





