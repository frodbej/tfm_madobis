# 🔬💻 TFM MADOBIS - An image analysis pipeline for detecting different cell phenotypes in breast cancer cells 🔬💻
This is the readme file that contains the information about the code and data used in the TFM that are available in this repository.

The processing pipeline developed in CellProfiler for this project is available in cppipe format. This pipeline has been used to generate the output data used in the analysis from image raw data.

The code and data corresponding to gene network analysis are available in the folder "gene_network_analysis".

The code and a sample of per_image and per_object data due to the large size of the files are available in the folder "image_data_analysis". Complete raw data can be accessed upon request to the authors.

The code and data corresponding to the analysis of the classifier results are available in the folder "classifier_results_analysis".

## Description of the project

Cell shape provides valuable information about cellular phenotype and the cell’s physiological state due to the connection between morphology and phenotype. Thus, morphological analysis could allow us to understand relevant underlying mechanisms in diseases with an important morphological component such as cancer, particularly in the metastatic process.

Furthermore, advances in automated microscopy have made it possible to develop image-based high-throughput cell profiling assays, that allow defining complex cellular phenotypes through the extraction of multiple features for each cell of the population. Thus, with this approach, we can generate large quantities of data at single-cell resolution to identify genes involved in a biological process and heterogeneous cell behaviors.

Based on the importance of cell shape dynamics on metastasis, we used a small interfering RNA library to specifically suppress more than 500 proteins related to the cytoskeleton on triple-negative breast cancer cells. Automated fluorescence images were taken providing a large amount of image data that was processed using CellProfiler software to extract relevant cell parameters from thousands of images with hundreds of cells per image.

Due to the high dimensionality and complexity of the data, machine learning algorithms were used to identify discordant phenotypes indicating new potential targets or effectors and concordant phenotypes showing a relation between different genes.                                                                                                                                  
## Contact

ferrodbej@alum.us.es
