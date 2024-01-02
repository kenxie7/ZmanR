
<p align = "center">
    <img width="300px" src="./pic/TEST3.png" align="left" alt="ZmanR" />
</p>

<p>ZmanR is an R package for the assignment of time labels for temporally-resolved single-cell sequencing (Zman-seq - זמן). We utilized plate-based MARS-seq and FACS sorting to reveal the temporal dynamics of immune cells by measuring the time that immune cells enter the tumor. We also developed novel computational methods for time-informed trajectory analysis. 
</p> 


Overview
========
<p align = "center">
<img width="500px" src="./pic/time_labeling.png" alt="Overview" title="Overview" align="center">
</p>
The concept is to label all immune cells in the peripheral blood at different time points with different fluorescence markers and computationally determine the time bin that each immune cell belongs to. Since the cells could potentially obtain multiple time bins, we implement a generalized linear model on unstained cells to classify the stained cells for each individual fluorescence. Finally, we assign the time bin to a cell based on the 'most recent' time bin that the cell gets upon infiltration to the tumor. 
</p> 

### Table of content
- [Installation](#Installation)
- [Example Usage](#example-usage)
    - [Preprocessing](#--Preprocessing)
    - [TimeAssignment](#--Time-Assignment)
    - [TrajectoryAnalysis](#--Trajectory-Analysis)
- [Citation](#citation-&-references)
- [Maintenance](#Maintenance)
### Installation
	git clone https://github.com/kenxie7/ZmanR
	R
	require(devtools)
	devtools::install("ZmanR")
### Preprocessing
We provide the example data for the peripheral and treatment data. The treatment data will be updated as example data soon. </br>
### Time Assignment

    library(ZmanR)
    # First we check whether the fluorescence signal is normal according to unstained cells when we fit a GLM model.
    # For now, the input is a cell metadata dataframe with log-transformed fluorescence data as columns for each cell, we will update for using Seurat/sce object as input. 
    # This will output the Q-Q, fitted vs residual, etc plots for evaluation and setting the standard deviation parameter for timebin assignment.
    FACS_model_eval(well_fcs_mc_Blood,
    		    fluorophores = c("log_BV711.A", "log_BUV737.A", "log_BB515.A", "log_PE.A"))
    
<div align="center">
<img src="./pic/glm_fit.png" alt="glm_fit"/>
</div>

    # Visualize the distribution of FACS fluorophore for unstained vs stained cells
    plot_FACS_histogram(well_fcs_mc_Blood, 
    			fluorophores = c("log_BV711.A", "log_BUV737.A", "log_BB515.A", "log_PE.A"))
    
<div align="center">
<img src="./pic/fluo_hist.png" alt="fluo_hist"/>
</div>

    # Assign the timebin to cells, where sd_threshold corresponds to the standard deviation threshold for each fluorophore:
    well_fcs_mc_Blood = FACS_model(well_fcs_mc_Blood, sd_threshold = c(2.5, 1.5, 2.5, 2),
    				   fluorophores = c("log_BV711.A", "log_BUV737.A", "log_BB515.A", "log_PE.A"),
  				   timebins = c("12H", "24H", "36H", "48H"))

    # Evaluate the assigned timebin with stained and nonstained cells: plotting each timebin color with APC or other self-defined FACS channel.
    # The contour shows the actual classification by the GLM model
    plot_FACS_eval(well_fcs_mc_Blood, 
    		   fluorophores = c("log_BV711.A", "log_BUV737.A", "log_BB515.A", "log_PE.A"),
  		   timebins = c("12H", "24H", "36H", "48H"), y_axis = "log_APC.A")
    
<div align="center">
<img src="./pic/stain_asssignment.png" alt="stain_assignment"/>
</div>

### Trajectory Analysis
We use metacells to first unbiasedly assign cells to metacells and annotated the metacells by gene expression. Then using cells with the time label, we can derive the AUC time of each metacell and construct the Zman trajectory along time as shown schematically below.
<div align="center">
<img src="./pic/zman_schematic.jpg" alt="zman_schematic"/>
</div>


### Citation & References

If you liked or used our work, please cite:
###### Kirschenbaum, D., Xie, K., Ingelfinger, F., Katzenelenbogen, Y., Abadie, K., Look, T., ... & Amit, I. (2023). Time-resolved single-cell transcriptomics defines immune trajectories in glioblastoma. Cell. https://doi.org/10.1016/j.cell.2023.11.032


### Maintenance

If there's any questions / problems regarding ZmanR, please feel free to contact Ken Xie - kk.xie419@gmail.com. Thank you!

