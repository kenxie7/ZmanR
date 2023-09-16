
<p align = "center">
    <img width="300px" src="https://github.com/kenxie7/ZmanR/assets/16096351/d1aa50a1-000e-492a-bfef-052afba29924" align="left" alt="ZmanR" />
</p>

<p>ZmanR is an R package for the assignment of time labels for temporally-resolved single-cell sequencing (Zman-seq - זמן). We utilized plate-based MARS-seq and FACS sorting to reveal the temporal dynamics of immune cells by measuring the time that immune cells enter the tumor. We also developed novel computational methods for time-informed trajectory analysis. 
</p> 


Overview
========
<p align = "center">
<img width="500px" src="https://github.com/kenxie7/ZmanR/assets/16096351/34596fe4-925d-4a8b-94fe-2ef8a132781a" alt="Overview" title="Overview" align="center">
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
	require(devtools)
	devtools::install_github("kenxie7/ZmanR")
### Preprocessing
We provide the example data for the peripheral and treatment data. The treatment data will be updated as example data soon. </br>
The provided example data are cell metadata from the single cell datasets
### Time Assignment

    library(ZmanR)
    # First we check whether the fluorescence signal is normal according to unstained cells when we fit a GLM model.
    # This will output the Q-Q, fitted vs residual, etc plots for evaluation and setting the standard deviation parameter for timebin assignment.
    FACS_model_eval(well_fcs_mc_Blood)
    
<div align="center">
<img src="https://github.com/kenxie7/ZmanR1/assets/16096351/d0b4a6b9-5452-40bb-9ab2-b25fbbbb3175" alt="glm_fit"/>
</div>

    # Visualize the distribution of FACS fluorophore for unstained vs stained cells
    plot_FACS_histogram(well_fcs_mc_Blood)
    
<div align="center">
<img src="https://github.com/kenxie7/ZmanR1/assets/16096351/7c4950d1-3745-4657-b10b-32a913ca2a75" alt="fluo_hist"/>
</div>

    # Assign the timebin to cells, where sd_threshold corresponds to the standard deviation threshold for each fluorophore:
    well_fcs_mc_Blood = FACS_model(well_fcs_mc_Blood, sd_threshold = c(2.5, 1.5, 2.5, 2))

    # Evaluate the assigned timebin with stained and nonstained cells:
    # The contour shows the actual classification by the GLM model
    plot_FACS_eval(well_fcs_mc_Blood)
    
<div align="center">
<img src="https://github.com/kenxie7/ZmanR1/assets/16096351/a0d3417e-8645-48d7-8e28-bcb61a0aeda0" alt="stain_assignment"/>
</div>

### Trajectory Analysis
Functions to be uploaded for time-informed trajectory inference, visualizations of time-correlated genes and transcription factors, and ligand-receptor analysis.

### Citation & References

This work is currently under submission.

### Maintenance

If there's any questions / problems regarding ZmanR, please feel free to contact Ken Xie - kk.xie419@gmail.com. Thank you!

