# Experimental-Modal-Analysis

This repository contains the code and documentation for a comprehensive modal analysis project on a metallic plate. The primary objective was to characterize the structure's dynamic behavior by identifying its resonance frequencies, damping factors, and mode shapes. The project integrates Finite Element Method (FEM) modeling with experimental data analysis for robust cross-validation.

Key Features:

FEM Modeling: A finite element model of the plate, including its geometric dimensions and boundary conditions, to predict natural frequencies and mode shapes.

Experimental Data Acquisition: Data acquired through hammer impact tests on the structure. The data was pre-processed to calculate Frequency Response Functions (FRFs).

Modal Analysis (SDOF/MDOF): Implementation of both Single Degree of Freedom (SDOF) and Multiple Degree of Freedom (MDOF) modal analysis techniques to extract modal parameters from the experimental data.

Comparison and Validation: A detailed comparison between FEM and experimental results using the Modal Assurance Criterion (MAC) to evaluate the correlation between predicted and measured mode shapes.

Data Visualization: Graphs and plots to visualize FRFs, modal circles, and mode shapes, providing a clear understanding of the structure's dynamic behavior.

Tools Used:

MATLAB: For data pre-processing, modal analysis, and visualization.

FEM Software (e.g., Abaqus, ANSYS, etc.): For finite element modeling and analysis.

Who is this repository for:

This repository is ideal for students, researchers, and engineers interested in structural engineering, vibration analysis, experimental modal analysis, and sound engineering. It provides a complete workflow from theory to practice, demonstrating how to validate a numerical model with real-world data.
