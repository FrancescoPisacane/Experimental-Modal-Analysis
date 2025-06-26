# Experimental-Modal-Analysis

## A Deep Dive into Structural Dynamics: Bridging FEM and Experimental Data

This repository showcases a comprehensive modal analysis project on a metallic plate, a cornerstone of structural dynamics. My goal was to fully characterize the plate's dynamic behavior by identifying its natural frequencies, damping factors, and intricate mode shapes. This project goes beyond basic analysis by robustly integrating Finite Element Method (FEM) simulations with experimental impact hammer test data, providing a complete and validated workflow from theory to practice.

## Key Features and Methodologies

- **FEM Modeling**: A detailed finite element model of the plate was created based on its precise geometric dimensions and boundary conditions. This numerical model served as a predictive tool, giving me a baseline for the structure's natural frequencies and mode shapes.
- **Experimental Data Acquisition**: Using an instrumented hammer, I excited the plate at various points, capturing the structural response with accelerometers. The raw time-domain data was then meticulously processed to compute **Frequency Response Functions (FRFs)**, the cornerstone of my experimental analysis.
- **Advanced Modal Analysis Techniques**: This project implements several powerful modal analysis techniques to extract the modal parameters from the measured FRFs:
    - **Single Degree of Freedom (SDOF) Analysis**: I used this method to analyze individual modes by fitting circles to the FRF plots near resonance peaks, allowing for the quick extraction of frequency and damping for isolated modes.
    - **Multiple Degree of Freedom (MDOF) Analysis**: To address closely coupled modes and provide more robust results, I employed advanced MDOF methods, including:
        - **Ibrahim's Time Domain (ITD) Method**: A time-domain technique that identifies modal parameters from the free decay response.
        - **Prony's Method**: Another time-domain method used to model a signal as a sum of decaying sinusoids, directly yielding frequencies, damping ratios, and amplitudes.
- **Validation and Correlation**: The project culminates in a rigorous validation phase. I performed a detailed comparison between the mode shapes predicted by my FEM model and those extracted from the experimental data using the **Modal Assurance Criterion (MAC)**. The MAC matrix provided a powerful visual and numerical tool to quantify the correlation between the two sets of results, confirming the accuracy of my models and measurements.
- **Data Visualization**: The repository includes scripts to generate a variety of plots, including FRFs, Nyquist plots (modal circles), and animations of the identified mode shapes, making the complex dynamic behavior of the plate easy to understand.

## Tools Used

- **MATLAB**: All data processing, analysis, and visualization were performed using MATLAB, showcasing its powerful capabilities for system identification and signal processing.
- **FEM Software (Abaqus)**: Used for the initial finite element modeling of the plate.

## Who is this repository for

This repository is an invaluable resource for students, researchers, and engineers in mechanical, civil, and aerospace engineering. It provides a complete, hands-on workflow for anyone looking to learn or apply experimental modal analysis, structural dynamics, and system identification.
