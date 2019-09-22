# cryoEM-notebooks

This repo contains information that we deem useful in understanding the end-to-end computational problems of cryoEM and cryoET. Specifically, we develop jupyter notebooks that walks us through the mathematical and computational operations required to reconstruct density maps from micrographs. Work in progress...

## Image formation model 
*see notebook on [(github)](notebooks/ImageFormationModel/Image%20Formation%20Model.ipynb) [(colab)](https://colab.research.google.com/github/slaclab/cryoEM-notebooks/blob/master/notebooks/ImageFormationModel/Image%20Formation%20Model.ipynb)*
Here we provide details on the microscope architecture, useful to understand the image formation model and its theoretical limits. We then proceed to describe ways to simulate the model, either by projecting 3D volumes in a set of 2D particles, or by accurately forward modelling the imaging process.
### Microscope Architecture
We start at the electron gun, go through the sample and end at the detector.
#### Electron Gun (shot-noise, coherence, ...)
#### Sample delivery
#### Detector (DQE, MTF, ...)
### Theoretical Limits
Going through the events spanning the emission of electrons to their detection, we list the various aspects of the imaging process that limit the scope of what could be imaged in principle.
#### CTF
#### Origin of Noise
#### Beam Induced Motion
### Simulations
Two main simulation approaches can be distinguished. They are related but differ in their level of "realism" and the use one has of them.
#### Simulate particles 
*see notebook on [(github)](notebooks/ImageFormationModel/Projection.ipynb) [(colab)](https://colab.research.google.com/github/slaclab/cryoEM-notebooks/blob/master/notebooks/ImageFormationModel/Projection.ipynb)*
#### Simulate micrographs 
*see notebook on [(github)](notebooks/Simulating%20data.ipynb) [(colab)](https://colab.research.google.com/github/slaclab/cryoEM-notebooks/blob/master/notebooks/ImageFormationModel/Simulating%20data.ipynb)*
#### Tilt Series for Tomo

## Micrographs: movies to picture
### Beam Induction Motion Correction

## Particles
### Particle Picking
### 2D Classification and Particle Rejection

## 3D Reconstruction - Inverse Problem
### Pose Estimation
### Denoising Images
### Denoising Volumes
### Deconvolution
#### Wiener filtering
### Volume Reconstruction

## Heterogeneity
Here we survey the various approaches available on the market to address the heterogeneity problem in single particle analysis with cryoEM. We also provide convenience tools to help the user prepare inputs and interpret outputs from these tools.
### Multibody Refinement with RELION - prepare input 
*see notebook on [(github)](notebooks/VolumeHeterogeneity/Automated%20Body%20Definition.ipynb) [(colab)](https://colab.research.google.com/github/slaclab/cryoEM-notebooks/blob/master/notebooks/VolumeHeterogeneity/Automated%20Body%20Definition.ipynb)*
In this notebook, we show how to prepare the body masks necessary to run RELION's multibody refinement tool, given a local resolution map as input.

### 3D Variability with cryoSPARC
### Visualization
*see notebook on [(github)](notebooks/VolumeHeterogeneity/Latent%20Space%20Visualizer.ipynb) [(colab)](https://colab.research.google.com/github/slaclab/cryoEM-notebooks/blob/master/notebooks/VolumeHeterogeneity//Latent%20Space%20Visualizer.ipynb)*
In this notebook, we show how to visualize the outputs from RELION's multibody refinement tool and cryoSPARC 3D Variability tool.
### Heterogeneity Analysis
*see notebook on [(github)](notebooks/VolumeHeterogeneity/Latent%20Space%20Analyzer.ipynb) [(colab)](https://colab.research.google.com/github/slaclab/cryoEM-notebooks/blob/master/notebooks/VolumeHeterogeneity//Latent%20Space%20Analyzer.ipynb)*
In this notebook, we show how to analize the outputs from RELION's multibody refinement tool and cryoSPARC 3D Variability tool.

## Tomograms vs SPA (2d vs 3d)
### Tomogram Reconstruction Methods

## Interfaces and Pipelines
### Image formats and IO 
*see notebook on [(github)](notebooks/Interfaces/Simple%20IO%20and%20Visualisation.ipynb) [(colab)](https://colab.research.google.com/github/slaclab/cryoEM-notebooks/blob/master/notebooks/Interfaces/Simple%20IO%20and%20Visualisation.ipynb)*
### Create particle stacks 
*see notebook on [(github)](notebooks/Interfaces/Create%20Particle%20Stacks.ipynb) [(colab)](https://colab.research.google.com/github/slaclab/cryoEM-notebooks/blob/master/notebooks/Interfaces/Create%20Particle%20Stacks.ipynb)*
### Interfacing with OSF 
*see notebook on [(github)](notebooks/Interfaces/Using%20OSF%20to%20retrieve%20and%20store%20datasets.ipynb) [(colab)](https://colab.research.google.com/github/slaclab/cryoEM-notebooks/blob/master/notebooks/Interfaces/Using%20OSF%20to%20retrieve%20and%20store%20datasets.ipynb)*

