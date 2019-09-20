# cryoEM-notebooks

This repo contains information that we deem useful in understanding the end-to-end computational problems of cryoEM and cryoET. Specifically, we develop jupyter notebooks that walks us through the mathematical and computational operations required to reconstruct density maps from micrographs. Work in progress...

- [Image formation model](notebooks/ImageFormationModel/Image%20Formation%20Model.ipynb)
  - Microscope Architecture
  - Electron Gun (shot-noise, coherence, ...)
  - Sample delivery
  - Detector (DQE, MTF, ...)
  - Theoretical Limits
    - CTF
    - Origin of Noise
    - Beam Induced Motion
  - Simulations
    - [Simulate particles](notebooks/ImageFormationModel/Projection.ipynb)
    - [Simulate micrographs](notebooks/Simulating%20data.ipynb)
    - Tilt Series for Tomo

- Micrographs: movies to picture
  - Beam Induction Motion Correction

- Particles
  - Particle Picking
  - 2D Classification and Particle Rejection

- 3D Reconstruction - Inverse Problem
  - Pose Estimation
  - Denoising Images
  - Denoising Volumes
  - Deconvolution
    - Wiener filtering
  - Volume Reconstruction

- Heterogeneity
  - [Multibody Refinement with RELION - prepare input](notebooks/VolumeHeterogeneity/Automated%20Body%20Definition.ipynb)
  - 3D Variability with cryoSPARC
  - Visualization
  - Heterogeneity Analysis

- Tomograms vs SPA (2d vs 3d)
  - Tomogram Reconstruction Methods

- Interfaces and Pipelines
  - [Image formats and IO](notebooks/Interfaces/Simple%20IO%20and%20Visualisation.ipynb)
  - [Create particle stacks](notebooks/Interfaces/Create%20Particle%20Stacks.ipynb)
  - [Interfacing with OSF](notebooks/Interfaces/Using%20OSF%20to%20retrieve%20and%20store%20datasets.ipynb)


