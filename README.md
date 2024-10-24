# Classification of mental tasks Using EEG Signals with Common spatial Pattern (CSP) and Particle Swarm Optimization (PSO) 

This repository contains MATLAB code for classifying Mental tasks using EEG signals, enhanced through the application of PSO. The associated research paper was published in the Elixir International Journal under the title "Performance Analysis of PSO and GA Algorithms in Order to Classifying EEG Data". The only difference is that CSP is used instead of PCA in this code.

## Requirements

This code was written and run in MATLAB 2009 software

## Main Parameters

The bellow values were considered for main parameters. You can set other values:

### Num segments: win= 20
### Size of each segment= 125
### Num epochs: epoch= 30
### Num particles: num_chrom= 110    

## Warning

if there is an error about NaN data, search them in test and train matrix and set them to zero.

## Dataset used

The proposed algorithm was implemented on a dataset comprising five mental tasks recorded at the Colorado Electroencephalography and Brain-Computer Interfaces Laboratory (CEBL) by Keirn and Aunon. These tasks include relaxation, 3D shape rotation around an axis, mental multiplication, letter writing to a friend, and visual counting (refer to Fig. 6). These signals were captured at a sampling rate of 250 samples per second across six channels: C1, C2, O1, O2, P1, and P2, each of which lasts for 10 seconds and contains 2500 samples. The dataset encompasses a total of 325 signals.

## Acknowledgments

We would like to acknowledge the contributions of the researchers and participants involved in the development of the Keirn's dataset, as well as the community that supports open-source projects.
