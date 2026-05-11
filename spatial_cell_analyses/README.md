<img width="303" height="102" alt="neuroesc_long" src="https://github.com/user-attachments/assets/e185a933-e27d-4436-ab1f-52a4fe389e38" />

# Spatial cell analysis (Matlab)
A collection of MATLAB functions for simulating and analyzing the firing properties of spatially modulated neurons, including place cells, grid cells, and head direction (HD) cells. This toolkit provides robust methods for generating synthetic data, computing spatial and directional rate maps, extracting information-theoretic metrics, and performing statistical validation via bootstrapping and shuffling.

## Included Functions
### 1. Synthetic Data Generation
* **`get_synth_neurons.m`**
  Generates a realistic correlated random walk trajectory for a virtual animal in a square arena. It also simulates synthetic spike trains for place cells, grid cells (using an interference model), and head direction Cells (using Von Mises tuning) via an inhomogeneous Poisson process. This is useful for unit-testing downstream analysis pipelines. It is a basic spatial simulation and does not model things like speed modulation, theta modulation or phase-precession.
<img width="1000" alt="image" src="https://github.com/user-attachments/assets/e77f169b-b887-4ba5-b055-3fc0af0ce24d" />

<details>
  <summary>Click to view the code example</summary>

  ```matlab
figure('Color','w')
tiledlayout(3,5,"TileSpacing",'tight')
N = 5;
[trajectory,cells] = get_synth_neurons(N,0,0,"hpc_widths",[2 40]);
for ii = 1:N
    [rmap,dmap] = rate_mapper([trajectory.x(:) trajectory.y(:)].*10,[cells.PC(ii).x(:) cells.PC(ii).y(:)].*10);
    m = get_spatial_info(dmap,rmap);
    nexttile
    imagesc(rmap)
    daspect([1 1 1])
    title(sprintf('Place Cell %d\nSpatial info: %.2fb/s',ii,m.skaggs_si_bits_per_sec),'FontWeight','normal'); 
end

N = 5;
[trajectory,cells] = get_synth_neurons(0,N,0);
for ii = 1:N
    [rmap,dmap] = rate_mapper([trajectory.x(:) trajectory.y(:)].*10,[cells.Grid(ii).x(:) cells.Grid(ii).y(:)].*10);
    m = get_spatial_info(dmap,rmap);
    nexttile
    imagesc(rmap)
    daspect([1 1 1])
    title(sprintf('Grid Cell %d\nSpatial info: %.2fb/s',ii,m.skaggs_si_bits_per_sec),'FontWeight','normal'); 
end

N = 5;
[trajectory,cells] = get_synth_neurons(0,0,N);
for ii = 1:N
    [rmap,dmap] = rate_mapper([trajectory.x(:) trajectory.y(:)].*10,[cells.HD(ii).x(:) cells.HD(ii).y(:)].*10);
    m = get_spatial_info(dmap,rmap);
    nexttile
    imagesc(rmap)
    daspect([1 1 1])
    title(sprintf('HD Cell %d\nSpatial info: %.2fb/s',ii,m.skaggs_si_bits_per_sec),'FontWeight','normal'); 
end
```
</details>

### 2. Rate Map Generation
* **`rate_mapper_hd.m`**
  Computes 1D directional firing rate and dwell time maps. It supports both standard boxcar-smoothed histogram methods and Kernel Spatial Density Estimation (KSDE) to map spike and position data into the head direction domain. See also [rate_mapper](https://github.com/Neuroesc/rate_mapper) for the 2D spatial equivalent.

### 3. Spatial & Directional Metrics
* **`get_spatial_info.m`**
  Calculates classic spatial information theoretic measures from 1D, 2D, or 3D rate maps. Output metrics include Skaggs spatial information content (bits/sec and bits/spike), sparsity, Shannon's entropy, Mutual Information, and Kullback–Leibler divergence.
* **`get_directional_info.m`**
  Calculates standard head direction statistics from directional rate maps. Outputs include the Rayleigh mean resultant vector length, the preferred firing direction (PFD), circular mean, and circular standard deviation.
* **`get_grid_score.m`**
  Calculates grid score from 2D rate maps. Outputs include the grid score (based on a number of published methods), grid spacing and orientation. Also includes quantification of 'average field' values based on the central peak of the autocorrelogram, such as average field radius.
* **`ndautoCORR.m`**
  Calculates autocorrelogram, which is used to calculate grid score. This function has been adapted to work on three-dimensional firing rate maps.
  
### 4. Statistical Validation
* **`get_spatial_bootshuff.m`**
  Computes spatial, grid, and directional metrics for a single cell by generating a bootstrapped "observed" value (via spike resampling with replacement) to control for unstable activity. It then computes exact non-parametric p-values (and parametric z-scores) by comparing the bootstrapped value against a null distribution generated via circular time-shifting of the spike train. This is based on the Method described by [Savelli et al. (2017)](https://elifesciences.org/articles/21354), with some improvements.

<img width="1000" alt="image" src="https://github.com/user-attachments/assets/e4f7167f-023a-4532-a96c-260ff6aee113" />

<details>
  <summary>Click to view the code example</summary>

  ```matlab
% get synthetic data
[trajectory,cells] = get_synth_neurons(1,1,1,"hpc_widths",[2 8]);

% settings for map generation
rmset.ppm = 100;
rmset.method = 'histogram';
rmset.binsize = 25;
rmset.ssigma = 30;
rmset.hd_type = 'histogram';
rmset.hd_bins = 60;
rmset.hd_boxcar = 3;

% bootstrap vs shuffle
res = get_spatial_bootshuff(trajectory.x(:),trajectory.y(:),trajectory.t(:),rmset,'spt',cells.PC(1).t,'poh',trajectory.hd(:),'scores',{'all'});
disp(res)
```
</details>

## Requirements
* **[Circular Statistics Toolbox](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)** (Required for `get_directional_info.m` to compute circular standard deviation, means, and Rayleigh vectors).
* **[rate_mapper](https://github.com/Neuroesc/rate_mapper)** (Required for `get_spatial_info.m` to generate firing rate maps).

## Basic Usage
Most functions are designed with an `arguments` block and support Name-Value pairs for flexible processing. Detailed usage instructions, expected inputs, and outputs can be found in the documentation block at the top of each file. 

