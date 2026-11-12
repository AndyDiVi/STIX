# Laboratory 3 – Solar Flare X-ray Image Reconstruction using CLEAN

<img width="1897" height="523" alt="cleanloop" src="https://github.com/user-attachments/assets/9e9da842-fc9c-481f-969b-96f4fdbe7f8a" />


Authors: Di Via & Hiller-Szala

## Overview

This laboratory focuses on reconstructing solar flare X-ray images from visibility measurements acquired by the Spectrometer/Telescope for Imaging X-rays (STIX), onboard ESA's Solar Orbiter.

Unlike a conventional imaging system, STIX encodes spatial information through modulated X-ray signals. These measurements provide a limited set of Fourier components, called visibilities, from which the original source intensity distribution must be reconstructed.

The objective is to implement the CLEAN algorithm in MATLAB to recover an estimate of the X-ray emission distribution from sparsely sampled Fourier data.

The laboratory analyzes a solar flare observed on April 9, 2022, in the 6–10 keV energy range.

## Instrument Background

STIX uses pairs of front and rear grids to modulate incoming X-rays, producing moiré patterns on segmented detectors.

These patterns encode the amplitude and phase of specific Fourier components of the solar X-ray source.

The instrument includes 30 imaging grid-detector systems. Each system provides information about a particular spatial frequency.

Since only a limited number of Fourier components are available, direct back projection produces a blurred image with sidelobe artifacts. The CLEAN algorithm addresses this reconstruction problem by iteratively removing contributions from the instrument's point spread function.

## Dataset

The input dataset is stored in `data2.xlsx` and contains:

- 24 measured complex visibilities.
- 4 columns describing spatial frequencies and visibility components.
- Measurements in the 6–10 keV energy range.

### Input structure

| Column | Description |
|---|---|
| 1 | Spatial frequency $u_j$, in $\mathrm{arcsec}^{-1}$ |
| 2 | Spatial frequency $v_j$, in $\mathrm{arcsec}^{-1}$ |
| 3 | Real part of the visibility |
| 4 | Imaginary part of the visibility |

Each complex visibility is assembled as:

$$
V_j = \operatorname{Re}(V_j) + i\operatorname{Im}(V_j)
$$

The data is imported in MATLAB using:

```matlab
data = readmatrix("data2.xlsx");

u = data(:, 1);
v = data(:, 2);
Vis = data(:, 3) + 1i * data(:, 4);
```

## Mathematical Background

### Fourier visibilities

Let $\phi(x,y)$ denote the spatial intensity distribution of the X-ray source.

Using the Fourier sign convention adopted in this laboratory, the visibility function is:

$$
V(u,v)
=
\iint_{\mathbb{R}^2}
\phi(x,y)e^{2\pi i(ux+vy)}
\,dx\,dy
$$

With complete Fourier coverage, the original distribution could be recovered through:

$$
\phi(x,y)
=
\iint_{\mathbb{R}^2}
V(u,v)e^{-2\pi i(ux+vy)}
\,du\,dv
$$

In this laboratory, only a sparse set of visibility measurements is available, making reconstruction a deconvolution problem.

### Conjugate symmetry

For a real-valued intensity distribution, the visibility function satisfies:

$$
V(-u,-v)=\overline{V(u,v)}
$$

This property is used to extend the dataset from 24 measured points to 48 Fourier samples by adding their conjugate counterparts.

```matlab
u_total = [u; -u];
v_total = [v; -v];
Vis_total = [Vis; conj(Vis)];
```

The additional samples enforce conjugate symmetry; they are not independent measurements.

## Methodology

The reconstruction workflow consists of five main stages.

### 1. Visibility analysis

The sampled spatial frequencies are displayed in the $u$–$v$ plane to inspect Fourier coverage.

The implementation also visualizes:

- Real parts of the visibilities.
- Imaginary parts of the visibilities.
- Visibility magnitudes.

These plots describe the information available for image reconstruction.

### 2. Dirty Map reconstruction

The Dirty Map is obtained through a discrete back projection of the visibility measurements:

$$
DM(x_m,y_n)
=
\sum_{j=1}^{N_v}
V_j e^{-2\pi i(u_jx_m+v_jy_n)}
$$

where $N_v=48$ after conjugate extension.

The spatial grid contains $128\times128$ pixels with a sampling interval of 1 arcsecond:

$$
x_m=-64,\ldots,63
$$

$$
y_n=63,\ldots,-64
$$

The implementation rotates the resulting map by $90^\circ$ counterclockwise for the coordinate convention used in the laboratory.

Both real and imaginary parts are displayed. The real part contains the reconstructed signal, while the imaginary part is negligible.

### 3. Dirty Beam construction

The Dirty Beam represents the point spread function associated with the sampled Fourier frequencies.

Before normalization, it is computed as:

$$
DB_{\mathrm{raw}}(x,y)
=
\sum_{j=1}^{N_v}
e^{-2\pi i(u_jx+v_jy)}
$$

It is then normalized:

$$
DB(x,y)
=
\frac{DB_{\mathrm{raw}}(x,y)}
{\max_{x,y}\left|DB_{\mathrm{raw}}(x,y)\right|}
$$

The Dirty Beam is evaluated on a larger $256\times256$ grid:

$$
x=-128,\ldots,127
$$

$$
y=127,\ldots,-128
$$

This larger grid allows shifted beam sections to be extracted during CLEAN without wrapping the beam around the image boundaries.

The beam is also rotated by $90^\circ$ counterclockwise.

### 4. CLEAN deconvolution

CLEAN iteratively represents the source as a collection of point components.

At each iteration, the algorithm:

1. Finds the brightest point in the residual Dirty Map.
2. Aligns the Dirty Beam with that point.
3. Subtracts a fraction of the corresponding beam contribution.
4. Adds the extracted amplitude to the Cleaned Component Map.

#### Initialization

The residual map is initialized with the Dirty Map, while the component map is initialized to zero:

$$
R^{(0)}=DM
$$

$$
CC^{(0)}=0
$$

#### Peak identification

At iteration $k$, the algorithm identifies the maximum of the real residual:

$$
p_k
=
\max_{x,y}\operatorname{Re}\left(R^{(k)}(x,y)\right)
$$

$$
(x_k,y_k)
=
\operatorname*{arg\,max}_{x,y}
\operatorname{Re}\left(R^{(k)}(x,y)\right)
$$

#### Residual update

Since the Dirty Beam has a unit peak, the residual update is:

$$
R^{(k+1)}(x,y)
=
R^{(k)}(x,y)
-
\gamma p_k DB(x-x_k,y-y_k)
$$

where $\gamma$ is the gain factor.

#### Component update

The extracted component is accumulated at the selected pixel:

$$
CC^{(k+1)}(x,y)
=
CC^{(k)}(x,y)
+
\gamma p_k\delta_{x,x_k}\delta_{y,y_k}
$$

Here, $\delta$ denotes the Kronecker delta on the discrete image grid.

#### Stopping criterion

The algorithm stops when the maximum real residual is at or below the selected threshold:

$$
p_k\leq\tau
$$

Otherwise, it continues until the maximum number of iterations is reached.

### CLEAN parameters

| Parameter | Value |
|---|---:|
| Maximum iterations | 300 |
| Gain factor $\gamma$ | 0.01 |
| Residual threshold $\tau$ | 2000 |
| Dirty Map dimensions | $128\times128$ |
| Dirty Beam dimensions | $256\times256$ |

The algorithm is called using:

```matlab
[DM_final, CC_final, DB_shifted] = ...
    clean_algorithm(DM, DB, 300, 0.01, 2000);
```

The reported execution reaches the stopping threshold at iteration 258.

### 5. Clean Beam and final image

The Cleaned Component Map is restored using an isotropic Gaussian Clean Beam centered at the origin:

$$
CB_{\mathrm{raw}}(x,y)
=
\frac{1}{2\pi\sigma^2}
\exp\left(
-\frac{x^2+y^2}{2\sigma^2}
\right)
$$

The standard deviation is derived from the selected Gaussian full width at half maximum:

$$
\sigma
=
\frac{\mathrm{FWHM}}
{2\sqrt{2\ln 2}}
$$

The implementation uses a restoring-beam width parameter of 28 arcseconds:

$$
\mathrm{FWHM}=28\ \mathrm{arcsec}
$$

The sampled Clean Beam is normalized to have a unit peak:

$$
CB(x,y)
=
\frac{CB_{\mathrm{raw}}(x,y)}
{\max_{x,y}CB_{\mathrm{raw}}(x,y)}
$$

Finally, the Clean Map is obtained by convolving the Cleaned Component Map with the Clean Beam:

$$
CM=CC*CB
$$

In MATLAB:

```matlab
FWHM = 28;
sigma = FWHM / (2 * sqrt(2 * log(2)));

CB = (1 / (2 * pi * sigma^2)) * ...
    exp(-(X1.^2 + Y1.^2) / (2 * sigma^2));

CB = CB / max(CB(:));

CM = conv2(CC_final, CB, 'same');
```

This restoration step replaces the discrete point components with smooth Gaussian contributions.

The implementation constructs the final image from the restored component map without adding the final residual map.

## Results

The reported CLEAN execution produces:

| Result | Value or description |
|---|---|
| Reported stopping iteration | 258 |
| Gain factor | 0.01 |
| Residual stopping threshold | 2000 |
| Cleaned Component Map | Loop-like structure |
| Final image dimensions | $128\times128$ |
| Restoring-beam width parameter | 28 arcseconds |

The reconstruction identifies a loop-like emission structure in the Cleaned Component Map and produces a smooth final image after Gaussian restoration.

The stopping message indicates that the selected residual threshold has been reached; it does not imply an exact recovery of the original source.

## Implementation

The laboratory is implemented in MATLAB.

The main components are:

- Data loading from `data2.xlsx`.
- Complex visibility construction.
- Conjugate-symmetric extension of the Fourier samples.
- Visibility and sampling visualization.
- Dirty Map back projection.
- Dirty Beam computation and normalization.
- Iterative CLEAN deconvolution.
- Gaussian Clean Beam construction.
- Final image restoration through convolution.

### Custom CLEAN function

The custom function has the following interface:

```matlab
function [DM, CC, DB_shifted] = ...
    clean_algorithm(DM, DB, max_iterations, gamma, th)
```

#### Inputs

| Input | Description |
|---|---|
| `DM` | Initial Dirty Map |
| `DB` | Normalized Dirty Beam |
| `max_iterations` | Maximum number of iterations |
| `gamma` | Gain factor |
| `th` | Residual stopping threshold |

#### Outputs

| Output | Description |
|---|---|
| `DM` | Final residual Dirty Map |
| `CC` | Final Cleaned Component Map |
| `DB_shifted` | Shifted Dirty Beam used in the last subtraction |

The function validates the input matrices and scalar parameters, initializes the component map, and performs the iterative beam subtraction.

## Visualizations

The implementation generates the following plots:

- Sampling distribution in the $u$–$v$ plane.
- Real parts of the visibilities.
- Imaginary parts of the visibilities.
- Visibility magnitudes.
- Real and imaginary parts of the Dirty Map.
- Sampled transfer function.
- Real and imaginary parts of the Dirty Beam.
- Three-dimensional Dirty Beam surface.
- Shifted Dirty Beam from the final subtraction.
- Final residual Dirty Map.
- Cleaned Component Map.
- Central Dirty Beam sections along $x=0$ and $y=0$.
- Gaussian Clean Beam.
- Final Clean Map.

## Requirements

- MATLAB.
- Image Processing Toolbox for the `imrotate` calls used in the script.
- The input dataset `data2.xlsx`.

If the CLEAN function is stored separately, `clean_algorithm.m` must be available in the working directory or on the MATLAB path.

## How to run

1. Clone or download this repository.
2. Place `data2.xlsx` in the MATLAB working directory.
3. Ensure that the main script and the `clean_algorithm` function are available.
4. Open the main MATLAB script.
5. Execute the script.

The script loads the visibility data, displays the Fourier sampling, computes the Dirty Map and Dirty Beam, runs CLEAN, and generates the final restored X-ray image.

## Conclusions

This laboratory demonstrates an image reconstruction workflow for sparsely sampled Fourier measurements from STIX.

The CLEAN algorithm progressively removes shifted Dirty Beam contributions from the back-projected image and accumulates a model of the source in the Cleaned Component Map.

With a gain factor of $\gamma=0.01$ and a residual threshold of 2000, the reported execution reaches its stopping criterion at iteration 258.

Convolving the extracted components with a Gaussian Clean Beam produces a smooth representation of the reconstructed solar flare emission.

The final reconstruction remains dependent on the available Fourier coverage, the stopping threshold, the gain factor, and the selected restoring beam.

---

Laboratory 3 – Di Via & Hiller-Szala
