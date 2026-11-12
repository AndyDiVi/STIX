## Mathematical Background
<img width="1897" height="523" alt="cleanloop" src="https://github.com/user-attachments/assets/242f0099-8e37-4db3-91a4-453fb8aa2fa4" />

### Fourier visibilities

Let $\phi(x,y)$ denote the spatial intensity distribution of the X-ray source.

Using the Fourier sign convention adopted in this laboratory, the visibility function is:

```math
V(u,v) = \iint_{\mathbb{R}^2} \phi(x,y)e^{2\pi i(ux+vy)}\,dx\,dy
```

With complete Fourier coverage, the original distribution could be recovered through:

```math
\phi(x,y) = \iint_{\mathbb{R}^2} V(u,v)e^{-2\pi i(ux+vy)}\,du\,dv
```

Each complex visibility is assembled as:

```math
V_j = \mathrm{Re}(V_j) + i\,\mathrm{Im}(V_j)
```

### Conjugate symmetry

For a real-valued intensity distribution, the visibility function satisfies:

```math
V(-u,-v) = \overline{V(u,v)}
```

This property extends the dataset from 24 measured points to 48 Fourier samples by adding their conjugate counterparts.

The additional samples enforce conjugate symmetry; they are not independent measurements.

## Methodology

### 1. Dirty Map reconstruction

The Dirty Map is obtained through a discrete back projection of the visibility measurements:

```math
DM(x_m,y_n) = \sum_{j=1}^{N_v} V_j e^{-2\pi i(u_jx_m+v_jy_n)}
```

Here, $N_v = 48$ after conjugate extension.

The spatial grid contains $128\times128$ pixels with a sampling interval of 1 arcsecond:

```math
x_m = -64,\ldots,63
```

```math
y_n = 63,\ldots,-64
```

The implementation rotates the resulting map by $90^\circ$ counterclockwise for the coordinate convention used in the laboratory.

### 2. Dirty Beam construction

Before normalization, the Dirty Beam is computed as:

```math
DB_{\mathrm{raw}}(x,y) = \sum_{j=1}^{N_v} e^{-2\pi i(u_jx+v_jy)}
```

It is then normalized:

```math
DB(x,y) = \frac{DB_{\mathrm{raw}}(x,y)}{\max_{x,y}\left|DB_{\mathrm{raw}}(x,y)\right|}
```

The Dirty Beam is evaluated on a larger $256\times256$ grid:

```math
x = -128,\ldots,127
```

```math
y = 127,\ldots,-128
```

This larger grid allows shifted beam sections to be extracted during CLEAN without wrapping the beam around the image boundaries.

### 3. CLEAN deconvolution

#### Initialization

The residual map is initialized with the Dirty Map, while the component map is initialized to zero:

```math
R^{(0)} = DM
```

```math
CC^{(0)} = 0
```

#### Peak identification

At iteration $k$, the algorithm identifies the maximum of the real residual:

```math
p_k = \max_{x,y}\mathrm{Re}\left(R^{(k)}(x,y)\right)
```

```math
(x_k,y_k) = \arg\max_{(x,y)}\mathrm{Re}\left(R^{(k)}(x,y)\right)
```

#### Residual update

Since the Dirty Beam has a unit peak, the residual update is:

```math
R^{(k+1)}(x,y) = R^{(k)}(x,y) - \gamma p_k DB(x-x_k,y-y_k)
```

Here, $\gamma$ is the gain factor.

#### Component update

The extracted component is accumulated at the selected pixel:

```math
CC^{(k+1)}(x,y) = CC^{(k)}(x,y) + \gamma p_k\delta_{x,x_k}\delta_{y,y_k}
```

Here, $\delta$ denotes the Kronecker delta on the discrete image grid.

#### Stopping criterion

The algorithm stops when the maximum real residual is at or below the selected threshold:

```math
p_k \leq \tau
```

Otherwise, it continues until the maximum number of iterations is reached.

### 4. Clean Beam and final image

The Cleaned Component Map is restored using an isotropic Gaussian Clean Beam centered at the origin:

```math
CB_{\mathrm{raw}}(x,y) = \frac{1}{2\pi\sigma^2}\exp\left(-\frac{x^2+y^2}{2\sigma^2}\right)
```

The standard deviation is derived from the selected Gaussian full width at half maximum:

```math
\sigma = \frac{\mathrm{FWHM}}{2\sqrt{2\ln 2}}
```

The implementation uses a restoring-beam width parameter of 28 arcseconds:

```math
\mathrm{FWHM} = 28\ \mathrm{arcsec}
```

The sampled Clean Beam is normalized to have a unit peak:

```math
CB(x,y) = \frac{CB_{\mathrm{raw}}(x,y)}{\max_{x,y}CB_{\mathrm{raw}}(x,y)}
```

Finally, the Clean Map is obtained by convolving the Cleaned Component Map with the Clean Beam:

```math
CM = CC * CB
```

Here, $*$ denotes two-dimensional convolution.
