Usage
  ./bin/multiscatter [options] data_in.dat > data_out.dat
  ./bin/multiscatter [options] < data_in.dat > data_out.dat

General options
  -help           Display this message
  -repeat n       Repeat algorithm n times (for benchmarking)
  -quiet          Don't report activity to stderr
  -v1             Use version 1.x interpretation of first line of input file
  -auto           Automatically select algorithm settings from wavelength etc
  -radar          Use settings appropriate for radar
  -lidar          Use settings appropriate for lidar
  -algorithms <small_angle_algorithm> <wide_angle_algorithm>
                  Manually select algorithm, where <small_angle_algorithm> is:
      none     - No single or small-angle scattering
      single   - Single scattering (no small-angle scattering)
      original - Original Hogan (2006) algorithm: speed O(N^2)
      fast     - Faster Hogan (2008) algorithm: speed O(N) DEFAULT
      explicit - Eloranta-like explicit representation of each order of scattering
      lag      - Hogan (2008) but also with lag calculation (m)
                  ...and where <wide_angle_algorithm> is:
      none     - No wide-angle scattering
      tdts     - Time-dependent two stream (Hogan and Battaglia 2008)
      lidar    - TDTS assuming forward lobe DEFAULT
      radar    - TDTS assuming no forward lobe
  -hsrl           Output particulate and air backscatter separately
  -gaussian-receiver
                  Receiver is Gaussian rather than top-hat shaped
  -jacobian       Output the approximate but fast Jacobian
  -numerical-jacobian
                  Output the Jacobian calculated (slowly) using finite
                  differences
  -ext-only       Only calculate the Jacobian with respect to extinction
  -adjoint        Output the adjoint as well
  -check-adjoint  Calculate the adjoint and check it with a numerical Jacobian
  -annular-detectors
                  2nd and subsequent detectors are ring shaped

Options for ORIGINAL small-angle (SA) algorithm
  -simple-optical-depth
                  Use simple optical depth integration
  -crude-optical-depth
                  Use crude optical depth integration
  -crude-integration
                  Use crude double/multiple scattering integration
  -no-multiscat-within-gate
                  Photon cannot be forward scattered more than once within a
                  single gate; note that this has the opposite effect to
                  "-crude-integration"
  -double-scattering-only
                  Don't include triple and higher-order scatterings
  -no-molecular-extinction
                  Molecules do not extinguish the beam
  -wide-angle-cutoff <theta>
                  Forward scattering at angles greater than <theta> radians
                  are deemed to escape, a crude way to deal with a problem
                  associated with aerosols
  -crude-double-scattering
                  Don't use Eloranta's slow but accurate double scattering
                 formulation
  -approx-exp     Appriximate the exp() function call for speed

Options for EXPLICIT small-angle (SA) algorithm
  -explicit n     Use an explicit model with n orders of scattering
  -output-distribution n dx
                  Output horizontal photon distributions at n points spaced
                  dx apart, starting at dx/2

Options for wide-angle multiple-scattering
  -ignore-source-widening
                  Ignore widening effect of small-angle scattering on source beam for wide-angle calc
  -optimize-wide-angle-gates
                  Skip optically thin gates in wide-angle calculation
  -simple-2s-coeffts
                  Use the simple upwind Euler formulation (can be unstable
                  for high optical depth)
  -ssa-scales-forward-lobe
                  Single-scattering albedo less than unity reduces the
                  forward scattering lobe as well as wide-angle scattering
  -num-samples m  Output m samples, allowing sampling of signals appearing
                  to originate below ground
  -propagation-to-stderr
                  Print the diffuse photon energies and variances to
                  standard error

Input data
  First line: 5 values
    1: Number of range gates in profile
    2: Wavelength (m)
    3: Altitude of instrument (m)
    4: Transmitter divergence, 1/e half-width (radians)
    5: Receiver field-of-view, half-width of first receiver (radians)
   6+: Fields-of-view for any other channels (radians)
  (note that with the "-v1" option, the version 1.x file format is assumed,
   in which the order is (1) number of range gates, (2) wavelength, (3) divergence,
   (4) field-of-view, (5) altitude, and only one field of view is allowed)
  Subsequent lines: 4, 5, 8, 10 or more values (all 8 required for the wide-angle
  calculation; all 10 to represent anisotropic phase functions near 180 deg);
  more to calculate the adjoint
    1: Range of gate above ground starting with nearest gate to instrument (m)
    2: Extinction coefficient of cloud/aerosol only (m-1)
    3: Equivalent-area cloud/aerosol radius (m)
    4: Extinction-to-backscatter ratio of cloud/aerosol (sr)
    5: Air extinction coefficient (m-1); optional (zero if missing)
    6: Single scattering albedo of cloud/aerosol
    7: Scattering asymmetry factor of cloud/aerosol
    8: Single scattering albedo of air
    9: Fraction of cloud/aerosol backscatter due to droplets (real number
         between 0 and 1), for representing anisotropic near-180 phase
         functions in the SA algorithm
   10: Pristine ice fraction, for anisotropic near-180 phase functions
         from Yang et al. (2000); note that irregular ice particles
         are believed to have isotropic near-180 phase functions, in
         which case set this variable to zero
  Subsequent lines (only with "-adjoint" option)
 11 to 10+n: Adjoint input for particulate backscatter of the n fields-of-view
 11+n: Adjoint input for air backscatter of the first field-of-view (for HSRL)

Output data (without -jacobian or -numerical-jacobian option)
  Default: 5 values per range gate
    1: Index of range gate
    2: Apparent range above ground (m)
    3: Extinction coefficient (m-1)
    4: Equivalent-area particle radius (microns)
    5: Apparent backscatter (m-1 sr-1), for particulates only if with "-hsrl"
  With the "-hsrl" option
    6: Apparent backscatter of the air only (m-1 sr-1)
  With the "-adjoint" option (where n=6 if "-hsrl", n=5 otherwise)
  n+1: Adjoint of extinction coefficient (m)
  n+2: Adjoint of single-scattering albedo
  n+3: Adjoint of asymmetry factor
  n+4: Adjoint of extinction-to-backscatter ratio (m sr)

Output data (with -jacobian or -numerical-jacobian option)
   n x m matrix: d(attenuated backscatter) / d(extinction coefficient)
and if "-ext-only" is not set:
   n x m matrix: d(attenuated backscatter) / d(single-scattering albedo)
   n x m matrix: d(attenuated backscatter) / d(asymmetry factor)
   n x m matrix: d(attenuated backscatter) / d(particle radius)
   1 x m vector: d(attenuated backscatter) / d(ext-to-bscat ratio)