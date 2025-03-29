#
# conda create --prefix=~/sw/circ_test --override-channels --channel defaults python==3.11
# conda activate ~/sw/circ_test
# pip install matplotlib numpy scipy ngsolve h5py psutil
#

import sys, os

from ngsolve import *

from netgen.geom2d import SplineGeometry

import numpy as np

import scipy.special as scs
from scipy.sparse.linalg import lsqr, cg
from scipy.linalg import svd
from scipy.linalg import diagsvd

import h5py
import logging
import gc
import psutil
import time
from concurrent.futures import ProcessPoolExecutor
import argparse
import fem
import custom_mesh
import matrix_form

#----------------------HELPER FUNCTIONS---------------------------------------------------------
# Format byte count to human-readable form:
def bytes_hr(byte_count, ignore_non_prefixed=True):
    if byte_count <= 1024 and ignore_non_prefixed:
        return ''
    prefix_strs = ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi']
    prefix_idx = 0
    byte_count = float(byte_count)
    while byte_count > 1024.0 and prefix_idx < len(prefix_strs):
        prefix_idx += 1
        byte_count /= 1024.0
    return '{:.1f} {:s}B'.format(byte_count, prefix_strs[prefix_idx])

# Log sense of memory usage:
def log_memory_usage():
    rss = psutil.Process().memory_info().rss
    rss_hr = bytes_hr(rss)
    if rss_hr:
        logging.debug('Current RSS:  %d (%s)', rss, rss_hr)
    else:
        logging.debug('Current RSS:  %d', rss)

# Coefficient function for mesh
def create_ncoef(mesh, nval): # CHANGE FOR GEO
    ncoef=CoefficientFunction(
            [nval[0] if mat=="circle1" else 
             nval[1] if mat=="circle2" else 
             nval[2] if mat=="circle3" else 
             nval[3] if mat=="circle4" else 
             1 for mat in mesh.GetMaterials()])
    scatter=CoefficientFunction(
            [6 if mat=="circle1" else 
             5 if mat=="circle2" else 
             4 if mat=="circle3" else 
             3 if mat=="circle4" 
             else 1 if mat=="air" 
             else 2 for mat in mesh.GetMaterials()])
    return ncoef, scatter

# Generate random data for circle
def generate_circles(numcircle): # CHANGE FOR GEO
    circles = []
    while len(circles) < numcircle:
        # Randomly generate radius
        radius = np.random.uniform(0.1, 0.3)  # Ensuring reasonable sizes
        # Randomly generate center coordinates within [-1, 1] x [-1, 1]
        xcen = (np.random.uniform(-1 + radius, 1 - radius), np.random.uniform(-1 + radius, 1 - radius))
        # Check for overlaps with existing circles
        overlaps = False
        for c in circles:
            dx = abs(xcen[0] - c['xcen'][0])
            dy = abs(xcen[1] - c['xcen'][1])
            distance = np.sqrt(dx**2 + dy**2)
            if distance < (radius + c['radius']):
                overlaps = True
                break

        if not overlaps:
            circles.append({
                'radius': radius,
                'xcen': xcen,
                'ind': len(circles) + 3
            })
    # Prepare data in the desired format
    circle_data = {'numcircle': numcircle}
    for i, c in enumerate(circles):
        circle_data.update({
            f'R{i+1}': c['radius'],
            f'xcen{i+1}': c['xcen'],
            f'ind{i+1}': c['ind']
        })
    return circle_data

# Main Loop
def do_iteration(iter_args):
    # Unpack the argument tuple:
    (i, params) = iter_args

    logging.info("Circle: %d", i)
    start_time = int(time.time())
    np.random.seed(start_time + i)
    
    numcircle = np.random.randint(1,4) # Number of circles
    circle_data = generate_circles(numcircle) # Get random scatterer
    nval = np.random.uniform(1.1, 1.8, 3) # Refractive index
    
    logging.info('Nval: %s', nval)
    logging.info('Data: %s', circle_data)

    # See matrix_form.py for available geometries
    # image and mesh functions must me modified for each geometry
    # create_ncoef must be internally modified for each geometry
    image = matrix_form.circle_matrix(params.xlim, params.Ngrid, numcircle, circle_data, nval)
    image = image.flatten() # For storage
    
    hmax_s = params.hmax_a / np.sqrt(np.max(nval)) # Mesh size in scatterer
    hmax = np.max([hmax_s, params.hmax_a]) # Global mesh size
    
    mesh = custom_mesh.circle(hmax_a=params.hmax_a,
                  hmax_s=hmax_s,
                  order=params.porder,
                  pml_rad=params.pml_rad,
                  pml_delta=params.pml_delta,
                  circle_data=circle_data
                  )

    mesh.SetPML(pml.Radial(rad=params.pml_rad, alpha=params.pml_parameter, origin=(0,0)), "pmlregion")
    ncoef, scatter = create_ncoef(mesh, nval)

    logging.info('    Mesh and Matrix Done')

    #------------SOLVE FORWARD PROBLEM FOR FAR FIELD PATTERN----------------------------
    uinf, phi, theta = fem.helmsol(mesh, params.porder, ncoef, params.kappa, params.inc_p, params.far_p)
    uinf = uinf.flatten() # Flatten for storage and inversion

    logging.info('    Solved Forward Problem')

    #------------SOLVE INVERSE PROBLEM FOR CONTRAST-------------------------------------
    born_operator_list = [fem.discretize_born(
                            params.kappa, 
                            params.xlim, 
                            inc_field, 
                            params.Ngrid, 
                            theta
                            ) 
                          for inc_field in phi
                        ]
    born = np.vstack(born_operator_list)
    m_approx = lsqr(born, uinf, damp=1e0)[0] # Dampening can be changed

    logging.info('    Solved Inverse Problem')

    log_memory_usage()
    gc.collect()
    log_memory_usage()

    #------------RETURN THE RESULTS----------------------------------
    return (i, image, m_approx, uinf)


#----------------------------INITIALIZE PARAMETERS--------------------------------------------------------

class Parameters:
    """Instances of the Parameters class collate the parameters that control the model and methods."""

    def __init__(self, **kwargs):
        """Initialize the receiver with default parameter values overridden by values present in the keyword argument list to this method.  E.g.

    p = Parameters(N_sample=200, Idx_sample=0, ninc=100, nfar=100)

"""
        self.N_sample_max = kwargs.get('N_sample_max', 20000)   # Maximum number of samples (in HDF5 file)
        self.N_sample = kwargs.get('N_sample', 3000)             # Number of samples to create/read
        self.Idx_sample = kwargs.get('Idx_sample', 0)           # Base index of samples (in HDF5 file)

        # Domain parameters:
        self.xlim = kwargs.get('xlim', 1)           # Size of grid
        self.Ngrid = kwargs.get('Ngrid', 100)       # Refinement of approximation grid

        # Wave parameters:
        self.ninc = kwargs.get('ninc', 100)          # Number of incident fields, should be >2kR
        self.nfar = kwargs.get('nfar', 100)          # Number of far fields
        self.inc_p = kwargs.get('inc_p',
                    { "n":self.ninc,                # Parameters for incident wave
                      "app":2*np.pi,
                      "cent":0
                    })
        self.far_p = kwargs.get('far_p',            # Parameters for measurement
                    { "n":self.nfar,
                      "app":self.inc_p["app"],
                      "cent":0
                    })
        self.kappa = kwargs.get('kappa', 16)        # Wavenumber
        self.porder = kwargs.get('porder', 4)       # Order of the polynmials in the FEM

        # Computed constants across all iterations:
        self.hmax_a = kwargs.get('hmax_a',          # Mesh size in air
                                2 * np.pi / self.kappa / 8)
        self.pml_rad = kwargs.get('pml_rad',        # Radius of inner edge of the PML
                                1 + 2 * 2 * np.pi / self.kappa)
        self.pml_delta = kwargs.get('pml_delta',    # Thickness of the PML
                                2 * np.pi / self.kappa)
        self.pml_width = kwargs.get('pml_delta',    # Thickness of the PML (again)
                                2 * np.pi / self.kappa)
        self.pml_parameter = kwargs.get('pml_parameter', 1j)

        # Image generation:
        self.should_plot = bool(kwargs.get('should_plot', False))
        if self.should_plot:
            try:
                import matplotlib
            except:
                self.should_plot = False
        self.plot_dir = kwargs.get('plot_dir', './plots')
        self._color_map = None

    def summarize(self):
        summary_str = f"""
----------------------------------
P A R A M E T E R    S U M M A R Y
----------------------------------

    Sampling:

        N_sample_max = {self.N_sample_max}, N_sample = {self.N_sample}, Idx_sample = {self.Idx_sample}

    Domain:

        xlim = {self.xlim}, Ngrid = {self.Ngrid}

    Wave profiles:

        kappa = {self.kappa}, porder = {self.porder}

        Incident:
            ninc = {self.ninc},
            n = {self.inc_p['n']}, app = {self.inc_p['app']}, cent = {self.inc_p['cent']}

        Far field:
            nfar = {self.nfar},
            n = {self.far_p['n']}, app = {self.far_p['app']}, cent = {self.far_p['cent']}

    Computed constants:

        hmax_a = {self.hmax_a},
        pml_rad = {self.pml_rad},
        pml_delta = {self.pml_delta},
        pml_parameter = {np.real(self.pml_parameter)} + {np.imag(self.pml_parameter)} * i

"""
        for log_line in summary_str.splitlines():
            logging.info('%s', log_line)

    def read_from_hdf5_group(self, hdf5_grp):
        """Import receiver's parameter values from an HDF5 dataset group.  Any missing attributes in the file will yield an exception."""
        self.N_sample_max = hdf5_grp.attrs['N_sample_max']

        self.xlim = hdf5_grp.attrs['xlim']
        self.Ngrid = hdf5_grp.attrs['Ngrid']

        self.ninc = hdf5_grp.attrs['ninc']
        self.nfar = hdf5_grp.attrs['nfar']

        self.inc_p['n'] = hdf5_grp.attrs['inc_param.n']
        self.inc_p['app'] = hdf5_grp.attrs['inc_param.app']
        self.inc_p['cent'] = hdf5_grp.attrs['inc_param.cent']

        self.far_p['n'] = hdf5_grp.attrs['far_param.n']
        self.far_p['app'] = hdf5_grp.attrs['far_param.app']
        self.far_p['cent'] = hdf5_grp.attrs['far_param.cent']

        self.kappa = hdf5_grp.attrs['kappa']
        self.porder = hdf5_grp.attrs['porder']

        self.hmax_a = hdf5_grp.attrs['hmax_a']
        self.pml_rad = hdf5_grp.attrs['pml_rad']
        self.pml_delta = hdf5_grp.attrs['pml_delta']
        self.pml_parameter = hdf5_grp.attrs['pml_parameter.real'] + hdf5_grp.attrs['pml_parameter.imag'] * 1j


    def write_to_hdf5_group(self, hdf5_grp):
        """Export receiver's parameter values to an HDF5 dataset group."""
        hdf5_grp.attrs['N_sample_max'] = self.N_sample_max

        hdf5_grp.attrs['xlim'] = self.xlim
        hdf5_grp.attrs['Ngrid'] = self.Ngrid

        hdf5_grp.attrs['ninc'] = self.ninc
        hdf5_grp.attrs['nfar'] = self.nfar

        hdf5_grp.attrs['inc_param.n'] = self.inc_p['n']
        hdf5_grp.attrs['inc_param.app'] = self.inc_p['app']
        hdf5_grp.attrs['inc_param.cent'] = self.inc_p['cent']

        hdf5_grp.attrs['far_param.n'] = self.far_p['n']
        hdf5_grp.attrs['far_param.app'] = self.far_p['app']
        hdf5_grp.attrs['far_param.cent'] = self.far_p['cent']

        hdf5_grp.attrs['kappa'] = self.kappa
        hdf5_grp.attrs['porder'] = self.porder

        hdf5_grp.attrs['hmax_a'] = self.hmax_a
        hdf5_grp.attrs['pml_rad'] = self.pml_rad
        hdf5_grp.attrs['pml_delta'] = self.pml_delta
        hdf5_grp.attrs['pml_parameter.real'] = np.real(self.pml_parameter)
        hdf5_grp.attrs['pml_parameter.imag'] = np.imag(self.pml_parameter)

    def do_plot(self, i, image, m_approx):
        import matplotlib

        # The Agg renderer should be threadsafe etc.
        matplotlib.use('Agg')

        import matplotlib.pyplot as plt
        import cmocean as cmo

        n_matrix = image.reshape(self.Ngrid, self.Ngrid) + 1 # True Solution
        n_born = m_approx.reshape(self.Ngrid, self.Ngrid) + 1# Approximated Solution

        # Calculate the common min and max values for the color scale
        min_val = 1
        max_val = 2

        # Create the combined figure with two subplots
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

        # Plot the approximation of n(x) in the first subplot with a common color scale
        im0 = axs[0].imshow(
                np.real(n_born), 
                extent=(-1, 1, -1, 1), 
                origin='lower', 
                cmap='cmo.dense', 
                vmin=min_val, 
                vmax=max_val)
        axs[0].set_title(f'Born Approximation - Sample {i}')
        axs[0].set_xlabel('x')
        axs[0].set_ylabel('y')

        # Plot the true n(x) in the second subplot with a common color scale
        im1 = axs[1].imshow(
                np.real(n_matrix), 
                extent=(-1, 1, -1, 1), 
                origin='lower', 
                cmap='cmo.dense', 
                vmin=min_val, 
                vmax=max_val)
        axs[1].set_title(f'Ground Truth - Sample {i}')
        axs[1].set_xlabel('x')
        axs[1].set_ylabel('y')

        # Create a single colorbar for both plots with the common color scale
        cbar = fig.colorbar(im0, ax=axs, orientation='vertical', fraction=0.046, pad=0.04)
        cbar.set_label('Refractive Index Function Value')

        # Save the figure
        fig.savefig(os.path.join(self.plot_dir, f'testing_sample_{i}.png'), bbox_inches='tight')
        plt.close(fig)


#
# For multiprocessing to work properly (either from the multiprocessing module or from the
# concurrent futures) the code must be separated so that all work that happens in the original
# serial context is constrained to execute in the __main__ context.  This prevents the additional
# processes from running setup code when they are spawned.
#
# This is why the body of the computational loop has been turned into the function do_iteration():
# to keep it separated from the __main__ serial context.
#
if __name__ == '__main__':
    #
    # Build our command-line argument options and parse whatever the user
    # handed to this script:
    #
    cli_parser = argparse.ArgumentParser()
    cli_parser.add_argument('--debug', '-d',
        dest='is_debug',
        action='store_true',
        default=False,
        help='Enable debugging output')
    cli_parser.add_argument('--generate-plots', '-p',
        dest='should_plot',
        action='store_true',
        default=False,
        help='Enable generation of figure plots')
    cli_parser.add_argument('--data-file', '-f', metavar='<file-path>',
        dest='data_file',
        default='circ_data.hdf5',
        help='Path to the data file to which results are written (default: circ_data.hdf5)')
    cli_parser.add_argument('--plot-dir', '-D', metavar='<dir-path>',
        dest='plot_dir',
        default='./plots',
        help='Directory to which figure plots are written (default: ./plots)')
    cli_parser.add_argument('--base-index', '-i', metavar='<int>',
        dest='base_index',
        type=int,
        default=0,
        help='Initial result index in the data file')
    cli_parser.add_argument('--n-sample', '-n', metavar='<int>',
        dest='N_sample',
        type=int,
        default=100,
        help='Sample range in the data file')
    cli_parser.add_argument('--n-sample-max', '-N', metavar='<int>',
        dest='N_sample_max',
        type=int,
        default=20000,
        help='Total possible samples in the data file')
    cli_args = cli_parser.parse_args()

    #
    # Include the timestamp, source module and line number; max logging
    # level is INFO — which is what we'll use for most of our messsages.
    #
    if cli_args.is_debug:
        logging.basicConfig(format='%(asctime)s [%(process)d][%(levelname)s] %(message)s (%(module)s:%(lineno)d)', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(asctime)s [%(process)d][%(levelname)s] %(message)s', level=logging.INFO)

    #------------------------INITIALIZE PARAMETERS----------------------------------------------------
    if cli_args.base_index + cli_args.N_sample > cli_args.N_sample_max:
        cli_args.N_sample = cli_args.N_sample_max - cli_args.base_index

    params = Parameters(should_plot=cli_args.should_plot,
                        plot_dir=cli_args.plot_dir,
                        Idx_sample=cli_args.base_index,
                        N_sample=cli_args.N_sample,
                        N_sample_max=cli_args.N_sample_max)


    #------------------------INITIALIZE HDF5 DATA FILE------------------------------------------------
    #
    # We use no caching since this is a write-only situation and we don't want a lot of memory to be
    # occupied unnecessarily.
    #
    # The critical parameters are ninc, nfar, Ngrid, N_sample, and N_sample_max.  The latter dictates
    # the maximum attainable size of the datasets.  For a set of 20000 samples handled in batches of
    # 100, N_sample_max would be 20000 and N_sample would be 100, with Idx_sample varied according to
    # the batch number (0, 100, 200, ..., 19900).
    #
    # If the file already exists, parameters are imported from it.  If the file does not exist, then
    # as part of the initialization the current parameter values (in params) are written to the file.
    #
    # Best-effort locking is requested in case this code will be run using a file on our Lustre file
    # system -- which does not support the necessary locking functionality.  Rather than throw an
    # exception, H5Py will attempt to work around it.  This shouldn't be a problem since this code
    # is structured with a single worker able to write to the file even when parallelized.
    #
    if os.path.isfile(cli_args.data_file):
        data_file = h5py.File(cli_args.data_file,
                    mode='r+',
                    rdcc_nbytes=0, rdcc_nslots=0, rdcc_w0=1,
                    locking='best-effort')
        params.read_from_hdf5_group(data_file)
        logging.info(f'HDF5 data file {cli_args.data_file} opened')
    else:
        # Create HDF5 data file for results.
        #
        data_file = h5py.File(cli_args.data_file,
                    mode='w',
                    rdcc_nbytes=0, rdcc_nslots=0, rdcc_w0=1,
                    locking='best-effort')
        params.write_to_hdf5_group(data_file)
        logging.info(f'HDF5 data file {cli_args.data_file} created')
    #
    # Create or load a reference to the existing datasets in the file.  Note that HDF5 does not
    # natively handle complex floating point, so a real and imaginary dataset must be created to
    # hold each half of the farfield data.
    #
    # All data are stored as native 64-bit (double-precision) floating point values.  The shape
    # starts out as a minimum (*, 1) out of (*,N_sample_max) but is immediately resized to match
    # the index and range of the current run.
    #
    farfield_real_ds = data_file.require_dataset(
                                'farfield.real', exact=True,
                                shape=(params.ninc*params.nfar, 1),
                                chunks=(params.ninc*params.nfar, 1),
                                maxshape=(params.ninc*params.nfar, params.N_sample_max),
                                dtype=np.float64, compression='gzip')
    if farfield_real_ds.shape[1] < params.Idx_sample + params.N_sample:
        farfield_real_ds.resize((params.ninc*params.nfar, params.Idx_sample + params.N_sample))
    farfield_imag_ds = data_file.require_dataset(
                                'farfield.imag', exact=True,
                                shape=(params.ninc*params.nfar, 1),
                                chunks=(params.ninc*params.nfar, 1),
                                maxshape=(params.ninc*params.nfar, params.N_sample_max),
                                dtype=np.float64, compression='gzip')
    if farfield_imag_ds.shape[1] < params.Idx_sample + params.N_sample:
        farfield_imag_ds.resize((params.ninc*params.nfar, params.Idx_sample + params.N_sample))

    approx_ds = data_file.require_dataset('approx', exact=True,
                                shape=(params.Ngrid**2, 1),
                                chunks=(params.Ngrid**2, 1),
                                maxshape=(params.Ngrid**2, params.N_sample_max),
                                dtype=np.float64)
    if approx_ds.shape[1] < params.Idx_sample + params.N_sample:
        approx_ds.resize((params.Ngrid**2, params.Idx_sample + params.N_sample))

    image_ds = data_file.require_dataset('image', exact=True,
                                shape=(params.Ngrid **2, 1),
                                chunks=(params.Ngrid **2, 1),
                                maxshape=(params.Ngrid **2, params.N_sample_max),
                                dtype=np.float64)
    if image_ds.shape[1] < params.Idx_sample + params.N_sample:
        image_ds.resize((params.Ngrid **2, params.Idx_sample + params.N_sample))

    data_file.flush()

    #
    # Summarize the current parameterization to stdout:
    #
    params.summarize()

    #
    # Make sure the plot directory exists if that is requested:
    #
    if params.should_plot:
        if not os.path.isdir(params.plot_dir):
            if os.path.exists(params.plot_dir):
                logging.error('Plotting directory %s exists but is not a directory', params.plot_dir)
            else:
                os.mkdir(params.plot_dir)
                logging.info('Create plotting directory %s', params.plot_dir)

    #
    # Allocate a multiprocessing pool of (by default) 1 process or however many CPUs Slurm
    # says it has allocated to us on this node.  This allows the job to easily adapt to
    # whatever --cpus-per-task was associated with the job submission -- without having to
    # modify this Python code or pass arguments to it.
    #
    n_workers = int(os.getenv('SLURM_CPUS_ON_NODE', 1))
    n_workers = 2 # Modify according to CPU power
    logging.info(f'Entering computational loop with {params.N_sample} sample(s) across {n_workers} worker(s)')

    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            #
            # We wish to have the worker processes step through iteration indices in the range
            # [0, N_sample).  The primary process will concurrently wait for the multiprocessing
            # pool to buffer and return the results of those function calls so that the data
            # can be written to the HDF5 file.
            #
            # Since the do_iteration function takes two arguments -- the iteration index and
            # a reference to the params object -- a "starmap" must be used.  The list of indices
            # and a list of N_sample params object references are merged into a single list of
            # tuples (index_i, params) using zip().
            #
            for (i, image, m_approx, uinf) in pool.map(do_iteration, zip(range(params.N_sample), [params]*params.N_sample)):
                # Write results to disk:
                file_idx = i + params.Idx_sample
                image_ds[:, file_idx] = image
                logging.info('    Wrote image[%d] to disk', file_idx)

                approx_ds[:, file_idx] = np.real(m_approx)
                logging.info('    Wrote approx[%d] to disk', file_idx)

                farfield_real_ds[:, file_idx] = np.real(uinf)
                farfield_imag_ds[:, file_idx] = np.imag(uinf)
                logging.info('    Wrote farfield[%d] to disk', file_idx)

                # Plot if necessary -- WARNING, this segfaults for me!
                if params.should_plot:
                    params.do_plot(file_idx, image, m_approx)

                # Purge any orphaned objects to keep memory usage under control:
                gc.collect()
    else:
        #
        # This is just a serial version of the process pool loop above; there's no reason to
        # add the overhead of multiprocess communications etc. for a serial execution:
        #
        for i in range(params.N_sample):
            # Purge any orphaned objects to keep memory usage under control:
            gc.collect()

            # Generate the next results set:
            (i, image, m_approx, uinf, epsilon, tau, uinf_born)= do_iteration((i, params))

            # Write results to disk:
            file_idx = i + params.Idx_sample
            image_ds[:, file_idx] = image
            logging.info('    Wrote image[%d] to disk', file_idx)

            approx_ds[:, file_idx] = np.real(m_approx)
            logging.info('    Wrote approx[%d] to disk', file_idx)

            farfield_real_ds[:, file_idx] = np.real(uinf)
            farfield_imag_ds[:, file_idx] = np.imag(uinf)
            logging.info('    Wrote farfield[%d] to disk', file_idx)

            # Plot if necessary -- WARNING, this segfaults for me!
            if params.should_plot:
                params.do_plot(file_idx, image, m_approx)

    logging.info('Computational loop completed.')

    log_memory_usage()
