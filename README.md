advection1d
===========

Getting the source
------------------

Clone from the git repo::

    git clone git@github.com:johnomotani/advection1d.git

Compiling
---------

### Option 1 - use system libraries

First install the BLAS and FFTW3 libraries on your system.

Then compile by running:

    ./compile_link_system_blas-fftw.sh

### Option 2 - download and compile BLAS and FFTW

[This was intended to test if setting extra optimisation options to FFTW and/or BLAS could improve performance. Did not see any improvement on my Linux Mint system /JTO]

Run the script:

    ./compile.sh

This will create sub-directories for BLAS, CBLAS and FFTW (under the ``external`` directory), download them and compile them; finally it will compile advection1d, linking to the locally installed BLAS and FFTW.

Running
-------

Run the executable ``advection1d``. It is suggested to add the ``advection1d`` directory to your ``$PATH`` to make this simple. Otherwise use the full (absolute or relative) path to the executable file.

To change the settings, create a file called ``input.toml`` in the directory where you want to run. See ``parameters.hxx`` for the possible parameters. The file is in [TOML](https://toml.io/en/) format: each line is in the format ``key = value``. Note that strings must be contained in quote marks and floats must have a number after the point, e.g. ``0.0`` not ``0.``.
