from setuptools import setup, find_packages

setup(
    name='episim',
    version='0.1.0',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'pandas',
        'numpy',
        'shutil',
        'uuid',
        'logging',
        'subprocess',
        'json',
    ],
    author='Lewis Knox',
    author_email='lknoxstr@bsc.es',
    description='A wrapper for EpiSim.jl',
    url='https://github.com/Epi-Sim/EpiSim.jl.git',  # Update with your repo URL
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)