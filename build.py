import os
from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from hatchling.build import build_wheel

class CustomBuildHook(BuildHookInterface):
    def initialize(self, version, build_data):
        build_data['force_include'] = {
            'build/bin/EpiSim': 'epi_sim/EpiSim'
        }
        return super().initialize(version, build_data)

# Remove the build function and directly use a dictionary for configuration
config = {
    "project": {
        "name": "epi_sim",
        "version": "0.1.0",
        "description": "A wrapper for EpiSim.jl",
        "readme": "README.md",
        "requires-python": ">=3.8",
        "license": "BSD-3-Clause",
        "authors": [
            {
                "name": "Lewis Knox",
                "email": "lknoxstr@bsc.es"
            }
        ],
        "dependencies": [
            "pandas"
        ],
        "classifiers": [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: BSD License",
            "Operating System :: OS Independent",
        ],
    },
    "build-system": {
        "requires": ["hatchling"],
        "build-backend": "hatchling.build"
    },
    "tool": {
        "hatch": {
            "build": {
                "hooks": {
                    "custom": {
                        "path": __file__,
                        "class": "CustomBuildHook"
                    }
                }
            }
        }
    }
}

if __name__ == "__main__":
    build_wheel(wheel_directory="dist", config_settings=config)