#!/bin/bash

# Check for conda installation
if ! command -v conda &> /dev/null
then
    echo "Conda not found. Please install conda before running this script."
    exit
fi

# Create environment from the YAML file
echo "Creating the 'FocalSV' conda environment from environment.yml..."
conda env create -f environment.yml

# Check if environment was successfully created
if [ $? -eq 0 ]; then
    echo "Environment 'FocalSV' created successfully."
else
    echo "Failed to create environment 'FocalSV'. Please check for errors."
    exit 1
fi

# Activate the environment
echo "Activating the 'FocalSV' environment..."
conda activate FocalSV

echo "Installation complete. You can now use the 'FocalSV' environment."
