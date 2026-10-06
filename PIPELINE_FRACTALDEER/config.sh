#!/bin/bash
conda env update --file ./environment_fractalDeer_PRD.yml --prune
~/anaconda3/envs/pipeline_fractalDeer_prd/bin/python -m ipykernel install --user --name pipeline_fractalDeer_prd --display-name "Conda (pipeline_fractalDeer)"

echo "Configuration terminée avec succès !"
