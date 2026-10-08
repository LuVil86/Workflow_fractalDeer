# Pipeline fractalDeer

## Description

This repository contains the scripts and functions to run a full fractal analyses of GPS trajectory

## Getting Started

In order for the pipeline to run, you need to set up the conda environment that is specified in environment_fractalDeer_PRD.yml . You can use the following command :

```         
conda env create --file ./environment_fractalDeer_PRD.yml 
```

then activate the environment

```         
conda activate pipeline_fractalDeer_prd
```

If you want to use the jupyter-notebook demo, you need to activate the environment and build the ipykernel within the conda environment. NOTE the path to conda environements are probably different in your machine, so change accordingly

```         
python -m ipykernel install --user --name pipeline_fractalDeer_prd --display-name "Conda (pipeline_fractalDeer)"
```
then run jupyter-notebook
```         
jupyter-notebook ./PIPELINE_FRACTALDEER/interactive_demo.ipynb
```


### Dependencies

- pipeline is (should be ..) fully operationnal if you use the conda environment provided.

### Installing

- How/where to download your program
- Any modifications needed to be made to files/folders

### Executing program

- How to run the program
- Step-by-step bullets

```         
code blocks for commands
```

## Help

Any advise for common problems or issues.

```         
command to run if program contains helper info
```

## Authors

Lucas Villard : coder and maintainer

## Version History

- 0.2
  - Various bug fixes and optimizations
  - See [commit change]() or See [release history]()
- 0.1
  - Initial Release

## License

No license for the moment

## Acknowledgments
