# Temporal integration is a robust feature of perceptual decisions

This repository contains the Matlab code to carry out the data analysis for the manuscript entitled "Temporal integration is a robust feature of perceptual decisions" by A. Hyafil, J. de la Rocha, C. Pericas,  LN. Katz,  AC. Huk & JW. Pillow ([see preprint](https://www.biorxiv.org/content/10.1101/2022.10.25.513647v1)).

Please report any bug / missing function.

## Set-up
- add the directory with subfolders to your Matlab path
- add also the [GUM toolbox](github.com/ahyafil/gum/) to your Matlab path
- edit `TemporalIntegrationRootDirectory.m` to add the path to local directory

## Data

All datasets to reproduce the results (behavioral data in monkeys, humans and rats; LIP spike count in monkeys) are provided as csv files in `data` folder. There is one row for each trial. The variables in the csv are straightforward:
- `subject` is the subject ID
- `session` is the session ID (for monkey and rat data, the human experiment is single session
- `target` denotes the target response (i.e. correct response). 0: left response, 1: right response
- `resp` is the actual response of the subject
- `stimulus_i` is the stimulus evidence of the $i$th sample in the stimulus sequence. This is changed to vector form using `preprocess_table.m`.

## Usage

The script `temporal_integration` produces all model-fitting and complementary analysis (i.e. disagree trials). Edit the first lines of code to select which analysis to run, on what dataset (monkey/human/rat), what variant of the non-integration model, etc. All results are stored in the `modelfits` directory, which is already populated (so you can produce the figures without running these analyses).

Figures can then be produced using the following scripts:

| Figure | Script | 
|---------------|---------------|
|Fig. 2 | /compare_model_fits.m |
|Fig. 3 | /plot_disagree_trials_figure.m |
|Fig. 4 | /plot_integration_map_figure.m |
|Fig. 5 | /compare_model_fits.m |
|Fig. 6 | /compare_model_fits.m |
|Fig. 2 supp 1 | /supplementary/plot_model_parameters.m |
|Fig. 2 supp 2 | /supplementary/compare_snapshot_models.m |
|Fig. 2 supp 3 | /supplementary/compare_extrema_detection_models.m |
|Fig. 3 supp 1 | /supplementary/plot_subjective_weights.m |
|Fig. 4 supp 1 | /plot_integration_map_figure.m |
|Fig. 4 supp 2 | /supplementary/plot_integration_maps_simulated_data.m |
|Fig. 4 supp 3 | /supplementary/temporal_integration_LIP.m |
|Fig. 5 supp 1 | /compare_model_fits.m |
|Fig. 6 supp 1 | /compare_model_fits.m |

