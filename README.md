# HSR_Tucker: Hyperspectral Super-Resolution with Coupled Tucker Approximation

This repository contains the implementation for the methods proposed in the preprint:
```
@unpublished{prevost:hal-01911969,
 TITLE = {{Hyperspectral super-resolution with coupled Tucker approximation: Identifiability and SVD-based algorithms}},
 AUTHOR = {Pr{\'e}vost, Cl{\'e}mence and Usevich, Konstantin and Comon, Pierre and Brie, David},
 URL = {https://hal.archives-ouvertes.fr/hal-01911969},
 YEAR = {2019},
 HAL_ID = {hal-01911969},
 HAL_VERSION = {v3},
}
```
It also reproduces the results of the experiments from the preprint.

## Installation 
[TensorLab v.3.0](https://www.tensorlab.net) should be installed and should be on the MATLAB path.

The datasets used in the examples can be downloaded [here](http://www.ehu.eus/ccwintco/index.php/Hyperspectral_Remote_Sensing_Scenes).

For some of the experiments, the following packages need to be downloaded:
 * [Codes for the STEREO method](https://github.com/marhar19/HSR_via_tensor_decomposition) for some  comparisons with comparisons with  the STEREO methods
 * [HySure package](https://github.com/alfaiate/HySure/tree/master/src) for comparisons with HySure
 * [Export Fig](https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig) to save some figures in PDF

## Contents
 * `src` contains the source codes for the methods
 * `demos` contains the codes for the reproducivle examples
 
