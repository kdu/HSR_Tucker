% HSR_Tucker: a collection of methods for the HSR problem
%
% Files
% test_exp1               - runs SCOTT for various ranks, saves files and figures
% test_exp2               - experiment on the recovery of spectral signatures
% test_exp3_pavia         - compares blind algorithms for Pavia University
% test_exp3               - compares non-blind algorithms for Salinas/IP
% test_exp4               - same as test_exp1 for the pansharpening problem
% test_exp3_ku            - runs BSCOTT for various ranks, Indian Pines
% bscott                  - runs BSCOTT algorithm
% scott                   - runs SCOTT algorithm
% scott_opti              - runs SCOTT algorithm with TensorLab optimization
% scuba                   - runs SCUBA algorithm
% stereo                  - runs STEREO algorithm
% stereo_init             - initializes the factor matrices for STEREO
% cc                      - computes Pearson cross-correlation btw two tensors
% compute_metrics         - returns array of metrics
% ergas                   - computes relative global error btw two tensors
% nmse                    - computes NMSE btw two tensors
% r_snr                   - computes R-SnR btw two tensors
% sam                     - computes spectral angle mapper btw two tensors
% crop                    - crops array to the specified dimensions
% gauss_kernel            - computes Gaussian blurring kernel
% solveC_normal           - solves C in stereo_init with normal equations
% solveC_qr               - solves C in stereo_init with QR factorization
% spatial_deg             - computes spatial degradation matrices from SRI
% spectral_deg            - computes spectral degradation matrix from SRI
