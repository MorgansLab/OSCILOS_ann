# OSCILOS-ann
## What is OSCILOS ?
Combustion instabilities are generally not desirable because they may lead to an early ageing of the combustion chamber or even to severe structural damage. Knowledge of this complex mechanism is necessary in the development of control strategies. The open source combustion instability low-order simulator (OSCILOS) is an open source code for simulating combustion instability. It is written in Matlab / Simulink and is very straightforward to run and edit. It can simulate both longitudinal and annular combustor geometries. It represents a combustor as a network of connected modules. The acoustic waves are modeled as either 1-D plane waves (longitudinal combustors) or 2-D plane/circumferential waves (annular combustors). A variety of inlet and exit acoustic boundary conditions are possible, including open, closed, choked and user defined boundary conditions. The response of the flame to acoustic waves is captured via a flame model; flame models ranging from linear n-tau models to non-linear flame describing functions, either prescribed analytically or loaded from experiment / CFD data, can be prescribed. The mean flow is calculated simply by assuming 1-D flow conditions, with changes only across module interfaces or flames. This current version is for annular modes. 

## How to contribute
For guidelines on contributing and reporting issues see [CONTRIBUTING.md](CONTRIBUTING.md).

## Licensing
OSCILOS is freely available under an [open source license](LICENSE.md).