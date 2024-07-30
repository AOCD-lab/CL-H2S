# High-Throughput Computational Screening of Metal Sulfides for the Chemical Looping Elemental Decomposition of H2S
<p align="center">
  <img src="https://github.com/user-attachments/assets/654bfc6e-cfd4-426b-907f-3824d9051b99" alt="image">
</p>

This repository contains scripts and files related to constructing sulfide databases and their analysis for the chemical looping decomposition of H2S.

The research is organized into five directories, each dedicated to a specific task:

## 01_MP_sulfide_database
- Script for retrieving mono- and bimetallic sulfide databases from Materials Project (MP) with predicted Gibbs free energies.
- Mono- and bimetallic sulfide databases

## 02_MP_sulfides_analysis
- Scripts for constructing redox sulfide pairs and performing thermodynamic analysis.
- `.json` files with detailed information about each pair, including the chemical loop, free energy changes at temperatures ranging from 400K to 1200K, Gibbs formation energies, and formation enthalpies.
- Datasets of material pairs from limiting reaction analysis, sorted by the lowest Gibbs limiting reaction energy.

## 03_BOWSR_error
- Scripts for benchmarking formation enthalpies and decomposition energies using the MEGNet/BOWSR approach against MP values.
- Formation enthalpies is behcnmarked for bimetallic sulfides, oxides, nitrides, and hydrides.
- Datasets of corresponding materials used for benchmarking.

## 04_Element_substitution
- A script for elemental substitution to create new bimetallic sulfide structures.
- Scripts for running BOWSR/MEGNet geometry optimization and energy prediction.
- A collection of scripts for post-processing the results of optimization.

## 05_Expansion_analysis
- Thermodynamic analysis of newly predicted materials and the redox pairs derived from these materials.

## Authors
- Polina Tolstova
- Luigi Cavallo
