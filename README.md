# statistical-thermodynamics
Project written in R developed for the statistical thermodynamic analysis of protein stability in complex conditions

## What are Biologics?
Biologics are a class of drugs that are derived from natural sources whose pharmaceutical component is too complex to be produced _ _en masse_ _ by organic synthesis
A lot of these drugs have a protein component and, as such, their stability is affected under different chemical conditions, i.e the stabilising or destabilising potentials of surrounding cosolvents can modify the point at which a protein will change its conformation.

## How do we measure the (de)stabilising effects of a cosolvent?
In the past, different models have been produced to provide insight into how different chemicals affect a protein's stability. However, up to this point, none had been able to model stabilising effects with regards to what actually occurs at the molecular level. And if they could model at the molecular level, they could only do so for binary systems of solvent-cosolvent systems (i.e one cosolvent and protein in a solution).
A statistical thermodynamic approach, developed by Prof. Seishi Shimizu at the University of York and Prof. Nobuyuki Matubayasi of Osaka University, models the stabilising effects of **multiple** cosolvents on a protein's stability in a ternary system.

## Okay, so how is that useful?
We can use the mathematical model to make predictions on the stability of a protein as long as we know the chemical potential and concentrations of cosolvents. This allows us to simulate the chemical conditions of a system and apply it to predict the stability of novel biopharmaceuticals.
Conditions for the uses of biologics vary widely from highly acidic conditions in the gut biome and very hypoxic conditions in cancers and tumours, to oxidising conditions of muscular tissue. Pharmaceuticals need to be stable under any condition that a target site is in, as well as the conditions of any entry points if they are ectopic to the target.

## So what have I done?
I've developed a script in R that analyses data drawn from previous experiments to extract the stabilising and destabilising effects of a cosolvent on a protein, as well as the melting temperature of a protein under the conditions produced by cosolvents
This script takes in .csv data files and processes them to produce graphs and data illustrating the effects of cosolvents on protein stability
