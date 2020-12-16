# MED27_analysis
Analysis of variants of MED27 (no mutant data is stored here)

MED27 is part of the mediator complex. This is a huge complex.
![complex](complex.png)

## Mammalian MED27

A [preprint](https://www.biorxiv.org/content/10.1101/2020.10.05.326918v1.full) has the full mammalian complex.

### Humanise
It is mouse, but there are only a few missenses with human for MED27 (L5I, S6N, G8S, N55H, S220N, S235N).

I tried to convert the whole complex automatedly, but got bogged down and never finished: see [humanise](humanise.md).
Especially given how minor the differences are.

### Hydrate
Being a big complex of small protein it is very solvent exposed.
This is problematic for calculating interface energy when using implicit solvent —cf. barnase-barnstar, G-protein &alpha;-&beta;.
Hydration however took too long (does not appear to work linearly) and the variants of interest were not interface residues.

## Model building

Prior to finding the preprint, I [attempted to make a model](model_attempts.md).

## Consurf

[Consurf](https://consurf.tau.ac.il/) was then used to add conservation as b-factors.

## gnomAD

gnomAD has weird transcript issues again. The correct transcript is used in v3: `MED27-201`/`ENST00000292035.10`, which is 311.
While the longer (326aa) `MED27-207`/`ENST00000651950.1` is wrong.
There is a single homozygous variant, p.Val203Ile, 
and only a few heterozygous variants in the `gnomad_r3_controls_and_biobanks` dataset
(p.Ala2Val, p.Ile5Val, p.Ser26Phe, p.Lys36Arg, p.Ala150Thr, p.Pro174Ala, p.Val203Ile, p.Pro280Leu)
that a frequency greater than 5e-5.
