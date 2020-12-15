# MED27_analysis
Analysis of variants of MED27 (no mutant data is stored here)

## State of the art

MED27 is part of the mediator complex. This is a huge complex.
![complex](complex.png)

A lot of research has gone into this, as [summarised in a recent review](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6393861/).

Most structures are from either brewer's year or fission yeast, none from humans.

* There are several parts ('modules'), the head, the middle and the tail, the stalk (MED14) and the CDK8 module.
* The stalk connects the head, middle and tail.
* The head is connected to the top of the middle so actually is on the side.
* The head may or may not join to the tail (MED17 and MED27)
* The tail is composed in brewer's yeast by MED2, MED3, MED5, MED15 and MED16.
* The [potential homologues](https://academic.oup.com/nar/article/36/12/3993/1135622) MED29, MED27 and MED24 take the role of MED2, MED3 and MED5 in metazoa
* Humans also have additionally MED25 in the tail
* Humans also have additionally MED26 in the middle
* Humans also have additionally MED27, MED28, MED29 and MED30 in the head
* Humans also have two copies of MED12 (MED12/MED12L) and MED13 (MED13/MED13L), which are specific for CDK8/19
* MED2, MED3 and MED15 form a subgroup as seen [XL-MS](https://www.cell.com/cell/fulltext/S0092-8674(16)31147-3) and [cyroEM](https://www.cell.com/cell/fulltext/S0092-8674(14)00609-6).
* In the [cyroEM](https://www.cell.com/cell/fulltext/S0092-8674(14)00609-6) MED3/MED27+MED2/29+MED15 form an L shape that binds to MED16 and MED14.
* There is a [4.4 Å cryoEM model](https://www.nature.com/articles/nature21393) (PDB:5U0P) with three short chains of MED3, but MED2 and MED15 aren't resolved as aren't other members of the tail.
* There are no available human nuclear cross-linking mass spec studies with a resolution sufficient to find contacts between MED29, MED27 and MED15.
* As human MED29 is much smaller than yeast MED2 and human MED27 larger than yeast MED3, it cannot be said whether MED29, MED27 and MED15 from an L shape
* the MED29, MED27 and MED15 subcomplex may have one lone helix forming a bundle with a helix from another, therefore any theoretical model may be biased.
* the XL-MS interaction in MED3 may, or may not be accurate for MED27.

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
