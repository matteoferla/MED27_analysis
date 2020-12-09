## Model attempts

Phyre2 and ITasser models were made and energy minimised with Pyrosetta (1 cycle with `ref2015_cart` and 15 with `ref2015`).


| model           |   total_score |   fa_atr |   fa_rep |   fa_sol |   fa_intra_rep |   fa_elec |   fa_dun |   rama_prepro |   omega |
|----------------:|--------------:|---------:|---------:|---------:|---------------:|----------:|---------:|--------------:|--------:|
| weights (ref2015)         |              |     1.0  |  0.55  |  1.0  |  0.005  |  1.0  |  0.7  |  0.45  |  0.4   |
|  MED27_ITasser_4 |        -666.9 |  -1826.4 |    427.8 |   1208.5 |          738.7 |    -581.2 |    528.1 |         240.6 |    88.1 |
|  MED27_ITasser_1 |        -614.5 |  -1942.7 |    481.9 |   1303.5 |          700.6 |    -604.1 |    545.2 |         402   |   102.6 |
|  MED27_ITasser_2 |        -585   |  -1883.9 |    438.7 |   1278.2 |          725.9 |    -593.4 |    545.7 |         404.3 |   122.9 |
|  MED27_ITasser_3 |        -546   |  -1825.9 |    439.5 |   1193.2 |          758.9 |    -534.2 |    515.9 |         380.2 |   136.7 |
|  MED27_ITasser_5 |        -537.9 |  -1741.2 |    412.3 |   1125.7 |          846.3 |    -503.7 |    526.1 |         289   |   107.7 |
|  MED27_phyre     |        -318.6 |  -1626.9 |    415.6 |   1040.2 |         1187.1 |    -451.6 |    525.9 |         562.1 |    56.7 |


The top ITasser model has a projected RSMD from a hypothetical solution of 15 Å (bad).
It has a slightly higher ∆G than ITasser model 4, 
but the attractive term of the Lennard-Jones equation and the electrostatic term are better.

### CryoEM
The only structure available is the low resolution _S. pombe_ highly disordered fragment, which is binding to MED17 and MED22 (head)
and not MED14, the stalk.

![MED27_5U0P](MED27_5U0P.jpg)

The homology to 5U0P_2 was so low it wasn't a big player in the models. None of the models resolved remotely similarly to it.
ITasser model 3 has two antiparallel sheets, not 4 and the numbering does not match, nor the helices.
So this is best ignored, although the real structure could be similar (or the density is parts of MED2+MED3+MED15)

### XL-MS
Aligning `MED3_YEAST` and `MED27_HUMAN` with Muscle gives the following:

    MED3  MDSIIPAGVKLDDLQVILAKNENETRDKVCKQINEARDEILPLRLQFNEFIQIMANIDQEGSKQADRMAKYLHIRDK--ILQLNDRFQTLSSHLEALQPLFSTVPEYLKTADNRDRSFQLLEPLSTYNKNGNAVCSTATVVSTNHSAAASTPTTTATPHANPITHAHSLSNPNSTATMQHNPLAGKRGPKSGSTMGTPTVHNSTAAAPIAAPKKPRKPRQTKKAKAQAQAQAQAQAQVYAQQSTVQTPITASMAAALPNPTPSMINSVSPTNVMGTPLTNMMSPMGNAYSMGAQNQGGQVSMSQFNGSGNGSNPNTNTNSNNTPLQSQL-NLNNLTPANILN--MSMNNDFQQQQQQQQQ--QQQPQPQYNMNMGMNNMNNGGKELDSLDLNNLELGGLNMDFL
    MED27 MADVINVSVNLEAF-------------------SQAISAIQALRSSVSRVFDCL----KDGMRNKETLEG----REKAFIAHFQDNLHSVNRDLNELERLSNLVG---KPSENH--------PLH------------------NSGLLSLDPVQDKTPLYSQLLQAYKWSN-----KLQYHA-----GLASG-LLNQQSLKRSANQMGVSAKRRPKAQPTTLVLPPQYVDDVISRIDRMFPEMSIHLSRPNGTSAMLLVTLGKVLKVIV---VMRSLFIDRTIVKG--YNENVYTEDGKLDI--WSKSNYQVFQKVTDHATTALLHYQLPQMPDVVVRSFMTWLRSYIKLFQAPCQRCGKFLQDGLPPTWRDFRTLEAFHDTCRQ-------------------

Using this and the data from [the yeast XL-MS dataset](https://www.cell.com/cell/fulltext/S0092-8674(16)31147-3),
allows some tentative tests. This assumes that the sections are structurally conserved and the pairwise alignment is correct
—which it does not.

Of note is that the MED3 residue 20–27 interaction had to be removed as as these flank a gap (both correspond to residue 14)
None of the models have a loop at residue 14 —ITasser 5 has a unstructured loop as opposed to a first helix. Whereas the residue 14 on a loop link is consistent with 5U0P_2 two short helices.
They all have a straight helix in the models.
Here are the distances of the CA atoms in Å.


| name                        |   10-79 |   40-50 |   40-47 |   50-47 |   50-124 |
|-----------------------------:|--------:|--------:|--------:|--------:|---------:|
|  MED27_ITasser_4 |    12   |     8.8 |    11.1 |     4.9 |     30.3 |
|  MED27_ITasser_1 |    17.8 |    17.4 |    15.5 |     8.2 |     30.2 |
|  MED27_ITasser_2 |    10.4 |    13.4 |    14.6 |     5.3 |     12.9 |
|  MED27_ITasser_3 |    14.3 |     8.5 |     9.4 |     5.6 |     10.3 |
|  MED27_ITasser_5 |    35.9 |    17.4 |    16.4 |     5.1 |     45.4 |
|  MED27_phyre    |    13.2 |    16.3 |    13.5 |     7   |     56.5 |

In terms of distances, the residue pairs 10–79 and 50–124 are distant.
The former pair, although closer in ITasser4, is problematic in both ITasser 4 and ITasser 1 (core residue pointing away).
Same for the latter.

Therefore, this dataset is not helpful. Consequently going with ITasser's top ranked pose is the best option.