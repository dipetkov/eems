
## FST matrix as a distance matrix for the .diffs input file

This section is in response to the following question from Àlex Mas, PhD student at David Coma's lab, in the Institute of Evolutionary Biology, in Barcelona:

**I wonder which conditions the matrix must accomplish and how a FST matrix could work as the .diffs file.**

(Thank you for the interesting question.)

It is not appropriate to plug in the matrix of pairwise $F_{ST}$s for the matrix of average genetic differences `diffs`. Instead, I have written a modified version `runeems_Fsts`.

* `FSTs-as-diffs.pdf` are my notes on `runeems_Fsts`. These serve as documentation, so read them first. The notes also explain why it was necessary to modify the original implementation.
* `runeems_Fsts/data` contains two datasets I have used for testing.
* `runeems_Fsts/src` contains the modified EEMS program to use with the $F_{ST}$ matrix.
* `runeems_Fsts/src0` contains an earlier version of the modified EEMS program, included here for completeness.
