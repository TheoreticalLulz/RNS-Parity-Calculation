# RNS-Parity-Calculation
Residue number systems present a method for representing large numbers as a system of smaller numbers, so that arithmetic computation may be performed more efficiently. Though there are numerous methods of accomplishing this task, none are more generic or deceptively simple than the system inherent to the Chinese Remainder Theorem. By constructing an arbitrary set of coprime moduli, one may distinctly represent any number less than the modulus product as a mutually independent set of residues. However, while arithmetic operations such as multiplication and addition are greatly improved, other operations such as parity calculation are remarkably difficult.

Within this repository, three different parity calculation techniques are compared. These techniques include the brute force approach, in which the number is converted to traditional decimal format, the Mi Lu & Jen-Shuin Chiang algorithm, which utilizes pre-computed hash tables for each modulus, and an algorithm of my own design. My algorithm presented the fastest method of parity calculation within the generic residue number system at the time, and is further evaluated within the document "RNS-Parity-Calculation.pdf", a report created as part of a class on the design and analysis of algorithms. Notice that the algorithm, evaluation, and code were all created by me.
