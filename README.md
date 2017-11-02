# bedtack

The name was chosen principally because it "tackles" BED files and the name BEDtools and BEDutils imply a much bigger platform of tools while this project is rather small.

The name also has a fuzzy evocation of "BED attack" or "BED tuck" which all imply a certain handling of the BED file format which actually varies quite a bit, principally on the number of columns it has.

If you do not have a BED file, but rather have a GFF file, there are tools out there to do the conversion, so it's not part of the current scope.

## 0 indexing
GFF and BED files use 0 indexing, as does (by default) the size files.
This means semi-open intervals. The first number is inclusive of the range, while the second is one above the final included position.

## chromosome order 4, 9, 5
Lexicographic ordering means 4,9,5 because Roman numerals are used for chromosome names
