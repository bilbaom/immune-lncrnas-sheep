#Bedtools bedGraph production 
#!/bin/bash
for f in ./*.bam;do
    NAME=$(basename -s .bam $f)
    bedtools genomecov -ibam $f -d -strand + | awk -v width=1 '!($1~/^NW/)&&($3!=0) {print $1,$2,$2+width,$3}' > ${NAME}.plus.bedGraph &
    bedtools genomecov -ibam $f -d -strand - | awk -v width=1 '!($1~/^NW/)&&($3!=0) {print $1,$2,$2+width,$3}' > ${NAME}.minus.bedGraph
done;

#BedGraphToBigWig (UCSC tools)
#!/bin/bash
for f in ./*.plus.bedGraph;do
	NAME=$(basename -s .plus.bedGraph $f)
	sort -k 1,1 -k2,2n ${f} > tmp_plus &&
	bedGraphToBigWig tmp_plus chromsizes.txt ${NAME}.plus.bw
	rm -f tmp_plus
done;

for f in ./*.minus.bedGraph;do
	NAME=$(basename -s .minus.bedGraph $f)
	sort -k 1,1 -k2,2n ${f} > tmp_minus &&
	bedGraphToBigWig tmp_minus chromsizes.txt ${NAME}.minus.bw
	rm -f tmp_minus
done;
