python GetAnnotation.py --input QTL.regions.list

echo "Draw only peak interval"
perl GetAnnotation_Peak.py --input QTL.regions.list.peak > log 2> log2 &

echo "Draw 1.5 LOD interval"
perl GetAnnotation_Region.py --input QTL.regions.list

