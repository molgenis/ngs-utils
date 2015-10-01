
module load bedtools/2.19.1
module list

HOME="/gcc/groups/gcc/tmp01/rkanninga/vcf_compare/"
TRUTH="NextGene"
COMPARE="3.1.2"

THISDIR=`pwd`

if [ ! -d $HOME/$TRUTH/intersect/ ]
then
	mkdir $HOME/$TRUTH/intersect/
fi

if [ ! -d $HOME/$COMPARE/intersect/ ]
then
	mkdir $HOME/$COMPARE/intersect/
fi

cd ${HOME}/${TRUTH}
for i in $(ls *.vcf)
do
bedtools intersect -a $i -b ${HOME}/captured.bed > $HOME/$TRUTH/intersect/${i}.intersect.vcf
echo "$i done"
done
echo ""
echo "TRUTH done"
echo ""

cd ${HOME}/${COMPARE}
for i in $(ls *.vcf)
do
bedtools intersect -a ${i} -b ${HOME}/captured.bed > $HOME/$COMPARE/intersect/$i.intersect.vcf
echo "$i done"
done

echo ""
echo "COMPARE done"
