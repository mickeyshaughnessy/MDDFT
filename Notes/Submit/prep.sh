#!/bin/bash
stem="db_aimd"
tex=`ls ../${stem}.tex`
echo "found $tex"
bbl=`ls ../${stem}.bbl`
echo "found $bbl"
echo "combining $tex and $bbl"
sed 'N;$!P;$!D;$d' $tex | sed "s/Figures/./g" > tmp
echo "\\end{document}" > end
cat tmp $bbl end > ${stem}.tex

raw=`grep includegraphics $tex`
i=1
for str in $raw; do
  str=${str##\{*\{}
  str=${str%\}*\}}
  str=${str##*\/}
  echo "$i linking  $str"
  i=$((i+1))
  cp -f  ../Figures/$str .
done

tex=${stem}.tex
echo "TeXing once"
pdflatex $tex > tex.stdout
echo "TeXing twice"
pdflatex $tex > tex.stdout

echo "use preview with \" reduce filesize \" " 
