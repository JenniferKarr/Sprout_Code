
FILENAME=$1
count=0

while read LINE
do
    let count++
    new=$(printf "%04d.jpg" ${count}) #04 pad to length of 4
    \cp ${LINE} ${new}  #COPY TO NEW FILE
done < $FILENAME

convert -antialias -delay 1x4 ????.jpg $2  #MAKE MOVIE
\rm ????.jpg  #REMOVE TEMPORARY FILES

