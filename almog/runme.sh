for i in 0 47 65 77 83 89 93 102 105
do
  date
  echo "python calcifcurves.py $i"
  python calcifcurves.py $i
  date
  echo "python findDCshortthreshold.py $i"
  python findDCshortthreshold.py $i
  date
  echo "python calcsteadystate.py $i"
  python calcsteadystate.py $i
  date
done
