source parameters.sh

num_sad_faces=0

# samtools must be installed - samtools_path must be to executable
if [ -x $samtools_path ]; then
  echo 🙌 Samtools is correctly installed
else
  echo 😢 Check samtools installation
  num_sad_faces=$((num_sad_faces+1))
fi

# python 2 must be installed
# python path in the parameters file must be python executable
if [ -x $python_path ]; then
  ret=`$python_path -c 'import sys; print("%i" % (sys.hexversion<0x03000000))'`
  if [ $ret -eq 0 ]; then
      echo 😢 Python 2 must be installed
      num_sad_faces=$((num_sad_faces+1))
  else
      echo 🙌 Python 2 is correctly installed
  fi
else
  echo 😢 Check Python installation
  num_sad_faces=$((num_sad_faces+1))
fi

# STAR mapper must be installed with accompanying genomes
if [ -x $STAR ]; then
  echo 🙌 STAR is correctly installed
else
  echo 😢 Check STAR installation
  num_sad_faces=$((num_sad_faces+1))
fi

# check if the STAR index exists
if [ -f $index ]; then
  echo 🙌 $index is found
else
  echo 😢 Please check STAR index
  num_sad_faces=$((num_sad_faces+1))
fi

# Check TrimGalore installation
if [ -x $trim_galore_path ]; then
  echo 🙌 TrimGalore is correctly installed
else
  echo 😢 Check TrimGalore installation
  num_sad_faces=$((num_sad_faces+1))
fi

# report the number of failed package installations
if ! [ $num_sad_faces -eq 0 ]; then
  echo 🆘 $num_sad_faces packages not installed correctly. Please check again.
else
  echo ★ All packages installed correctly! You can now run sciRNAseq
fi
