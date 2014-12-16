
# Part of software for the ECML'14 submission "Reliability maps: A tool to enhance probability estimates and improve classification accuracy"
# Author: Meelis Kull, meelis.kull@bristol.ac.uk, meelis.kull@gmail.com
# Written February-April 2014.
# All rights reserved.



sudo Rscript code_clean.R "diabetes" 0.01 0.10 200 1 &
sudo Rscript code_clean.R "australian" 0.01 0.10 200 1 &
sudo Rscript code_clean.R "breast-cancer" 0.01 0.10 200 1 &
sudo Rscript code_clean.R "german.numer" 0.01 0.10 200 1 &
sudo Rscript code_clean.R "mushrooms" 0.01 0.10 200 1 &
sudo Rscript code_clean.R "svmguide1" 0.01 0.10 200 1 &
sudo Rscript code_clean.R "diabetes" 0.01 0.10 200 2 &
sudo Rscript code_clean.R "australian" 0.01 0.10 200 2 &
sudo Rscript code_clean.R "breast-cancer" 0.01 0.10 200 2 &
sudo Rscript code_clean.R "german.numer" 0.01 0.10 200 2 &
sudo Rscript code_clean.R "mushrooms" 0.01 0.10 200 2 &
sudo Rscript code_clean.R "svmguide1" 0.01 0.10 200 2 &
sudo Rscript code_clean.R "diabetes" 0.01 0.10 200 3 &
sudo Rscript code_clean.R "australian" 0.01 0.10 200 3 &
sudo Rscript code_clean.R "breast-cancer" 0.01 0.10 200 3 &
sudo Rscript code_clean.R "german.numer" 0.01 0.10 200 3 &
sudo Rscript code_clean.R "mushrooms" 0.01 0.10 200 3 &
sudo Rscript code_clean.R "svmguide1" 0.01 0.10 200 3 &
sudo Rscript code_clean.R "diabetes" 0.02 0.20 200 1 &
sudo Rscript code_clean.R "australian" 0.02 0.20 200 1 &
sudo Rscript code_clean.R "breast-cancer" 0.02 0.20 200 1 &
sudo Rscript code_clean.R "german.numer" 0.02 0.20 200 1 &
sudo Rscript code_clean.R "mushrooms" 0.02 0.20 200 1 &
sudo Rscript code_clean.R "svmguide1" 0.02 0.20 200 1 &


sudo Rscript code_clean.R "cod-rna" 0.01 0.10 200 2 &
sudo Rscript code_clean.R "splice" 0.01 0.10 200 2 &
sudo Rscript code_clean.R "liver-disorders" 0.01 0.10 200 2 &

sudo Rscript code_clean.R "liver-disorders" 0.05 0.20 200 2
sudo Rscript code_clean_parallel.R "liver-disorders" 0.05 0.20 200 2

sudo Rscript code_clean_parallel.R "splice" 0.05 0.20 200 2
sudo Rscript code_clean_parallel.R "liver-disorders" 0.05 0.20 200 2
sudo Rscript code_clean_parallel.R "svmguide1" 0.05 0.20 200 2
sudo Rscript code_clean_parallel.R "australian" 0.05 0.20 200 2
sudo Rscript code_clean_parallel.R "breast-cancer" 0.05 0.20 200 2
sudo Rscript code_clean_parallel.R "diabetes" 0.05 0.20 200 2