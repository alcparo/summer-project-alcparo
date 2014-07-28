
# Part of software for the ECML'14 submission "Reliability maps: A tool to enhance probability estimates and improve classification accuracy"
# Author: Meelis Kull, meelis.kull@bristol.ac.uk, meelis.kull@gmail.com
# Written February-April 2014.
# All rights reserved.

Rscript binary_script_real.r "data_banknote_authentication" 0.01 0.10 200
sudo Rscript binary_script_real.r "diabetes_num" 0.01 0.10 200


Rscript ecml14real.r "vehicle"         0.01 0.10 200
Rscript ecml14real.r "yeast"         0.01 0.10 200
Rscript ecml14real.r "segment"         0.01 0.10 200
Rscript ecml14real.r "page-blocks" 0.01 0.10 200
Rscript ecml14real.r "sat"         0.01 0.10 200
Rscript ecml14real.r "shuttle"     0.01 0.10 200
Rscript ecml14wang.r 1  400 0.01 0.10 200
Rscript ecml14wang.r 2  400 0.01 0.10 200
Rscript ecml14wang.r 3  400 0.01 0.10 200
Rscript ecml14wang.r 4  400 0.01 0.10 200
Rscript ecml14wang.r 5  400 0.01 0.10 200
Rscript ecml14wang.r 6  400 0.01 0.10 200
Rscript ecml14wang.r 7  400 0.01 0.10 200
Rscript ecml14wang.r 8  400 0.01 0.10 200
Rscript ecml14wang.r 9  400 0.01 0.10 200
Rscript ecml14wang.r 10 400 0.01 0.10 200
Rscript ecml14wang.r 1  2000 0.01 0.10 200
Rscript ecml14wang.r 2  2000 0.01 0.10 200
Rscript ecml14wang.r 3  2000 0.01 0.10 200
Rscript ecml14wang.r 4  2000 0.01 0.10 200
Rscript ecml14wang.r 5  2000 0.01 0.10 200
Rscript ecml14wang.r 6  2000 0.01 0.10 200
Rscript ecml14wang.r 7  2000 0.01 0.10 200
Rscript ecml14wang.r 8  2000 0.01 0.10 200
Rscript ecml14wang.r 9  2000 0.01 0.10 200
Rscript ecml14wang.r 10 2000 0.01 0.10 200
Rscript ecml14wang.r 1  10000 0.01 0.10 200
Rscript ecml14wang.r 2  10000 0.01 0.10 200
Rscript ecml14wang.r 3  10000 0.01 0.10 200
Rscript ecml14wang.r 4  10000 0.01 0.10 200
Rscript ecml14wang.r 5  10000 0.01 0.10 200
Rscript ecml14wang.r 6  10000 0.01 0.10 200
Rscript ecml14wang.r 7  10000 0.01 0.10 200
Rscript ecml14wang.r 8  10000 0.01 0.10 200
Rscript ecml14wang.r 9  10000 0.01 0.10 200
Rscript ecml14wang.r 10 10000 0.01 0.10 200
