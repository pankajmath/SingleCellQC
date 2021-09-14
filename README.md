# SingleCellQC
When we have a number of biological experiments conducted on different dates with same controls in each experiment and we wish to know which controls in these experiments might be less than ideal, this code helps flag such experiemnts/controls. In summary, we follow these steps:
 
1. obtain the distribution of single cell expression for the vehicle condition
2. normalize it using median/MAD
3. get the quantile values (1st to 99th)
4. align the quantile curve to the standard normal distribution using the dynamic time warping (DTW).

Next, we compute the pairwise EMD distance between controls from different experiemnts and cluster them to get the experiements where controls group together as well as any outliers.


## Reference
If you use this code for your research/work, please cite the following article:

Fabio Stossi, Pankaj K. Singh, Ragini M. Mistry, Hannah L. Johnson, Radhika D. Dandekar, Maureen G. Mancini, Arvind U. Rao, Michael A. Mancini
_Quality Control for Single Cell Imaging Analytics Using Endocrine Disruptor-induced Changes in Estrogen Receptor Expression_ 


