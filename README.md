# GEM
GEM (**G**enomic organization reconstructor based on conformational **E**nergy and **M**anifold learning) is a manifold learning framework for reconstructing spatial organizations of chromosomes.

## CITATION
Guangxiang Zhu, Wenxuan Deng, Hailin Hu, Rui Ma, Sai Zhang, Tommy Kaplan, Jian Peng, and Jianyang Zeng*. ``A manifold learning based framework for reconstructing spatial organizations of chromosomes'', *Under review*, 2017.


## INSTALLATION
Put all the MATLAB script files in your MATLAB path. 

## USAGE
Directly run the m-file **GEM.m** with parameters in the directory. 

* Example

    do the following at the MATLAB command line: 
    

		 GEM('./HiC.txt', './loci.txt', 1E4, 4, 5E12, 0)

* Parameters

		GEM(HiC_file, loci_file, max_iter, M, lambdaE, infer_latent)

    HiC_file: File name of Hi-C map. 

    loci_file: File name of genomic loci.

    max_iter: Maximum number of iterations. Default is 1E4.
    
    M: Number of conformations. Default is 4.
    
    lambdaE: Coefficient of energy term. Default is 5E12.
    
    infer_latent: Whether to infer the latent function (1/0).


* Input file

    File of Hi-C map: A N × N (N is the number of genomic loci) symmetric matrix separated by the table delimiter. The elements of it represent interaction frequencies of the Hi-C map. Example: HiC.txt.
    
    File of genomic loci: A N × 1 matrix separated by the table delimiter. The elements of it represent the sequence position of the genomic loci. Example: loci.txt.
    
	The example files (HiC.txt and loci.txt) contain normalized the Hi-C map and the genomic loci of chromosome 14 for 1Mb bins, which are derived from Yaffe et al. (http://compgenomics.weizmann.ac.il/tanay/?page_id=283).

* Output information

	Final total cost, data cost, energy cost and inferred latent funtion (optional).

* Output file

    conformation[1-M].txt: The reconstructed chromatin conformation 1-M. Each conformation is a N × 3 matrix. 
    
    proportions.txt: The corresponding weights of conformations (M × 1 matrix).

* Parameter selection

    The coefficient of energy term depends on how much the users concentrate on energy stability. It is a trade-off between spatial constraint from Hi-C data and energy restriction. Users can set the parameter according to their emphasized aspect. Additionally, there are alternative ways to select the parameter automatically, such as Bayesian approach and TOPSIS. We provide the implement of Bayesian approach here. If you desire better parameters, implement Bayesian parameter selection by inputting the following at the MATLAB command line:
    
	 	 BayesParaSelect(begin_para, end_para, real_volume, HiC_file, loci_file, max_iter, M, infer_latent)
    

	begin_para and end_para: Select parameters in the range of [5×10^begin_para, 5×10^end_para]. Default are 8 and 16. For example, if you set begin_para 8 and end_para 16, GEM will select the best parameter from [5E8,5E9,5E10,5E11,5E12,5E13,5E14,5E15,5E16].
	
	real_volume: Real volume of the chromatin. If you do not have the priori information of the real volume, you can set real_volume -1 or -2 to use the estimated value provided by GEM (-1 for human cell, -2 for yeast cell).

	HiC_file: File name of Hi-C map. 
	
	loci_file: File name of genomic loci.
	
	max_iter: Maximum number of iterations. Default is 1E4.
    	
	M: Number of conformations. Default is 4.
	
	infer_latent: Whether to infer the latent function (1/0).
    
    Considering that the parameter selection is time-consuming, intact automatical parameter selection is not always necessary. Fortunately, the default setting is good enough in general, which is argued in the paper of GEM. Also, users can fine-tune the default setting of parameters according to the output data cost and energy cost.

## NOTES
This software was developed and tested on MATLAB R2010b/R2014b/2016a and Windows/Linux operating systems.


## CONTACTS
Comments and bug-reports are higly appreciated. 

Guangxiang Zhu, Tsinghua University

insmileworld@gmail.com
