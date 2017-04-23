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
    

        ```GEM('./X.txt','./loci.txt',5E12,4,1E4,0)```

* Parameters

    First: File name of Hi-C map. 

    Second: File name of genomic loci.

    Third: Energy coefficient. Default is 5E12.
    
    Forth: Number of conformations. Default is 4.
    
    Fifth: Maximum number of iterations. Default is 1E4.
    
    Sixth: Whether to infer the latent function (1/0).


* Input file

    File of Hi-C map: A N × N (N is the number of genomic loci) symmetric matrix separated by the table delimiter. The elements of this matrix represent interaction frequencies of Hi-C map. Example: X.txt.
    
    File of genomic loci: A 𝑁 × 1 matrix separated by the table delimiter. The elements of this matrix represent the sequence position of the genomic loci. Example: loci.txt.

* Output information

    The final total cost, data cost, energy cost and inferred latent funtion (optional) are shown.

* Output file

    structure.txt: The reconstructed chromatin structure, i.e., an ensemble of conformations (N × 3 × M matrix). M is the number of conformations.

* Parameter selection

    The energy coefficient depends on how much the users concentrate on energy stability. It is a trade-off between spatial constraint from Hi-C data and energy restriction. Users can set the parameter according to their emphasized aspect. Additionally, there are alternative ways to select the parameter automatically, such as Bayesian approach and TOPSIS. We provide the implement of Bayesian approach here. If you desire better parameters, implement Bayesian parameter selection by inputting the following at the MATLAB command line:
    
    ```BayesParaSelect(beginpara,endpara,real_volume)```
    
    beginpara & endpara: Select parameters in the range of [beginpara, endpara]. Default is [5E8, 5E16].
    
    real_volume: Real volume of the chromatin. If you do not have the priori information of the real volume, you can set real_volume -1 to use the estimated value provided by GEM. 
    
    Note that, the parameter selection is time-consuming. Fortunately, the default setting is enough in general, which was argued in the paper of GEM.

## NOTES
This software was developed and tested on MATLAB R2010b/R2014b/2016a and Windows/Linux operating systems.


## CONTACTS
Comments and bug-reports are higly appreciated. 

Guangxiang Zhu, Tsinghua University

insmileworld@gmail.com
