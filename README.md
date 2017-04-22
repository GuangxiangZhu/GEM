# GEM
GEM (**G**enomic organization reconstructor based on conformational **E**nergy and **M**anifold learning) is a manifold learning framework for reconstructing spatial organizations of chromosomes.

## CITATION
Guangxiang Zhu, Wenxuan Deng, Hailin Hu, Rui Ma, Sai Zhang, Tommy Kaplan, Jian Peng, and Jianyang Zeng*. ``A manifold learning based framework for reconstructing spatial organizations of chromosomes'', *Under review*, 2017.


## INSTALLATION
Put all the MATLAB script files in your MATLAB path. 

## USAGE
Directly run the m-file **GEM.m** with parameters in the directory. 

* Example

    do the following at the MATLAB command line:  $$a_1$$

    ```GEM(5E12,4,1E4)```

* Parameters

    First： Energy coefficient. Default is 5E12.

    Second： Number of conformations. Default is 4.

    Third： Maximum number of iterations. Default is 1E4.



* Input data

* Output data

* Parameter selection

    The energy coefficient depends on how much the users concentrate on energy stability. It is a trade-off between spatial constraint from Hi-C data and energy restriction. Users can set the parameter according to their emphasized aspect. Additionally, there are alternative ways to select the parameter automatically, such as Bayesian approach and TOPSIS. We provide the implement of Bayesian approach here. If you desire better parameters, implement Bayesian parameter selection by inputting the following at the MATLAB command line:
    
    ```BayesParaSelect```
    
    Note that, the parameter selection is time-consuming. Fortunately, the default setting is enough in general, which was argued in the paper of GEM.

## NOTES
This software was developed and tested on MATLAB R2010b/R2014b/2016a and Windows/Linux operating systems.


## CONTACTS
Comments and bug-reports are higly appreciated. 

Guangxiang Zhu, Tsinghua University

insmileworld@gmail.com
