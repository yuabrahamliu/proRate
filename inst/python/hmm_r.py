# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 19:40:41 2023

@author: yuabr
"""

#%% Interface

def scmerfunction(fname2load, trainsub, anno2load, y2load, K, k):
    
    #Config
    import rpy2.robjects as robjects
    
    import scanpy as sc
    
    import pandas as pd
    import numpy as np
    
    import scmer
    
    #Interface between R and Python
    betasfile = fname2load
    
    robjects.r['load'](betasfile)
    betas = robjects.r[trainsub]
    
    
    
    #Load anno dataframe from R RData file as robjects.vectors.DataFrame
    robjects.r['load'](anno2load)
    totalanno = robjects.r['anno']
    
    tst = pd.DataFrame(totalanno)
    tst = tst.T
    tst.columns = totalanno.colnames
    anno = tst.set_index('sentrix', drop=False)
    anno = anno.loc[betas.rownames,]
        
    #Load y factor from R RData file as robjects.vectors.FactorVector
    robjects.r['load'](y2load)
    totaly = robjects.r['y']
    
    totalyvals = list(totaly)
    totalylabels = [totaly.levels[i-1] for i in totalyvals]
    totallabels = pd.DataFrame({'label': totalylabels, 'val': totalyvals})
    totallabels.index = tst.sentrix
    
    y = totallabels.loc[betas.rownames,]
    
    del tst, totalyvals, totalylabels, totallabels

    y.index = betas.rownames
    
    #Transfer the rojbects.vectors.FloatMatrix to Scanpy as AnnData
    tst = np.asarray(betas)
    adata = sc.AnnData(tst)
    
    del tst
    
    adata.obs = y
    adata.obs_names = betas.rownames
    adata.var_names = betas.colnames
    
    #Scanpy analysis
    sc.settings.set_figure_params(dpi = 50, facecolor = 'white')
    sc.pl.highest_expr_genes(adata, n_top = 20)
    
    sc.tl.pca(adata, svd_solver = 'arpack')
    sc.pl.pca(adata, color = 'label', title = 'Before', 
              save = 'before.' + K + '.' + k)
    sc.pl.pca_variance_ratio(adata, log = False)
        
    sc.pp.neighbors(adata, n_neighbors = 30, n_pcs = 30)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color = 'label', title = 'Before', 
               save = 'before.' + K + '.' + k)
    
    #SCMER feature selection
    model = scmer.UmapL1(lasso = 5.5e-6, ridge = 0, 
                         n_pcs = 100, perplexity = 30., 
                         use_beta_in_Q = False, n_threads = 6, pca_seed = 1234)
    model.fit(adata.X)
    
    print(*adata.var_names[model.get_mask()])
        
    #PCA and UMAP validation
    new_adata = model.transform(adata)
    
    sc.tl.pca(new_adata, svd_solver = 'arpack')
    sc.pl.pca(new_adata, color = 'label', title = 'After', 
              save = 'after.' + K + '.' + k)
    
    sc.pp.neighbors(new_adata, n_pcs = 30, use_rep = "X_pca", n_neighbors = 30)
    sc.tl.umap(new_adata)
    sc.pl.umap(new_adata, color = 'label', title = 'After', 
               save = 'after.' + K + '.' + k)
    
    features = list(adata.var_names[model.get_mask()])
    
    savename = 'features.' + K + '.' + k
    robjects.r.assign(savename, features)
    robjects.r("save(" + savename + ", file = '" + savename + ".RData')")
    
        
    return features

#%% Config

import os
import pandas as pd
import numpy as np

import rpy2.robjects as robjects

from hmmlearn import hmm

#wkdir = 'C:\\Users\\yuabr\\Desktop\\Transfer\\codetransfer\\proRate\\proRate_V2\\'

wkdir = 'C:\\Users\\Yu Liu\\Desktop\\proRate\\data\\'

os.chdir(wkdir)

ratiosfile = '.\\ratios.RData'

robjects.r['load'](ratiosfile)

ratios = robjects.r['ratios']

ratios = pd.DataFrame(ratios)

ratios = ratios.T

ratios.columns = ['ratio_adj']

def gaussion_prior(x):
    
    x = x[:-1]
    
    mean_prior = np.apply_along_axis(np.mean, 1, np.array([x[0:3], x[-3:]]))
    
    var_prior = np.apply_along_axis(np.var, 1, np.array([x[0:3], x[-3:]]))
    
    if mean_prior[0] >= mean_prior[1]: 
        
        x = x[np.logical_and(x > np.quantile(x, 0.1), x < np.quantile(x, 0.9))]
        #logical_and(x1, x2)
        #Compute the truth value of x1 AND x2 element-wise.
        
        x = np.sort(x)
        
        mean_prior = np.apply_along_axis(np.mean, 1, np.array([x[0:3], x[-3:]]))
        
        var_prior = np.apply_along_axis(np.var, 1, np.array([x[0:3], x[-3:]]))
        
    return mean_prior.reshape(2, -1), var_prior.reshape(2, -1)

startprob_prior = np.array([1, 0])

expand_ratio = ratios.iloc[:,0].values

expand_ratio = expand_ratio.reshape(-1, 1)

transmat_prior = np.array([[0.9, 0.1], 
                           [0, 1]])

means_prior, vars_prior = gaussion_prior(x = expand_ratio)



model = hmm.GaussianHMM(n_components=2, 
                        n_iter=100, 
                        random_state=1234, 
                        
                        init_params='tmcs', 
                        params='tmc', 
                        
                        transmat_prior = transmat_prior, 
                        
                        means_prior = means_prior, 
                        covars_prior = vars_prior, 
                        
                        startprob_prior = startprob_prior)

#model.startprob_ = startprob_prior

model.fit(expand_ratio[:-1])

for i in range(99):
    
    if np.sum(np.abs(model.transmat_[1:,] - transmat_prior[1:,])) > 10**-3: 
        
        model.transmat_[1:,] = transmat_prior[1:,]
        
        transmat_prior = model.transmat_
        
        model = hmm.GaussianHMM(n_components=2, 
                                n_iter=100, 
                                random_state=1234, 
                                
                                init_params='tmcs', 
                                params='tmc', 
                                
                                transmat_prior = transmat_prior, 
                                
                                means_prior = means_prior, 
                                covars_prior = vars_prior, 
                                
                                startprob_prior = startprob_prior)

        #model.startprob_ = startprob_prior

        model.fit(expand_ratio[:-1])
    
    else: 
        
        break
    
    print(i)
    

model.transmat_[1:,] = transmat_prior[1:,]


try: 
    
    pres = model.predict(expand_ratio)
    
    if np.unique(pres).shape[0] == 1: 
        
        res = None
    
    else: 
        
        final_point = np.where(pres == 0)[0][-1]
    
        fronts = expand_ratio[:final_point + 1,0]
        latters = expand_ratio[final_point + 1:,0]
        
        final_point_r = final_point + 1
        
        if len(latters) == 0 or len(fronts) == 0: 
            
            res = None
            
        else: 
            
            res = {'point': final_point_r, 
                   'fronts': fronts, 
                   'latters': latters}
        
except ValueError: 
    
    res = None




def hmm_r(ratios, 
          hmmseed=2023):
    
    #import pandas as pd
    import numpy as np

    from hmmlearn import hmm
    
    hmmseed = int(hmmseed)
    
    def gaussion_prior(x):
        
        x = x[:-1]
        
        mean_prior = np.apply_along_axis(np.mean, 1, np.array([x[0:3], x[-3:]]))
        
        var_prior = np.apply_along_axis(np.var, 1, np.array([x[0:3], x[-3:]]))
        
        if mean_prior[0] >= mean_prior[1]: 
            
            x = x[np.logical_and(x > np.quantile(x, 0.1), x < np.quantile(x, 0.9))]
            #logical_and(x1, x2)
            #Compute the truth value of x1 AND x2 element-wise.
            
            x = np.sort(x)
            
            mean_prior = np.apply_along_axis(np.mean, 1, np.array([x[0:3], x[-3:]]))
            
            var_prior = np.apply_along_axis(np.var, 1, np.array([x[0:3], x[-3:]]))
            
        return mean_prior.reshape(2, -1), var_prior.reshape(2, -1)

    startprob_prior = np.array([1, 0])

    expand_ratio = ratios.iloc[:,0].values

    expand_ratio = expand_ratio.reshape(-1, 1)
    
    transmat_prior = np.array([[0.9, 0.1], 
                               [0, 1]])

    means_prior, vars_prior = gaussion_prior(x = expand_ratio)



    model = hmm.GaussianHMM(n_components=2, 
                            n_iter=100, 
                            random_state=hmmseed, 
                            
                            #init_params='tmcs', 
                            init_params='tmc', 
                            params='tmc', 
                            
                            transmat_prior = transmat_prior, 
                            
                            means_prior = means_prior, 
                            covars_prior = vars_prior) 
                            
                            #startprob_prior = startprob_prior)

    model.startprob_ = startprob_prior

    model.fit(expand_ratio[:-1])

    for i in range(99):
        
        if np.sum(np.abs(model.transmat_[1:,] - transmat_prior[1:,])) > 10**-3: 
            
            model.transmat_[1:,] = transmat_prior[1:,]
            
            transmat_prior = model.transmat_
            
            model = hmm.GaussianHMM(n_components=2, 
                                    n_iter=100, 
                                    random_state=hmmseed, 
                                    
                                    #init_params='tmcs', 
                                    init_params='tmc', 
                                    params='tmc', 
                                    
                                    transmat_prior = transmat_prior, 
                                    
                                    means_prior = means_prior, 
                                    covars_prior = vars_prior) 
                                    
                                    #startprob_prior = startprob_prior)

            model.startprob_ = startprob_prior

            model.fit(expand_ratio[:-1])
        
        else: 
            
            break
        
        print(i)
        

    model.transmat_[1:,] = transmat_prior[1:,]
    

    try: 
        
        pres = model.predict(expand_ratio)
        
        if np.unique(pres).shape[0] == 1: 
            
            res = None
        
        else: 
            
            final_point = np.where(pres == 0)[0][-1]
    
            fronts = expand_ratio[:final_point + 1,0]
            latters = expand_ratio[final_point + 1:,0]
            
            final_point_r = final_point + 1
            
            if len(latters) == 0 or len(fronts) == 0: 
                
                res = None
                
            else: 
                
                res = {'point': final_point_r, 
                       'fronts': fronts, 
                       'latters': latters}
            
    except ValueError: 
        
        res = None
    
    return res

    
#hmm_r(ratios=ratios)










