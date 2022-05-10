import pickle
import numpy as np
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
from sepia.SepiaPredict import SepiaFullPrediction

with open('data/mvcDataSimStudy_50.pkl','rb') as f:
    mvcData = pickle.load(f)
with open('data/tStar.pkl','rb') as f:
    tStar = pickle.load(f)
with open('data/Xpred_es.pkl','rb') as f:
    Xpred = pickle.load(f)
    
nSamp = 2500
nBurn = 1000

N = 1000
nPredSamp = 100
n = Xpred.shape[0]
nY = mvcData[1]['Ydata']['obs']['orig'].shape[0]
predy = np.zeros((N,nPredSamp,n,nY))
for i in range(N):
    if i%100==0:
        print(i/10,'%')
    data = SepiaData(x_sim = mvcData[i]['XTdata']['sim']['X']['orig'],
                 t_sim = mvcData[i]['XTdata']['sim']['T']['orig'],
                 y_sim = mvcData[i]['Ydata']['sim']['orig'].T,
                 y_ind_sim = mvcData[i]['YindSim'].squeeze(),
                 x_obs = mvcData[i]['XTdata']['obs']['X']['orig'],
                 y_obs = mvcData[i]['Ydata']['obs']['orig'].T,
                 y_ind_obs= mvcData[i]['YindObs'].squeeze())
    data.transform_xt()
    data.standardize_y(scale='scalar')
    data.create_K_basis(n_pc=.99)
    model = SepiaModel(data)
    model.tune_step_sizes(50, 20,prog=False,verbose=False)
    model.do_mcmc(nSamp,prog=False)
    pred_samples = model.get_samples(nPredSamp,nburn=nBurn)
    pred=SepiaFullPrediction(x_pred=Xpred, samples=pred_samples, model=model)
    predy[i,:,:,:] = pred.get_yobs(as_obs=True)
    
    
with open('data/sepia_pred/sepia_ypred_es_m50.pkl','wb') as f:
    pickle.dump(predy,f)
