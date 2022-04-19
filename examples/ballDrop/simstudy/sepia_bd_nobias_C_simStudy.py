import pickle
import numpy as np
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData

with open('mvcDataSimStudy.pkl','rb') as f:
    mvcData = pickle.load(f)

for i in range(len(mvcData)):
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
    model.do_mcmc(2500,prog=False)
    mvcData[i]['tSamp'] = model.get_samples(nburn=1000)['theta']
    
with open('mvcDataSimStudy.pkl','wb') as f:
    pickle.dump(mvcData,f)
