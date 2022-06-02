import pickle
import numpy as np
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
import sys

# pass arg
min = int(sys.argv[1])
max = int(sys.argv[2])
n = max-min
m = int(sys.argv[3])

load_file = 'data/mvcDataSimStudy_'+str(m)+'.pkl'
print('loading data from',load_file,'\n')
with open(load_file,'rb') as f:
    mvcData = pickle.load(f)
with open('data/tStar.pkl','rb') as f:
    tStar = pickle.load(f)
save_file = 'data/bd_ub_sepia_sim_results_m'+str(m)+'.pkl'
print('results will be saved to',save_file,'\n')

tSamp = np.zeros((1500,n))

for i in range(min,max):
    if i%(n/10)==0: print(i,'\n')
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
    tSamp[:,(i-min)] = model.get_samples(nburn=1000)['theta'].squeeze()

print('saving...\n')

with open('data/sepia_'+str(m)+'_'+str(min)+'-'+str(max)+'.pkl','wb') as f:
    pickle.dump(tSamp,f)
#with open(save_file,'wb') as f:
#    pickle.dump(tSamp,f)
print('done.\n')
