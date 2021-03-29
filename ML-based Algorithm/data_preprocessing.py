import numpy as np
import pandas as pd
import sys
from sklearn.preprocessing import OneHotEncoder
def data_process(path):
    dataset=pd.read_csv(path)
    x=dataset.iloc[1:,:-1].values
    x=pd.DataFrame(x)
    num_samp=np.size(x,0)
    num_col=np.size(x,1)
    y = dataset.iloc[1:, num_col].values
    #model_encode=x[:,0].reshape(num_samp,1)
    #x_ohe=OneHotEncoder().fit(model_encode)
    #model_encoded = x_ohe.transform(model_encode).toarray()
    #x_encoded = np.hstack((model_encoded,x[:,1:]))
    y=y.reshape(num_samp,1)
    y_ohe=OneHotEncoder().fit(y)
    y_encoded = y_ohe.transform(y).toarray()
    y_encoded = pd.DataFrame(y_encoded)
    
    return x,y_encoded

if __name__ == '__main__':
    path = sys.argv[1]
    X,Y=data_process(path)
    X.to_csv('sample.csv')
    Y.to_csv('target.csv')
    np.savetxt('sample.csv', X, delimiter=',')
    np.savetxt('target.csv', Y, delimiter=',')
    print("complete")
