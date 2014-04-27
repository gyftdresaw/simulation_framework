
# routines for saving and loading simulations 
import cPickle
import numpy as np

# making arguments optional so ordering doesn't matter
def save(filename,t=None,cdata=None,centers=None,sim=None):
    to_file = open(filename,'wb') # writing in binary mode
    # so we don't need to worry about ordering, save as dictionary
    dump_data = {'t':t,'cdata':cdata,'centers':centers,'sim':sim}
    cPickle.dump(dump_data,to_file)
    to_file.close()

def load(filename):
    from_file = open(filename,'rb')
    dump_data = cPickle.load(from_file)
    # for some reason the numpy arrays are segfaulting unless they get reinitialized
    dump_data['t'] = np.array(dump_data['t'])
    dump_data['cdata'] = np.array(dump_data['cdata'])
    dump_data['centers'] = np.array(dump_data['centers'])
    return dump_data # dictionary with 4 keys defined above in save
