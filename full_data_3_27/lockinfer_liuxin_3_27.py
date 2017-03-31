#-*- encoding:utf-8
import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import svds
import math
import random
import pickle

dir_name = '/home/ubuntu/liuxin/anti-fraud-experiment/full_data_3_27/'
data_name = "/home/ubuntu/data/twitter_rv.net"
singular_k = 5
split_chrac = "\t"
#stop_n = 90000000#used for data reading

def select_seeds(u,ii,jj):
    krou,ktheta = 5,5
    vectorij = []
    for arr in u:
        vectorij.append([float(arr[ii-1]),float(arr[jj-1])])
    N = len(vectorij)
    routhetaij = []
    minrou,maxrou = 1.0,-1.0
    mintheta,maxtheta = -90.0,90.0
    for [entryi,entryj] in vectorij:
        rou = (entryi**2+entryj**2)**0.5
        if entryi == 0:
            if entryj > 0: theta = 90.0
            if entryj < 0: theta = -90.0
            if entryj == 0: theta = 0.0
        else:
            theta = math.atan(entryj/entryi)/math.pi*180.0
        routhetaij.append([rou,theta])
        minrou,maxrou = min(minrou,rou),max(maxrou,rou)
    gaprou,gaptheta = (maxrou-minrou)/krou,(maxtheta-mintheta)/ktheta
    freqrou,freqtheta = [0 for i in range(0,krou)],[0 for i in range(0,ktheta)]
    freq = [[0 for j in range(0,ktheta)] for i in range(0,krou)]
    u2irou,u2itheta = {},{}
    for u in range(0,N):
        [rou,theta] = routhetaij[u]
        irou = max(min(int((rou-minrou)/gaprou),krou-1),0)
        itheta = max(min(int((theta-mintheta)/gaptheta),ktheta-1),0)
        freqrou[irou] += 1
        u2irou[u] = irou
        if irou > 0:
            freqtheta[itheta] += 1
            u2itheta[u] = itheta
        freq[irou][itheta] += 1
    posrou = 0
    seed_start = False
    for irou in range(0,krou):
        if seed_start and freqrou[irou] > 0:
            posrou = irou
            break
        if freqrou[irou] == 0: seed_start = True
    seedrou,seedtheta = set(),set()
    for u in u2irou:
        if u2irou[u] >= posrou:
            seedrou.add(u)
    for u in u2itheta:
        seedtheta.add(u)
    return seedrou & seedtheta

def blame_fan_idol(nodeij,seeds,N):#edge list,seed list, node number
    D = 1.0*len(nodeij)/N/N
    n = 10
    density = (1.0/n*math.log(1.0*n/N)*2)/math.log(D)
    injectedfans,injectedidols = set(),set()
    blamedfans,blamedidols = set(seeds),set()
    lastnum_fan,lastnum_idol = 0,0
    iter_num = 0
    while not (len(blamedfans) == lastnum_fan and len(blamedidols) == lastnum_idol):
        iter_num += 1
        lastnum_fan = len(blamedfans)
        lastnum_idol = len(blamedidols)
        nodej2num = {}
        for [nodei,nodej] in nodeij:
            if nodei in blamedfans:
                if not nodej in nodej2num:
                    nodej2num[nodej] = 0
                nodej2num[nodej] += 1
        blamedidols = set()
        numedge = 0
        sort_nodej2num = sorted(nodej2num.items(),key=lambda x:-x[1])
        for item in sort_nodej2num:
            nodej,num = item[0],item[1]
            blamedidols.add(nodej)
            numedge += num
            if 1.0*numedge < density*len(blamedfans)*len(blamedidols):
                break
        nodei2num = {}
        for [nodei,nodej] in nodeij:
            if nodej in blamedidols:
                if not nodei in nodei2num:
                    nodei2num[nodei] = 0
                nodei2num[nodei] += 1
        blamedfans = set()
        numedge = 0
        sort_nodei2num = sorted(nodei2num.items(),key=lambda x:-x[1])
        for item in sort_nodei2num:
            nodei,num = item[0],item[1]
            blamedfans.add(nodei)
            numedge += num
            if 1.0*numedge < density*len(blamedfans)*len(blamedidols):
                break
        print "Blame Fans&Idols Iteration: ",iter_num,"; Blamed fans: ",len(blamedfans),"; Blamed idols: ", len(blamedidols)
    print "迭代了"+str(iter_num)+"轮"
#     print "precisionfan, recallfan, precisionidol, recallidol"
#     print precisionfan,recallfan,precisionidol,recallidol
    return blamedfans,blamedidols


edge_i = []
edge_j = []
N = 0# number of nodes
iter_n = 0

f = open(data_name,"r")
for line in f.xreadlines():
    tmp = line.strip().split(split_chrac)
    edge_i.append(int(tmp[0]))
    edge_j.append(int(tmp[1]))
    N = max(N,int(tmp[0]))
    N = max(N,int(tmp[1]))
    iter_n += 1
    if (iter_n%10000000==0):
        print "Reading Data Iteration:",iter_n/1000000
#     if (iter_n==stop_n):
#         break
print 
print "Reading data done.\n\n"

#Sparse Matrix decomopistion
data = [1.0] * len(edge_i)
sparse_m = scipy.sparse.coo_matrix((data, (edge_i, edge_j)), shape=(N+1,N+1))
print "Sparse matrix done.\n\n"

with open(dir_name+'sparse_m', 'wb') as handle:
    pickle.dump(sparse_m, handle)

svd_result = svds(sparse_m, k=singular_k)
u,s,v = svd_result[0],svd_result[1],svd_result[2]
print "SVD done."

with open(dir_name+'u.pickle', 'wb') as handle:
    pickle.dump(u, handle)

with open(dir_name+'v.pickle', 'wb') as handle:
    pickle.dump(v, handle)
    
seeds = set()
for i in range(1,singular_k+1):
    for j in range(1,1+singular_k):
        if i==j:
            continue
        cur_seeds = select_seeds(u,i,j)
        seeds = seeds.union(cur_seeds)
    print "Seed Selection Iteration:",i
print "Seed selection done.\n\n"


print "seeds size:",len(seeds)
blamed_fans,blamed_idols = blame_fan_idol(nodeij,seeds,N)
print "Blamed_fans:",len(blamed_fans),"Blamed_idols:",len(blamed_idols)



f_name = "bm_fans.txt"
f = open(f_name,"w")
for i in blamed_fans:
    f.write(str(i)+"\n")
f.close()

f_name = "bm_idols.txt"
f = open(f_name,"w")
for i in blamed_idols:
    f.write(str(i)+"\n")
f.close()
print "Finally Done.\n\n"