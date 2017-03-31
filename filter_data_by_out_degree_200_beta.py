# -*- encoding:utf-8 -*-
f_raw_name = "/home/ubuntu/data/twitter_rv.net"
f_cleaned_name = "/home/ubuntu/data/twitter_rv_pro_200.net"
f_tmp_name = "/home/ubuntu/data/twitter_rv_tmp_200.net"
f_dict_name = "/home/ubuntu/data/twitter_rv_dict_200.net"

iter_n = 0
stop_n = 90000000
min_thre = 200#出度下限

f_raw = open(f_raw_name,"r")
f_cleaned = open(f_cleaned_name,"w")

#删除和左边有关的边
last_ID = 0
last_edges = []
remove_ids = []
for line in f_raw.xreadlines():
    cur_edge = line.strip("\n").split("	")
    left = cur_edge[0]
    right = cur_edge[1]
    if left==last_ID:
        last_edges.append(cur_edge)
    else:
        if len(last_edges)>min_thre:
            for edge in last_edges:
                f_cleaned.write(str(edge[0]+"	"+edge[1]+"\n"))
        else:
            remove_ids.append(last_ID)
        last_ID = left
        last_edges = [cur_edge]
    iter_n += 1
    if (iter_n%10000000==0):
        print iter_n/1000000
#     if (iter_n==stop_n):
#         break
if len(last_edges)>min_thre:
    for edge in last_edges:
        f_cleaned.write(str(edge[0]+"	"+edge[1]+"\n"))
else:
    remove_ids.append(last_ID)
f_raw.close()
f_cleaned.close()

#删除和右边有关的边
f_cleaned = open(f_cleaned_name,"r")
f_tmp = open(f_tmp_name,"w")
for line in f_cleaned.xreadlines():
    cur_edge = line.strip("\n").split("	")
    left = cur_edge[0]
    right = cur_edge[1]
    if right not in remove_ids:
        f_tmp.write(str(edge[0]+"	"+edge[1]+"\n"))  
print "Done"

#hash
f_cleaned = open(f_cleaned_name,"w")
f_tmp = open(f_tmp_name,"r")
d = {}
num = 0
for line in f_tmp.xreadlines():
    cur_edge = line.strip("\n").split("	")
    left = cur_edge[0]
    right = cur_edge[1]
    if left not in d.keys():
        d[left]=num
        num = num +1
    if right not in d.keys():
        d[right] = num
        num = num + 1
    f_tmp.write(str(d[left]+"	"+d[right]+"\n"))  
    
# !!!!!!!!!!!!!!
f_dict = open(f_dict_name,"w")
f_dict.write("{")