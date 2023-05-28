

# Python program to illustrate the intersection
# of two lists in most simple way
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    not_lst3 = [value for value in lst1 if value not in lst2]
    return lst3,not_lst3

for state in [0,1]:
    if state==0:state_name='Quiet'
    else:state_name='Movement'
    for cluster in [0,1]:
        print(state_name,'  \tCluster:',cluster,'\t', len(store_cluster_dictionary[state][cluster]))
    
print('\nIntersection:')
for cluster in [0,1]:
    # Compare clusters
    inters_list,remaining = intersection(store_cluster_dictionary[0][cluster], store_cluster_dictionary[1][cluster])
    inters_list_reverse,remaining_reverse = intersection(store_cluster_dictionary[1][cluster], store_cluster_dictionary[0][cluster])

    print('Cluster:',cluster,'\t',len(inters_list), '\tleftover:', len(remaining))
    print('Cluster:',cluster,'\t',len(inters_list_reverse), '\tleftover:', len(remaining_reverse),'\n\n')
