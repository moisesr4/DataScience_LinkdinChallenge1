#premir modification ricardo
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 16:09:11 2017

@author: cbothore
"""


import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import pylab
import numpy as np
import pickle

def list_of_different_attribute_values(d):
    return set([v for values in d.values() for v in values])
    
def draw_graph(g, node_attribute=None, list_of_values_of_attributes=None):
    """
    Draw the graph g.

    Parameters
    ----------
    g : graph
       A networkx graph
    node_attribute : string 
       The name of the node attribute used to assign colors to the drawing
    list_of_values_of_attributes : list
        A list of all the potential values of node_attribute to assign one color
        per value.
    """
    #initialze Figure
    plt.figure(num=None, figsize=(20, 20), dpi=80)
    plt.axis('off')
    fig = plt.figure(1)
    
    pos = nx.spring_layout(g, iterations=100)

    
    if node_attribute and list_of_values_of_attributes: 
        # To associate colors to nodes according to an attribute, here college
        # build a color_map, one for each college
        color_map={}
        i=0.0
        for s in list_of_values_of_attributes:
            color_map[s]=i
            i+=1/len(list_of_values_of_attributes)
        color_map[None]=1 # for nodes without values for the attribute node_attribute
        
        # The values supplied to node_color should be in the same order as the nodes 
        # listed in G.nodes(). We take an arbitrary mapping of values color_map and 
        # generate the values list in the correct order
        #values = [color_map[G.node[node].get(node_attribute)] for node in G.nodes()] # for attributes encoded in the graph
        values=[]        
        for node in G.nodes():
            if node in node_attribute:
                if node_attribute[node]:
                    # we arbitrarily take the first value 
                    values.append(color_map[node_attribute[node][0]])   
            else:
                values.append(1)
               
        nx.draw_networkx_nodes(g,pos, cmap=plt.get_cmap('jet'), node_color=values)
    else:
        nx.draw_networkx_nodes(g,pos)
       
    nx.draw_networkx_edges(g,pos)
    nx.draw_networkx_labels(g,pos)

    cut = 1.00
    xmax = cut * max(xx for xx, yy in pos.values())
    ymax = cut * max(yy for xx, yy in pos.values())
    plt.xlim(0, xmax)
    plt.ylim(0, ymax)
    plt.show()
    pylab.close()
    del fig

def properties(g):
    """
    Computes simple and classic graph metrics.

    Parameters
    ----------
    g : graph
       A networkx graph
    """
    # networkx short summary of information for the graph g
    print(nx.info(g))
    
    # Draw the degree distribution. Powerlow distribution for a real (complex) network
    plt.figure(num=None)
    fig = plt.figure(1)
    degree_sequence=[d for n, d in g.degree()] # degree sequence
    print("Degree sequence %s" % degree_sequence)
    plt.hist(degree_sequence, bins='auto')  
    plt.title("powerlaw degree distribution")
    plt.ylabel("# nodes")
    plt.xlabel("degree")
    plt.show()
    pylab.close()
    del fig
 
    precomputed_eccentricity = nx.eccentricity(g) # costly step, we save time here!
    print("Graph density %f" % nx.density(g))
    print("Diameter (maximum eccentricity): %d" % nx.diameter(g,precomputed_eccentricity))
    print("Radius (minimum eccentricity): %d" % nx.radius(g,precomputed_eccentricity)) #The radius is the minimum eccentricity.
    print("Mean eccentricity (eccentricity(v) = the maximum distance from v to all other nodes): %s" % np.mean(list(precomputed_eccentricity.values())))
    print("Center is composed of %d nodes (nodes with eccentricity equal to radius)" % len(nx.center(g, precomputed_eccentricity)))
    print("Periphery is composed of %d nodes (nodes with eccentricity equal to the diameter)" % len(nx.periphery(g,precomputed_eccentricity)))
    print("Mean clustering coefficient %f" % np.mean(list(nx.clustering(g).values())))
    total_triangles=sum(nx.triangles(g).values())/3    
    print("Total number of triangles in graph: %d" % total_triangles)



# --------------------- Let's take a quick NetworkX tour ------------------------------#
# load the graph
G = nx.read_gexf("mediumLinkedin.gexf")

# load the profiles. 3 files for each type of attribute
# Some nodes in G have no attributes, or just location, or all of same
# Some nodes may have 2 or more colleges or employers, so we
# use dictionaries to store the attributes  
college={}
location={}
employer={}
# The dictionaries are loaded as dictionaries from the disk (see pickle in Python doc)
with open('mediumCollege.pickle', 'rb') as handle:
    college = pickle.load(handle)
with open('mediumLocation.pickle', 'rb') as handle:
    location = pickle.load(handle)
with open('mediumEmployer.pickle', 'rb') as handle:
    employer = pickle.load(handle)
print("\n First a fiew informations ! \n")
print("Nb of users with one or more attribute college: %d" % len(college))
print("Nb of users with one or more attribute location: %d" % len(location))
print("Nb of users with one or more attribute employer: %d" % len(employer))

# Let's have a look on the graph
#draw_graph(G, node_attribute=employer, list_of_values_of_attributes=list_of_different_attribute_values(employer)and list_of_different_attribute_values(college)and list_of_different_attribute_values(location))
# print some properties to understand the type of graph
#properties(G)

print("\n\n")

# to see the nodes with the highest degree
# We use sorted with as parameter the list of items (key + value) from our dict 
# The sort key is a lambda function which gives the indice of the data to sort.
# t is a tuple (key, value) from the dict
# Lambda says to sort with the value of the tuple: t[1] 
# to sort the tuple on the key, choose t[0]
sorted_nodes_by_degree=sorted([(n,d) for n, d in G.degree()], key=lambda t: t[1])
print("The 5 nodes with the highest degree are: %s" % str(sorted_nodes_by_degree[-5:]))
print("5 of the nodes with the smallest degree are: %s" % str(sorted_nodes_by_degree[:5]))
print("\n\n")

# To see the neighbors of a node
print("U7912 has %d neighbors" % len([n for n in G.neighbors('U7912')]))
print("U7912 has %d neighbors" % G.degree('U7912'))

if 'U7912' in college:
    print("U7912's college(s): %s" % str(college['U7912']))
print("Her/his neighbors'colleges are:")
for n in G.neighbors('U7912'):
    if n in college:
        print("--- U7912's neighbor %s went in college(s): %s" % (n,str(college[n])))
    else:
        print("--- %s has not filled his/her college profil" % str(n))
print("\n\n")

# Fast examination of all edges is achieved using adjacency iterators.
# Useful to explore the neighbors of each node
# Here, we compute the homophily according to the attribute college
similar_neighbors=0
total_number_neighbors=0 # to verify the number of edges ;-)!!!
for n in G.nodes():
    for nbr in G.neighbors(n):
        total_number_neighbors+=1
        if n in college and nbr in college:
            if len([val for val in college[n] if val in college[nbr]]) > 0:
                similar_neighbors+=1
homophily_C=similar_neighbors/total_number_neighbors




print("Homophily ('college' attribute), i.e. total proportion of neighbors sharing attributes: %f" % homophily_C)
print('Is our number of edges (=%d) similar to networkX property (=%d)? No? Normal! why?' % (total_number_neighbors, len(G.edges())))

print("\n\nNow your turn!!! ")
print("------------------------------------ NOS REPONSES ----------------------------------------------------")

# --------------------- Now your turn -------------------------------------#
# Lest's us say you are U19886 and you want to work at Google.
# Explore, implement your strategy to be helped in getting this job at Google.


###Computing Homophilies
similar_neighbors=0
total_number_neighbors=0 # to verify the number of edges ;-)!!!
for n in G.nodes():
    for nbr in G.neighbors(n):
        total_number_neighbors+=1
        if n in location and nbr in location:
            if len([val for val in location[n] if val in location[nbr]]) > 0:
                similar_neighbors+=1
homophily_L=similar_neighbors/total_number_neighbors

similar_neighbors=0
total_number_neighbors=0 
for n in G.nodes():
    for nbr in G.neighbors(n):
        total_number_neighbors+=1
        if n in employer and nbr in employer:
            if len([val for val in employer[n] if val in employer[nbr]]) > 0:
                similar_neighbors+=1
homophily_E=similar_neighbors/total_number_neighbors
print("---------- ",homophily_E,homophily_L,homophily_C)



###Creating tables with the lists of schools, jobs, and location and who are working there

ListsOfSchools={}
for j in G.nodes():
    if j in college:
        for k in college[j]:
            if k not in ListsOfSchools:
                ListsOfSchools[k]=[]
            ListsOfSchools[k].append(j)
ListsOfJobs={}
for j in G.nodes:
    if j in employer:
        for k in employer[j]:
            if k not in ListsOfJobs:
                ListsOfJobs[k]=[]
            ListsOfJobs[k].append(j)

ListOfLocations={}
for j in G.nodes:
    if j in location:
        for k in location[j]:
            if k not in ListOfLocations:
                ListOfLocations[k]=[]
            ListOfLocations[k].append(j)

#print(ListsOfJobs,ListOfLocations,ListsOfSchools)



###About me and my neighbors 
def same_attribute(a,b,list):
    for i in list:
        if a in list[i]:
            if b in list[i]:
                return True ### here we might be able to quantify : if they were together in 2 different uni, they might have been closer
    return False 




print('I AM LOCATED IN : ',location["U19886"],'----')

for i in G.neighbors('U19886'):
    print("About my friend, ",i)
    if i in college:
        
        print(college[i])
    if i in employer:
        print(employer[i])
    print(location[i])
 #   print(same_attribute('U19842',i,ListOfLocations),"that they are in the same plane \n")
    print("\n")



### weighting the graph. 

for i in G.nodes():
    for j in G.neighbors(i):
        pij=1 #probability of helping -> initialized at 1 
        if same_attribute(i,j,ListOfLocations):
            pij*=(1-homophily_L)
        if same_attribute(i,j,ListsOfJobs):
            pij*=(1-homophily_E)
        if same_attribute(i,j,ListsOfSchools):
            pij*=(1-homophily_C)
        w=pij**(1.2)
        G[i][j]['weight']=w
        #print(G[i][j]['weight'])
        G.edges

degree_sequence=[d for n, d in G.degree()]
avg=np.mean(degree_sequence)

F=nx.DiGraph() # getting the graph directed for another approximation
for i in G.nodes():
    for j in G.neighbors(i):
        F.add_node(i)
        F.add_node(j)
        pij=1 #probability of helping -> initialized at 1 
        if same_attribute(i,j,ListOfLocations):
            pij*=(1-homophily_L)
        if same_attribute(i,j,ListsOfJobs):
            pij*=(1-homophily_E)
        if same_attribute(i,j,ListsOfSchools):
            pij*=(1-homophily_C)
        w=pij**(1/2)
        degj=G.degree(j)
        degi=G.degree(i)
        F.add_edge(i,j,weight=w*(degi/avg)**(1/2))
        F.add_edge(j,i,weight=w*(degj/avg)**(1/2))
       # print(F[i][j]['weight'],F[j][i]['weight'],"\n")

# But who works at google ? 

ListOfGooglers=ListsOfJobs['google']+ListsOfJobs['google inc']

print("These are the people who work at Google ",ListOfGooglers)


def fusion(T1,T2) :
    if T1==[] :
        return T2
    if T2==[] :
        return T1
    if T1[0]<T2[0] :
        return [T1[0]]+fusion(T1[1 :],T2)
    else :
        return [T2[0]]+fusion(T1,T2[1 :])

def trifusion(T) :
    if len(T)<=1 : 
        return T
    T1=[T[x] for x in range(len(T)//2)]
    T2=[T[x] for x in range(len(T)//2,len(T))]
    return fusion(trifusion(T1),trifusion(T2))
lengths=[]
lnaive=[]
for i in ListOfGooglers:
    length_of_the_shortest_path_to_i=nx.dijkstra_path_length(G,'U19886',i)
    lengths.append(length_of_the_shortest_path_to_i)
    lnaive.append(nx.dijkstra_path_length(G,'U19886',i,'a'))

print("\n Without ponderating the graph, the lengths are : ",trifusion(lnaive))


print ("\n ---------------------------Now with the 1st approximation------------------------------")
j=lengths.index(min(lengths))#shortest path
lj=lengths[j]
lengths[j]+=5 # On augmente de manière artificielle afin de trouver les deuxièmes et troisièmes sans perturber l'indice
k=lengths.index(min(lengths))#second shortest
lk=lengths[k]
lengths[k]+=5 #idem que là haut
l=lengths.index(min(lengths))#third shortest
ll=lengths[l]
lengths[j]-=5
lengths[k]-=5


m=lengths.index(max(lengths))

print("The closest is : ",ListOfGooglers[j]," the length to ",ListOfGooglers[j]," is ", lengths[j]," while the avg length is ", np.mean(lengths)," and the longest is ", lengths[m] ,". \n The second and third lengths are ",lk,ll)
#print("\n And here are lengths corresponding to these people", lengths)
print("\n here is the shortest path: ", nx.dijkstra_path(G,'U19886',ListOfGooglers[j]))
print("\n the second and third : ", nx.dijkstra_path(G,'U19886',ListOfGooglers[k]), nx.dijkstra_path(G,'U19886',ListOfGooglers[l]))

print("\n--------------------------Now we start considering the degree of someone---------------------------------------")
lengths=[]
for i in ListOfGooglers:
    length_of_the_shortest_path_to_i=nx.dijkstra_path_length(F,'U19886',i)
    lengths.append(length_of_the_shortest_path_to_i)
j=lengths.index(min(lengths))#shortest path
lj=lengths[j]
lengths[j]+=5 # On augmente de manière artificielle afin de trouver les deuxièmes et troisièmes sans perturber l'indice
k=lengths.index(min(lengths))#second shortest
lk=lengths[k]
lengths[k]+=5 #idem que là haut
l=lengths.index(min(lengths))#third shortest
ll=lengths[l]
lengths[j]-=5
lengths[k]-=5


m=lengths.index(max(lengths))

print("The closest is : ",ListOfGooglers[j]," the length to ",ListOfGooglers[j]," is ", lengths[j]," while the avg length is ", np.mean(lengths)," and the longest is ", lengths[m] ,". \n The second and third lengths are ",lk,ll)
#print("\n And here are lengths corresponding to these people", lengths)
print("\n here is the shortest path: ", nx.dijkstra_path(F,'U19886',ListOfGooglers[j]))
print("\n the second and third : ", nx.dijkstra_path(F,'U19886',ListOfGooglers[k]), nx.dijkstra_path(F,'U19886',ListOfGooglers[l]))


print("-------------------------------Now if we consider the number of people working at google known by our target at google----------")
lengths=[]
numberofpeopleknownatgoogle=[] #for each person, we will assign the number of people they know at google
for i in ListOfGooglers:
    length_of_the_shortest_path_to_i=nx.dijkstra_path_length(F,'U19886',i)
    c=0 #un compteur
    for j in G.neighbors(i):
        if j in ListOfGooglers: #we check if these people belong to google
            c+=1
    numberofpeopleknownatgoogle.append(c) 
    if c!=0:      
        lengths.append(length_of_the_shortest_path_to_i*(np.mean(numberofpeopleknownatgoogle)/c)**(1/5))
    else : 
        lengths.append(length_of_the_shortest_path_to_i)


print( " \n The number of people known at google is ",  trifusion(numberofpeopleknownatgoogle), " the average is " ,np.mean(numberofpeopleknownatgoogle), "and min and max are ", min(numberofpeopleknownatgoogle),",", max(numberofpeopleknownatgoogle))

j=lengths.index(min(lengths))#shortest path
lj=lengths[j]
lengths[j]+=5 # On augmente de manière artificielle afin de trouver les deuxièmes et troisièmes sans perturber l'indice
k=lengths.index(min(lengths))#second shortest
lk=lengths[k]
lengths[k]+=5 #idem que là haut
l=lengths.index(min(lengths))#third shortest
ll=lengths[l]
lengths[j]-=5
lengths[k]-=5


m=lengths.index(max(lengths))


print("The closest is : ",ListOfGooglers[j]," the length to ",ListOfGooglers[j]," is ", lengths[j]," while the avg length is ", np.mean(lengths)," and the longest is ", lengths[m] ,". \n The second and third lengths are ",lk,ll)
#print("\n And here are lengths corresponding to these people", lengths)
print("\n here is the shortest path: ", nx.dijkstra_path(F,'U19886',ListOfGooglers[j]))
print("\n the second and third : ", nx.dijkstra_path(F,'U19886',ListOfGooglers[k]), nx.dijkstra_path(F,'U19886',ListOfGooglers[l]))
print("\n Degree of U7024 :",G.degree('U7024'))