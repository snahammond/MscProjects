# This programme is to handle the Final for Advanced Algorithm and 
# computational models 1
# 
# Developed by:
#  Samuel Nii Armah Hammond
# For the supervision of:
# 	Prof. Giacomo Fiumara

##this import is used for the growth
##this is for the calculations
##this is for drawing the graph
import random
import math
import matplotlib.pylab as plt

#Create an initial network composed of a limited number of nodes m.0
#m.0 is also the initial number of links
#therefore m.0 nodes will have n.0 links (a simple way to make the network) 
#completely wired; no isolated node is allowed
#empty list to hold the edgeList representation of the graph
#return the graph in EdgeList format
def makeAwiredGraphInEdgeListFormat(m0):
 graphInEdgeListFormat = [] 
 linksCounter=0
 for x in list(range(m0)): 
  for y in list(range(x,m0)):
   if x != y and linksCounter < m0:
    linksCounter=linksCounter+1
    graphInEdgeListFormat.append([x,y])
 	
 return graphInEdgeListFormat


#now create the initial network in a dictionary format
def createInitialNetworkInDF(initialWiredNetwork):
 networkInDF = {}
 for edge in initialWiredNetwork:  
  if edge[0] not in networkInDF.keys():   
   key = edge[0]
   value = [edge[1]]
   networkInDF.update({key:value})
  else:   
   key = edge[0]
   value = edge[1]
   networkInDF[key].append(value)
  
  if edge[1] not in networkInDF.keys():   
   key = edge[1]
   value = [edge[0]]
   networkInDF.update({key:value})
  else:   
   key = edge[1]
   value = edge[0]
   networkInDF[key].append(value)
 return networkInDF


##insert a new edge only if the edge u,v does not exist
def insert_edge(networkInDF,u,v):
 
 if u in networkInDF.keys():  
  if v not in networkInDF[u]:   
   networkInDF[u].append(v)   
   if v in networkInDF.keys():
    networkInDF[v].append(u)
   else:
    networkInDF.update({v:[u]})   
   newlyCreatedEdge = [u,v]  
   return newlyCreatedEdge
  else:   
   return [None]
 elif v in networkInDF.keys():  
  if u not in networkInDF[v] :   
   networkInDF[v].append(u)   
   if u in networkInDF.keys() :
    networkInDF[u].append(v)
   else:       
    networkInDF.update({u:[v]})   
   newlyCreatedEdge = [u,v]
   return newlyCreatedEdge
  else:   
   return [None]
 else:  
  return [None]

#iterative Process to add new nodes till N  
def iterativeProcessAddingNewNode(initialNetworkInDF,m,N):
 counter = m
 networkToWorkOn = initialNetworkInDF
 
 timeAndAverageDegree = {}
 
 timeStep = 0
 
 keyI = timeStep
 valueI = round((2*m*(timeStep+1))/len(initialNetworkInDF.keys()))
 timeAndAverageDegree.update({keyI:valueI})
 
 while(counter<N):  
  timeStep = timeStep +1
  
  currentNodesInNetwork = networkToWorkOn.keys()  
  selectNodes = random.sample(currentNodesInNetwork,m)
  
  for luckyNode in selectNodes:
   insert_edge(networkToWorkOn,counter,luckyNode)   
  
  keyI = timeStep
  valueI = (2*m*(timeStep+1))/len(networkToWorkOn.keys())
  
  timeAndAverageDegree.update({keyI:valueI})   
  counter = counter+1
  
 return (networkToWorkOn,timeAndAverageDegree)
  
def getNodeAndDegreeDic(newlyCreatedBarabasiAlNWT):
 nodeDegreeDict = {}  
 for node in newlyCreatedBarabasiAlNWT.keys():
  key = node
  value = len(newlyCreatedBarabasiAlNWT[node])
  nodeDegreeDict.update({key:value})
 return nodeDegreeDict
 
 
def sortMyDictry (degreeNodeDistr):
 kl=list(degreeNodeDistr.keys())
 new_kl = kl + []
 new_kl.sort()
 new_d = {}
 for k in new_kl:
  new_d[k] = degreeNodeDistr[k]
 return new_d
 

def measureDegDistr(newlyCreatedBarabasiAlNWT,N,m):
  degreeProbNodeDistr = {}
  degreeNodeDistr = {}
  nodeDegreeDict = getNodeAndDegreeDic(newlyCreatedBarabasiAlNWT)
  
  for nodeDegree in nodeDegreeDict.keys():
   if nodeDegreeDict[nodeDegree] not in degreeNodeDistr.keys():
    key = nodeDegreeDict[nodeDegree]
    value = 1
    degreeNodeDistr.update({key:value})
   else:
    key = nodeDegreeDict[nodeDegree]
    value = 1
    degreeNodeDistr[key] = degreeNodeDistr[key] + value	
  
  degreeNodeDistr = sortMyDictry(degreeNodeDistr) 
  
  for degree in degreeNodeDistr.keys():
   key = degree
   f = degreeNodeDistr.get(key)   
   
   probOfDegree = f/N 
   value = probOfDegree
   degreeProbNodeDistr.update({key:value})
  
  return (degreeProbNodeDistr)

#function to calculate the theoritical Degree Dist line  
#find smallest y value
#find the smallest deg
def measureDegDistrT(nodeDegreeProbDistr,N,m): 
  degreeProbNodeDistr = {}
  
  allProb = list(nodeDegreeProbDistr.values())
  minProb = min(allProb)  
  
  allDeg = list(nodeDegreeProbDistr.keys())
  minDeg = min(allDeg)
  
  probOfDegree = 1
  deg = minDeg
  while probOfDegree > minProb:
   probOfDegree = math.exp(1)/m * (math.exp(-(deg/m)))
   degreeProbNodeDistr.update({deg:probOfDegree})
   deg += 1

  return (degreeProbNodeDistr)

#find the theoretical results of k(t) for node i such that this node i arrived at time 0
def measureAverageDegForIatT(m,N):
 timeAndAverageDegreeDictT = {}
 
 for t in range(N):
  key = t
  value = m * math.log( math.exp(1) * ((m+t-1)/(m+0-1)) )
   
  timeAndAverageDegreeDictT.update({key:value})
 
 return timeAndAverageDegreeDictT
   
  

##run main to show implemented functions
#main data we will use to draw our graph
netWorkInDF,newlyCreatedBarabasiAlNWT = {},{}
timeAndAverageDegreeList =[]
nodeDegreeProbDistrList=[]
m0 = []
for times in range(3):
 
 if(times == 0):
  m0.append(int(input("Please input the initial number of links: ")))
 else:
  m0.append(int(input("Please input another initial number of links: ")))

N= int(input("Please input the maximum number of nodes: ")) 
index = 0
if(N>1000):
 

	for m0i in m0:
	 
	 if ( m0i>2):
	  
	  print('---------------------Creating Initial Wired Network----------------------------------')
	  initialWiredNetwork = makeAwiredGraphInEdgeListFormat(m0i)
	  
	  print('---------------------Create network in a dictionary format we can work with easily----')
	  
	  netWorkInDF = createInitialNetworkInDF(initialWiredNetwork)
	  
	  print('---------------------End of Initial wired Network----------------------------------') 
	  
	  
	  print('---------------------Creating barabasi Albert Network with ', N , ' nodes added iteratively--------')
	  newlyCreatedBarabasiAlNWT,timeAndAverageDegree = iterativeProcessAddingNewNode(netWorkInDF,m0i,N)	  
	  
	  timeAndAverageDegreeList.append(timeAndAverageDegree)	  
	  
	  
	  nodeDegreeProbDistr = measureDegDistr(newlyCreatedBarabasiAlNWT,N,m0i)
	  nodeDegreeProbDistrList.append(nodeDegreeProbDistr)
	   
	  list_key_value = [ [k,v] for k, v in nodeDegreeProbDistr.items() ]
	  x, y = zip(*list_key_value)
	  
	  if index == 0:
	   
	   labeli = str(m0i)
	   style = 'dotted'
	   plt.plot(x, y, label =labeli,ls= style)
	   plt.yscale('log')
	   
	  elif index ==1:
	   
	   labeli = str(m0i)
	   style = 'dotted'
	   plt.plot(x, y, label =labeli,ls= style )
	   plt.yscale('log')
	  else:
	   
	   labeli = str(m0i)
	   style = 'dotted'
	   plt.plot(x, y, label =labeli,ls= style)
	   plt.yscale('log')   
	  
	 
	 else:
	  raise ValueError('invalid parameter provided, Initial links should be greater than 2') 

	 index = index + 1  
	
	nodeDegreeProbDistr = measureDegDistrT(nodeDegreeProbDistrList[0],N,m0[0],)
	 
	list_key_value = [ [k,v] for k, v in nodeDegreeProbDistr.items() ]
	x, y = zip(*list_key_value)
	labeli = str(m0[0]) +'T'
	style = 'dotted'
	plt.plot(x, y, label =labeli,ls= style )
	 
	plt.legend(loc='best') 
	plt.show()
	

	

	for timeAndAverageDegreeDict in timeAndAverageDegreeList:
	 
	 list_key_value = [ [k,v] for k, v in timeAndAverageDegreeDict.items() ]
	 x, y = zip(*list_key_value)
	 style = 'dotted'
	 plt.plot(x, y, label ='',ls= style )
	 plt.xscale('log')
	
    
	timeAndAverageDegreeDictT =measureAverageDegForIatT(m0[0],N)
	
	list_key_value = [ [k,v] for k, v in timeAndAverageDegreeDictT.items() ]
	x, y = zip(*list_key_value)
	style = 'dotted'
	plt.plot(x, y, label ='',ls= style )
	plt.xscale('log')
	 
	plt.legend(loc='best') 
	plt.show()
	
else:
 print('Sorry, please enter N > 1000')





