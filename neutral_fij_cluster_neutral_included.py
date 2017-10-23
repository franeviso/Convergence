
############################### INDIVIDUAL-BASED MODEL of NEUTRAL EVO-ECO MUTUALISTIC DYNAMICS ###################################
##################################################################################################################################
############################# Francisco Encinas-Viso 2012 
## Calculation of identical sites as in Melian et al (2012) PLoS Comp Biol.

import random 
import numpy as np
import copy
from sets import Set
import scipy.sparse as ss
from multiprocessing import Pool
import math
import networkx as nx
from sets import Set
import operator


random.seed(np.random.randint(1000))
def genome(gen):
    adn = [random.choice([-1,1]) for i in range(0,gen)]
    return adn 

def genome_with_neutral_gene(gen):
    adn = [random.choice([-1,1]) for i in range(0,gen)]
    adn.append(random.randint(10,20))
    return adn     


class indiv(object): 
  total_plant=0
  total_animal=0
  @staticmethod
  def reset():
      indiv.total_plant=0
      indiv.total_animal=0
  def __init__(self, sex, gen, bio):
      self.sex = sex
      self.gen = gen
      self.tag = False
      self._adn = genome_with_neutral_gene(gen)
      self.numgen = 1
      #self.indnumber = indiv.total
      self.trait = sum(self._adn[0:gen-1]) + random.gauss(0,1) + self.gen# alternative random.normalvariate(sum(self._adn),var) + random.random()
      self.bio = bio
      self.partners = []
      self.offspring = []
      self.species_n=0
      if bio is 'plant':
	indiv.total_plant += 1
	self.indnumber = indiv.total_plant
      else:
	indiv.total_animal +=1
	self.indnumber = indiv.total_animal
    
  def __call__(self):
      print self.adnid      
  def add_partner(self,x):
      self.partners.append(x)
  def add_offspring(self,x):
      self.offspring.append(x)         
  @property
  def adnid(self):
      return self._adn
  
 # lambda = 2, generate around 2-3 mutations
def mutate_pois(indgen,mutrate):
    x = np.random.poisson(mutrate*len(indgen),1)  # number of selected sites for mutation
    if x is not 0:
       adn = indgen
       y = len(adn)
       index_ = [random.randrange(0,y) for i in range(0,x)]
       # new_adn = [(adn[i] = random.randint(0,1)) for i in index_]
       last_ = y-1
       for i in index_:
		   if i < last_:
			   adn[i] = random.choice([-1,1])
		   elif i == last_:
			   adn[i] = random.randint(2,20)   
    else:
       adn=indgen
    return adn
    
  
def kill(ind):
  xx = range(0,len(ind)) # it should be N
  idx = random.choice(xx)  # pick an individual with respective index
  x = ind[idx]
  x.tag = True
  return x.sex, idx, x.indnumber

def gen_mean(Q,N):
    total = sum(sum(np.triu(Q)))
    total = total - len(np.diag(Q))
    mean = total/(((N*N)/2)-len(np.diag(Q)))
    return np.around(mean,3)

       
def genetic_matrices_fij(Q, Gdif):
    Q=np.where(Q>Gdif,1,0)
    return Q

def Q_matrices(vecs,genes):
    dime = len(vecs)
    Q = np.zeros((dime,dime))
    f = np.zeros((dime,dime))
    for i in range(len(vecs)):
      for j in range(len(vecs)):
		  Q[i,j] = float(np.dot(vecs[i]._adn,vecs[j]._adn))/genes
		  f[i,j] = (1+Q[i,j])/2
    return Q,f

def Q_matrices1(vecs,genes):
    dime = len(vecs)
    Q = np.zeros((dime,dime))
    f = np.zeros((dime,dime))
    for i in range(len(vecs)):
      for j in range(len(vecs)):
		  Q[i,j] = float(np.dot(vecs[i]._adn[0:genes-1],vecs[j]._adn[0:genes-1]))/genes
    return Q

def Q_matrices2(vecs,genes):
    dime = len(vecs)
    Q = [float(np.dot(vecs[i]._adn,vecs[j]._adn))/genes for i in range(len(vecs)) for j in np.arange(i,len(vecs)) ]  
    return np.array(Q)
    
    
def genetic_distance(ind1,ind2,genes):
  q = float(np.dot(ind1._adn[0:genes-1],ind2._adn[0:genes-1]))/genes
  f = (1+q)/2
  return f

def genetic_distance_bin(ind1,ind2,genes,Gdif):
  q = float(np.dot(ind1._adn,ind2._adn))/genes;f=0
  if (1+q)/2 > Gdif:
	  f=1
  return f  

def euclid(x1,y1,x2,y2):
	return math.sqrt( np.power(x1-x2,2) + np.power(y1-y2,2))
	
## Create matrix of euclidean distances
def matrix_eucl(X,Y):
	mat = np.zeros((len(X),len(Y)))
	for i in range(len(X)):
		for j in range(len(Y)):
			mat[i,j] = euclid(X[i,0],X[i,1],Y[j,0],Y[j,1])
	return mat/np.max(mat)  

## Matrix of distance intra-guild
def distance_intraguild(N,mean,sd):
  D = sd*np.random.randn(N,N) + mean     
  D = np.where(D<0,abs(random.gauss(mean,sd)),D)
  upper = np.triu(D)
  upper = np.where(upper>1,abs(random.gauss(mean,sd)),upper)
  for i in range(N):
    upper[i,i] =1
  flipped = np.rot90(np.fliplr(upper))
  d = flipped + upper
  d = np.where(d>1,1,d)
  d = np.around(d,3)
  return d

def distance_betweenguild(N,mean,sd):
  D = sd*np.random.randn(N,N) + mean     
  D = np.where(D<0,abs(random.gauss(mean,sd)),D)
  upper = np.triu(D)
  flipped = np.rot90(np.fliplr(upper))
  d = flipped + upper
  d=np.where(d>1,1,d)
  d = np.around(d,3)
  return d

def new_distance(distance,row,mean,sd):
  D = distance
  new = sd*np.random.randn(N,1) + mean     
  new = np.where(new<0,abs(random.gauss(mean,sd)),new)
  new = np.where(new>1,1,new)
  
  
  #### Calculation of speciation events 
def wraper(tup):
  return genetic_distance_bin(*tup)
  
def wraper2(tup):
  return genetic_distance(*tup)

def main2(xx,males,genes):
  pool = Pool(2)
  ans = pool.map(wraper2, [(xx,i,genes) for i in males])
  pool.close()
  pool.join()
  return ans
  
  
def main(vecs,N,genes,Gdif):
  pool = Pool(2)
  ans = [pool.map(wraper, [(vecs[i],vecs[k],genes,Gdif) for k in np.arange(i,N)]) for i in range(N) ]
  return ans
    
def matrix_transf(mat,n):
  mat2 = [item for sublist in mat for item in sublist]
  p=-1;mat1=ss.lil_matrix((n,n))
  for i in range(n):
    for j in np.arange(i,n):
       p +=1
       mat1[i,j] = mat2[p]
  return mat1

def vector_fij_males(xx,males,genes):
  vecf = main2(xx,males,genes)
  return vecf
  
def speciation(vecs,N,genes,Gdif):
  Qvec = main(vecs,N,genes,Gdif)
  Q = matrix_transf(Qvec,N)
  #Qt = ss.triu(Qs,k=1)               # obtain only upper triangular
  nspec,l = ss.cs_graph_components(Q)
  return nspec,l  
   

## Create N identical individuals
def clones(N,genes, bio,ptrait):
  # Clone DNA
  xdna = genome_with_neutral_gene(genes)
  while sum(xdna)+genes<ptrait:
    xdna= genome_with_neutral_gene(genes)
  # Create indiv vectors
  males = [indiv('male',genes, bio) for i in range(N/2)]
  females = [indiv('female',genes,bio) for i in range(N/2)]
  # Change ADN to clone's ADN
  for i in range(N/2):
    males[i]._adn = xdna
    males[i].trait = sum(males[i]._adn[0:genes-1]) + random.gauss(0,1) + genes
    females[i]._adn = xdna
    females[i].trait = sum(females[i]._adn[0:genes-1]) + random.gauss(0,1) + genes
  clones = males + females
  return clones

def clones_traits(N,dna,numgen, bio):
  # Clone DNA
  xdna = mutate(dna,numgen)
  genes = len(xdna)
  # Create indiv vectors
  males = [indiv('male',genes, bio) for i in range(N/2)]
  females = [indiv('female',genes,bio) for i in range(N/2)]
  # Change ADN to clone's ADN
  for i in range(N/2):
    males[i]._adn = xdna
    males[i].trait = sum(males[i]._adn) + random.gauss(0,1) + genes
    females[i]._adn = xdna
    females[i].trait = sum(females[i]._adn) + random.gauss(0,1) + genes
  clones = males + females
  return clones
  
    
## Asexual reproduction
def mating_asex(indivs, killed_index, killed_sex, killed_indnumber, numgen=1):
  indivs.pop(killed_index)
  x = random.choice(indivs)
  #numgen = 10
  parentadn = copy.copy(x._adn)
  babyadn = mutate(parentadn,numgen)
  p = indiv(killed_sex,len(babyadn))
  p._adn = babyadn
  p.indnumber = killed_indnumber
  return p
  


## Mating function for animals (includes mutation and recombination)

def matingsex_animal4(indivs, plants, killed_index, killed_sex, killed_indnumber, dif, N, numgen, Distance_aa, Distance_pa, genes,rec,dist):
  indivs.pop(killed_index)
  sample_again = False
  alls = [i for i in indivs if Distance_aa[killed_indnumber-1,i.indnumber-1]<dist ]
  mom = [i for i in alls if i.sex is 'female']
  while sample_again is False:
     x = random.choice(mom)    # choose a female
     males= [i for i in alls if i.sex is 'male'] # males vector
     if len(males) > 50:
        sample_males = random.sample(males,25)
     else:
        sample_males = males
     chosen_males = [i for i in sample_males if genetic_distance(x,i,genes)>dif]
     #print len(chosen_males)
     if len(chosen_males)==0:
        sample_again = False
     else:
        y = random.choice(chosen_males)
        sample_again = True
        break
  plant = [j for j in plants if Distance_pa[j.indnumber-1,x.indnumber-1] < dist]
  pp = random.choice(plant)
  #  Recombination - crossover: half male adn and half female adn
  cross = random.randint(0,genes-5) #crossover position
  neutral_x = x._adn[-1]
  neutral_y = y._adn[-1]
  halfy_1 = y._adn[0:cross]
  halfx_1 = x._adn[cross:genes-1]
  adn_1 = halfy_1 + halfx_1
  halfy_2 = y._adn[cross:genes-1]
  halfx_2 = x._adn[0:cross]
  adn_2 = halfx_2 + halfy_2
  neutral = random.choice([neutral_x,neutral_y])
  adnbaby_pre = random.choice([adn_1,adn_2])
  adnbaby = adnbaby_pre + [neutral]
  # Need to addd mutation to adnbaby
  babyadn = mutate_pois(adnbaby,numgen)
  p = indiv(killed_sex,len(babyadn),'animal')
  p._adn = babyadn
  p.trait = sum(babyadn) + random.gauss(0,1) + genes
  p.indnumber = killed_indnumber
  # Add mutualistic partners
  x.add_partner(pp.indnumber)
  y.add_partner(pp.indnumber)
  pp.add_partner(x.indnumber)
  x.add_offspring(p)
  y.add_offspring(p)
  # New distances in the animal-animal and plant-animal matrix
  # Distance_a[p.indnumber,] = 
  # Distance_pa[p.indnumber,] = 
  return p


# Mating function for plants (includes recombination, mutation)

def matingsex_plant4(indivs, animals, killed_index, killed_sex, killed_indnumber, dif, N, numgen, Distance_pp, Distance_pa,genes,rec,dist):
  indivs.pop(killed_index)
  sample_again = False
  alls = [i for i in indivs if Distance_pp[killed_indnumber-1,i.indnumber-1]<dist]
  mom = [i for i in alls if i.sex is 'female']
  while sample_again is False:
      x = random.choice(mom)    # choose a female
      males= [i for i in alls if i.sex is 'male'] # males vector
      if len(males) > 50:
        sample_males = random.sample(males,25)
      else:
        sample_males = males
      chosen_males = [i for i in sample_males if genetic_distance(x,i,genes)>dif]
      #print len(chosen_males)
      if len(chosen_males)==0:
        sample_again = False
      else:
        y = random.choice(chosen_males)
        sample_again=True
        break
  anim = [j for j in animals if Distance_pa[x.indnumber-1, j.indnumber-1] < dist]
  matching = [j for j in anim if j.trait >= x.trait]
  #print '\n animals matching: %d' % (len(matching))
  a = random.choice(anim)
  #  crossover: half male adn and half female adn
  cross = random.randint(0,genes-5) #crossover position
  neutral_x = x._adn[-1]
  neutral_y = y._adn[-1]
  halfy_1 = y._adn[0:cross]
  halfx_1 = x._adn[cross:genes-1]
  adn_1 = halfy_1 + halfx_1
  halfy_2 = y._adn[cross:genes-1]
  halfx_2 = x._adn[0:cross]
  adn_2 = halfx_2 + halfy_2
  neutral = random.choice([neutral_x,neutral_y])
  adnbaby_pre = random.choice([adn_1,adn_2])
  adnbaby = adnbaby_pre + [neutral]
  #rint "Neutral allele: %d" % neutral
  # Need to addd mutation to adnbaby
  babyadn = mutate_pois(adnbaby,numgen)
  p = indiv(killed_sex,len(babyadn),'plant')
  p._adn = babyadn
  p.trait = sum(babyadn) + random.gauss(0,1) + genes
  p.indnumber = killed_indnumber
  # Add mutualistic partners
  x.add_partner(a.indnumber)
  y.add_partner(a.indnumber)
  a.add_partner(x.indnumber)
  x.add_offspring(p)
  y.add_offspring(p)
  return p



def nodfun(M):
  n = np.shape(M)[0]; m = np.shape(M)[1]
  Mtrow = [sum(M[i,:]) for i in range(n)]
  Mtcol = [sum(M[:,j]) for j in range(m)]
  Npaired_row = np.zeros((n,n))
  Npaired_col = np.zeros((m,m))
  for i in range(n):
    for l in range(n):
        if Mtrow[i] > Mtrow[l] and Mtrow[l]!=0:
            Npaired_row[i,l] =  100*(float(sum(M[i,:]*M[l,:])/(Mtrow[l])))
        else:
            Npaired_row[i,l] = 0
  for i in range(m):
    for l in range(m):
        if Mtcol[i] > Mtcol[l] and Mtcol[l]!=0:
            Npaired_col[i,l] =  100*(float(sum(M[:,i]*M[:,l])/(Mtcol[l])))
        else:
            Npaired_col[i,l] = 0
  NODFr = [sum(Npaired_row[i,i:n]) for i in range(n-1)]
  NODFc = [sum(Npaired_col[j,j:m]) for j in range(m-1)]   
  Nsumr= sum(NODFr); Nsumc = sum(NODFc)
  nodf= float((Nsumr + Nsumc))/((n*(n-1)/2) + (m*(m-1)/2))
  return nodf   
    
def mean_traits(traits,nodes):
  meansp=[]
  for i in range(len(nodes)):
    vectraits = [traits[j] for j in nodes[i]]
    meansp.append(np.mean(vectraits))
  return meansp

def var_traits(traits,nodes):
  vari=[]
  for i in range(len(nodes)):
    vectraits = [traits[j] for j in nodes[i]]
    vari.append(np.std(vectraits))
  return vari

def partnersfun(partner_list,nodes):
   partners=[]; spx=[]
   for i in range(len(nodes)):  #species
     row=[]
     for j in nodes[i]:         # individuals
       row.append(partner_list[j])
     partners.append(row)    
   for i in range(len(nodes)):
     y = Set(partners[i][0])  
     for j in range(len(nodes[i])):
       x = Set(partners[i][j])
       y.update(y.union(x))
     spx.append(y)
   return spx
   
   
# Species interaction matrix 
def matrix_species(nodes_p,nodes_a,species_partners):
   matrix_ap = np.zeros((len(nodes_p),len(nodes_a)))  
   new = np.zeros((len(nodes_p),len(nodes_a)))  
   for i in range(len(nodes_p)):
     for j in range(len(nodes_a)):
        if species_partners[i] & Set(nodes_a[j]):
          matrix_ap[i,j]=1
          new[i,j] = len(species_partners[i] & Set(nodes_a[j]))
   return matrix_ap,new
   

def convergence(G,P,th):
	n = len(G);conv = 0
	for i in range(n):
		G[i,i] = 0
		P[i,i] = 0
	sister=[]
	for i in range(n):
		m = np.max(G[i,:])
		maxval = [k for k, j in enumerate(G[i,:]) if j == m]
		sister.append(maxval[0])
	phensim=[]
	for i in range(n):
		sim = [k for k, j in enumerate(P[i,:]) if j > th]
		phensim.append(sim)	
	for i in range(n):
		if phensim[i]:
			for k in phensim[i]:
				if k != sister[i]:
					conv += 1
	return conv
## @ ----> Important: species are sorted by abundance in descending order
# Provide number of convergence events per species
# This shows convergence in terms of percentage, so it is considered the GUILD size
	
def convergence2(G,P,th):
	n = len(G);convs=[]
	for i in range(n):
		G[i,i] = 0
		P[i,i] = 0
	sister=[]
	for i in range(n):
		m = np.max(G[i,:])
		maxval = [k for k, j in enumerate(G[i,:]) if j == m]
		sister.append(maxval[0])
	phensim=[]
	for i in range(n):
		sim = [k for k, j in enumerate(P[i,:]) if j > th]
		phensim.append(sim)
	for i in range(n):
		if phensim[i]:
			conv_per_sp=0
			for k in phensim[i]:
				if k != sister[i]:
					conv_per_sp+=1
			convs.append(float(conv_per_sp)/(n-2))
	return convs

def convergence3(G,P):
	n = len(G);convs=[]
	for i in range(n):
		G[i,i] = 0
		P[i,i] = 0
	sister=[]
	for i in range(n):
		m = np.max(G[i,:])
		maxval = [k for k, j in enumerate(G[i,:]) if j == m]
		sister.append(maxval[0])
	for i in range(n):
	    conv=0
	    for k in range(n):
	        if P[i,k] > P[i,sister[i]]:
	            conv+=1
	    print conv        
	    convs.append(float(conv)/(n-2))              	      
	return convs		
	
def species_nodes(Gm):
  Gmat = np.triu(Gm) 
  N = Gmat.shape[0]
  diagonal = zip(range(N),range(N))
  for i,j in diagonal:
     Gmat[i,j] = 0   
  indices  = np.where(Gmat>0)
  edges = zip(indices[0],indices[1])
  G = nx.Graph()
  G.add_nodes_from(range(N))
  G.add_edges_from(edges)
  node_sets=nx.connected_components(G)
  return node_sets	      


def species(Gmat,N):
  Gmat1 = np.triu(Gmat) 
  diagonal = zip(range(N),range(N))
  for i,j in diagonal:
    Gmat1[i,j] = 0  
  indices  = np.where(Gmat1>0)
  edges = zip(indices[0],indices[1])
  G = nx.Graph()
  G.add_nodes_from(range(N))
  G.add_edges_from(edges)
  sps = nx.number_connected_components(G)
  return sps


def phensim(traits):
  maxt= max(traits)
  N = len(traits)
  normtrait=np.zeros((N,N))
  for i in range(N):
    for j in range(N):
       normtrait[i,j]= 1 - abs(traits[i]-traits[j])/maxt # remove abs() to check dispersion
       #normtrait[i,j]= 1 - ((traits[i]-traits[j])/maxt)
  return normtrait  


def phensim_pa(traitsp,traitsa):
  maxtp= max(traitsp); maxta=max(traitsa)
  NP = len(traitsp); NA=len(traitsa)
  phensimpa=np.zeros((NP,NA))
  for i in range(NP):
    for j in range(NA):
       phensimpa[i,j]= 1 - abs(float(traitsp[i])/maxtp - float(traitsa[j])/maxta) # remove abs() to check dispersion
  return phensimpa

def complementarity(meantraitplants, meantraitanimals,th):
    matrixphensimpa = phensim_pa(meantraitplants,meantraitanimals)
    x = np.shape(matrixphensimpa);comp=[]
    for i in range(x[0]):
        y = np.where(matrixphensimpa[i,:]>th)
        comp.append(float(len(y[0]))/x[1])
    return comp

## Neutral marker functions

def cavallisf(x,y, number_of_alleles):
	sumfreq=0
	for i in range(number_of_alleles):
		sumfreq += np.sqrt(x[i]*y[i])
	return 2*np.sqrt(2)*(1 - sumfreq)/(math.pi)	


def frequency_sp(alleles_neutral_list,alleles_neutral):
	allele_neutral_freqs = []
	alleles_neutral_list2 = alleles_neutral_list
	for i in range(len(alleles_neutral)):#species
		alleles_freqs = np.zeros((len(alleles_neutral_list2),1))
		for al in range(len(alleles_neutral_list2)):
			counts=0 
			for k in alleles_neutral[i]:
				if k == alleles_neutral_list2[al]:
					counts += 1
			alleles_freqs[al] = float(counts)/len(alleles_neutral[i]);
		allele_neutral_freqs.append(alleles_freqs)
	return allele_neutral_freqs				
        
def matrix_geneticdist_neutral(allele_neutral_freqs,number_of_alleles):
	species = len(allele_neutral_freqs)
	genmat_neutral = np.zeros((species,species))
	for i in range(species):
		for j in range(species):
			genmat_neutral[i,j] = cavallisf(allele_neutral_freqs[i],allele_neutral_freqs[j],number_of_alleles)
	return genmat_neutral
	
	        
        
def build_neutral_alleles_list(individuals,sortednodes):
	alleles_neutral = []
	for k in range(len(sortednodes)):
		alleles_neutral_sp=[]
		for i in sortednodes[k]:
			alleles_neutral_sp.append(individuals[i]._adn[-1])
		alleles_neutral.append(alleles_neutral_sp)
	alleles_neutral_list = Set([alleles_neutral[0][0]])
	for i in alleles_neutral:
		for k in i:
			if alleles_neutral_list not in Set([k]):
				alleles_neutral_list.union_update(Set([k]))
	return alleles_neutral, list(alleles_neutral_list)
	
