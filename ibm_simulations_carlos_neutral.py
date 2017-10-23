############################### INDIVIDUAL-BASED MODEL of NEUTRAL EVO-ECO MUTUALISTIC DYNAMICS ###################################
##################################################################################################################################
############################# Francisco Encinas-Viso 2012 

## Main outputs
# first line will include parameter values
# Matrix Q of genetic distances
#   - write file with matrix: q_matrix.dat
# Individual initial and final trait values
#   - write matrix initial-final trait values
# Individual's mutualistic partners:
#   - write file mutualistic partners or use pickle

## import functions
from time import time, asctime
from  neutral_fij_cluster_neutral_included import *
import sys
import csv

############ simulations  ##############
# Initial population of N individuals 
N=2000
genes =150
#genetic distance
dif = 0.97
sims = 1
## Parameters
generations = 4000*N
# lambda poisson
#numgen = 0.15
# Threshold of convergence (genetic and phenotypic similarity)
th = 0.97



simn = sys.argv[1] # string id of the simulation
dif = float(sys.argv[2]) # genetic similarity minimum 
dist = float(sys.argv[3]) # maximum geographic distance
rec = float(sys.argv[4]) # recombination prob
numgen = float(sys.argv[5]) # Mutation rate (Poisson lambda)
sims = int(sys.argv[6]) # Number of replicates


out_filename = "sim_ibm_"
log_filename = out_filename + simn + ".log"
log_file = open(log_filename, "w")
log_file.write ("Mutation rate: %f\n" % numgen )
log_file.write ("Generations: %d\n" % generations )
log_file.write("Launched: %s\n" % asctime())
log_file.close()

# Files for output:

file_spec_abundance = open("file_spec_abundance_%s_%d.csv" % (simn,sims),"wb")
file_matgen_anim_spec = open("file_matgen_anim_spec_%s_%d.csv" % (simn,sims),"wb")
file_matphen_anim_spec = open("file_matphen_anim_spec_%s_%d.csv" % (simn,sims),"wb")
file_matgen_plants_spec = open("file_matgen_plants_spec_%s_%d.csv" % (simn,sims),"wb")
file_matphen_plants_spec = open("file_matphen_plants_spec_%s_%d.csv" % (simn,sims),"wb")
file_matgen_neutral_anim_spec = open("file_matgen_neutral_anim_spec_%s_%d.csv" % (simn,sims),"wb")
file_matgen_neutral_plants_spec = open("file_matgen_neutral_plants_spec_%s_%d.csv" % (simn,sims),"wb")
file_traits_p = open("file_traits_p_%s_%d.csv" % (simn,sims),"wb")
file_traits_a = open("file_traits_a_%s_%d.csv" % (simn,sims),"wb")
file_matrixpa = open("file_matrixpa_%s_%d.csv" % (simn,sims),"wb")
file_nestedness = open("file_nestedness_%s_%d.csv" % (simn,sims),"wb")
file_complementarity = open("file_complementarity_%s_%d.csv" % (simn,sims),"wb")
file_neutral_alleles_p = open("file_neutral_alleles_p_%s_%d.csv" % (simn,sims),"wb")
file_neutral_alleles_a = open("file_neutral_alleles_a_%s_%d.csv" % (simn,sims),"wb")

t_start=time()

for k in range(sims):
      if k>0:
         del plants[:]
         del animals[:]
         indiv.reset()       
      P_coords = np.random.uniform(0,1,(N,2))
      A_coords = np.random.uniform(0,1,(N,2))        
      Distance_pp = matrix_eucl(P_coords,P_coords)#distance_intraguild(N,mean,sd)
      Distance_aa = matrix_eucl(A_coords,A_coords)#distance_intraguild(N,mean,sd)
      Distance_pa = matrix_eucl(P_coords,A_coords)#distance_betweenguild(N,mean,sd)
      plants = clones(N,genes,'plant',0)
      animals = clones(N,genes,'animal',plants[0].trait)
      traits_plants = [plants[i].trait for i in range(N)]  
      traits_animals = [animals[i].trait for i in range(N)]
      animal_speciation=[]  
      # synchronic 
      
      for i in range(generations):
		  #start=time()
          dead_a = kill(animals)
          new_a = matingsex_animal4(animals, plants, dead_a[1], dead_a[0], dead_a[2], dif, N, numgen,Distance_aa,Distance_pa,genes,rec,dist)
          animals.append(new_a)
          dead = kill(plants) 
          new = matingsex_plant4(plants, animals, dead[1], dead[0], dead[2], dif, N, numgen,Distance_pp,Distance_pa,genes,rec,dist)
          plants.append(new)
          #print new._adn
          #print i
          #elapsed = time() - start
          #if i%500000 ==0:
             #animal_speciation.append(speciation(animals,N,genes,dif))
             #elapsed = time.time() - start
             #print animal_speciation
           
      traits_plants_evol = [plants[i].trait for i in range(len(plants))]  
      traits_animals_evol = [animals[i].trait for i in range(len(animals))]   
      ppartners = []
      apartners = []
      numberp = []
      numbera = []
      for i in range(N):
          ppartners.append(plants[i].partners)
          apartners.append(animals[i].partners)
          numberp.append(plants[i].indnumber)
          numbera.append(animals[i].indnumber)
      plant_genmat= Q_matrices1(plants,genes)
      animal_genmat = Q_matrices1(animals,genes)
      plant_genmat2= genetic_matrices_fij(plant_genmat,dif)
      animal_genmat2 = genetic_matrices_fij(animal_genmat,dif)
      qmatp = plant_genmat
      qmata = animal_genmat
      plant_species = species(plant_genmat2,N)
      animal_species = species(animal_genmat2,N)
      nodesp = species_nodes(plant_genmat2)
      nodesa = species_nodes(animal_genmat2)
      nodes_p=[];nodes_a=[]
      plants_abundance=[];animals_abundance=[]
      for i in nodesp:
          plants_abundance.append(len(i))
          nodes_p.append(i)
      for j in nodesa:
          animals_abundance.append(len(j))
          nodes_a.append(j)
      file_spec_abundance.write("Plants: \n")
      np.savetxt(file_spec_abundance, plants_abundance, delimiter=",",fmt= '%4.2f')
      file_spec_abundance.write("Animals: \n")
      np.savetxt(file_spec_abundance, animals_abundance, delimiter=",",fmt= '%4.2f')	  
      pmatf= phensim(traits_plants_evol)
      pmatfa = phensim(traits_animals_evol)
      sortedplants=sorted(zip(plants_abundance,range(len(nodes_p))),reverse=True)
      sortednodesplants=[]
      for i,j in sortedplants:
          sortednodesplants.append(nodes_p[j])
      sortedanimals=sorted(zip(animals_abundance,range(len(nodes_a))),reverse=True)
      sortednodesanimals=[]
      for i,j in sortedanimals:
          sortednodesanimals.append(nodes_a[j])
		  		  

      meangen = np.zeros((len(nodes_a),len(nodes_a)))
      meanphen = np.zeros((len(nodes_a),len(nodes_a)))
      for k in range(len(nodes_a)):
          for h in range(len(nodes_a)):
              meangen[k,h] = np.mean([qmata[i,j] for i in sortednodesanimals[k] for j in sortednodesanimals[h]])
              meanphen[k,h] = np.mean([pmatfa[i,j] for i in sortednodesanimals[k] for j in sortednodesanimals[h]])
      meangenp = np.zeros((len(nodes_p),len(nodes_p)))
      meanphenp = np.zeros((len(nodes_p),len(nodes_p)))
      for k in range(len(nodes_p)):
          for h in range(len(nodes_p)):
              meangenp[k,h] = np.mean([qmatp[i,j] for i in sortednodesplants[k] for j in sortednodesplants[h]])
              meanphenp[k,h] = np.mean([pmatf[i,j] for i in sortednodesplants[k] for j in sortednodesplants[h]])
              
      neutral_alleles_list_plants = build_neutral_alleles_list(plants,sortednodesplants)
      neutral_alleles_list_animals = build_neutral_alleles_list(animals,sortednodesanimals)
      
      allele_neutral_freqs_plants = frequency_sp(neutral_alleles_list_plants[1],neutral_alleles_list_plants[0])
      allele_neutral_freqs_animals = frequency_sp(neutral_alleles_list_animals[1],neutral_alleles_list_animals[0])
      
      genmat_neutral_plants = matrix_geneticdist_neutral(allele_neutral_freqs_plants,len(neutral_alleles_list_plants[1]))
      genmat_neutral_animals = matrix_geneticdist_neutral(allele_neutral_freqs_animals,len(neutral_alleles_list_animals[1]))        
			  
      np.savetxt(file_matgen_anim_spec, meangen, delimiter=",",fmt= '%4.2f')
      wr2 = csv.writer(file_matgen_anim_spec)
      wr2.writerow("\n")
      np.savetxt(file_matphen_anim_spec, meanphen, delimiter=",",fmt= '%4.2f')
      wr2 = csv.writer(file_matphen_anim_spec)
      wr2.writerow("\n")
      np.savetxt(file_matgen_neutral_anim_spec, genmat_neutral_animals, delimiter=",",fmt= '%4.2f')
      wr2 = csv.writer(file_matgen_neutral_anim_spec)
      wr2.writerow("\n")
	  
      np.savetxt(file_matgen_plants_spec, meangenp, delimiter=",",fmt= '%4.2f')
      wr2 = csv.writer(file_matgen_plants_spec)
      wr2.writerow("\n")
      np.savetxt(file_matphen_plants_spec, meanphenp, delimiter=",",fmt= '%4.2f')
      wr2 = csv.writer(file_matphen_plants_spec)
      wr2.writerow("\n")
      np.savetxt(file_matgen_neutral_plants_spec, genmat_neutral_plants, delimiter=",",fmt= '%4.2f')
      wr2 = csv.writer(file_matgen_neutral_plants_spec)
      wr2.writerow("\n")
	  
      convergence_events_animals = convergence3(genmat_neutral_animals,meanphen)
      convergence_events_plants = convergence3(genmat_neutral_plants,meanphenp)
      spp = partnersfun(ppartners,sortednodesplants)
      spa = partnersfun(apartners,sortednodesanimals)
      matrixpa = matrix_species(sortednodesplants,sortednodesanimals,spp)
      nestedness = nodfun(matrixpa[0])
      meantraitplants = mean_traits(traits_plants_evol,sortednodesplants)
      vartraitplants = var_traits(traits_plants_evol,sortednodesplants)
      traits_stat_plants = np.array((meantraitplants,vartraitplants))
      vartraitanimals = var_traits(traits_animals_evol,sortednodesanimals)
      meantraitanimals = mean_traits(traits_animals_evol,sortednodesanimals)
      complementarity_events = complementarity(meantraitplants,meantraitanimals,th)
      trait_stat_animals = np.array((meantraitanimals,vartraitanimals))
      np.savetxt(file_matrixpa, matrixpa[0], delimiter=",",fmt= '%4.2f')
      wr2 = csv.writer(file_matrixpa)
      wr2.writerow("\n")
      np.savetxt(file_traits_p, traits_stat_plants, delimiter=",",fmt= '%4.2f')
      file_traits_p.write("Convergence events: \n")
      np.savetxt(file_traits_p, convergence_events_plants, delimiter=",",fmt= '%4.2f')
      file_traits_p.write("\n")
      np.savetxt(file_traits_a, trait_stat_animals, delimiter=",",fmt= '%4.2f')
      file_traits_a.write("Convergence events: \n")
      np.savetxt(file_traits_a, convergence_events_animals, delimiter=",",fmt= '%4.2f')
      file_traits_a.write("\n")
      file_nestedness.write("Nestedness: ")
      file_nestedness.write(str(nestedness))
      file_nestedness.write("\n")
      file_complementarity.write("Complementarity events: \n")
      np.savetxt(file_complementarity, complementarity_events, delimiter=",",fmt= '%4.2f')
      file_complementarity.write("\n")	  
      file_neutral_alleles_p.write("List of neutral alleles - plants: \n")
      np.savetxt(file_neutral_alleles_p, neutral_alleles_list_plants[1], delimiter=",",fmt= '%4.2f')
      file_neutral_alleles_p.write("\n")		  
      file_neutral_alleles_a.write("List of neutral alleles - animals: \n")
      np.savetxt(file_neutral_alleles_a, neutral_alleles_list_animals[1], delimiter=",",fmt= '%4.2f')
      file_neutral_alleles_a.write("\n")      	      


      

      
   ## importante guardar los ind_number de cada individuo
file_spec_abundance.close()
file_matgen_anim_spec.close()
file_matphen_anim_spec.close()
file_matgen_plants_spec.close()
file_matphen_plants_spec.close()
file_matgen_neutral_plants_spec.close()
file_matgen_neutral_anim_spec.close()
file_traits_p.close()
file_traits_a.close()
file_matrixpa.close()
file_nestedness.close()
file_complementarity.close()
file_neutral_alleles_p.close()
file_neutral_alleles_a.close()
 
t_end = time()  
# Duration time
log_file = open(log_filename, "a")
log_file.write("Duration: %f" % int(t_end - t_start) )
log_file.close()


