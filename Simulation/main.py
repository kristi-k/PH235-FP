import csv
import math
import random
import pylab
from collections import defaultdict

R = 8.3144598E-3

gReactions= []

class Reaction(object):
    def __init__(self, row):
        self.reaction = row['Reaction']
        self.order = int(row['Order'])
        self.reactants = row['Rall'].split(' | ')
        self.reactants = [x.strip() for x in self.reactants] #remove whitespace
        self.products = row['Pall'].split(' | ')
        self.products = [x.strip() for x in self.products] #remove whitespace
        self.Tmin = float(row['Tmin'])
        self.Tmax = float(row['Tmax'])
        self.A = float(row['A'])
        try:        
            self.n = float(row['n'])
        except:
            self.n = 0.0 #Original equation n = 0
        self.Ea = float(row['Ea'])
        
    def rate_constant(self, temp):
        global R
        return self.A * (temp / 298.0)**self.n * math.exp(-self.Ea/(R*temp))
        

def read_reaction_data(filename):
    min_rate_constant = 1E-25
    max_rate_constant = 1E3    
    #Reaction	Order	Rall	R1	R2	Pall	P1	P2	P3	P4	Tmin	Tmax	A	n	Ea    
    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        global gReactions        
        for row in reader:
            r = Reaction(row)
            try:
                rc = r.rate_constant(r.Tmax)
            except OverflowError as e:
                continue
            if rc < min_rate_constant or rc > max_rate_constant:
                continue
            gReactions.append(r)            
            #print(row['Reaction'], row['Ea'])
        #end for loop
    #end with
    print('Loaded %d reactions' % len(gReactions))

def read_user_input():
    global gReactions
    
    r1 = input('Enter Reactant: ').strip()
    #r1 = 'ClO4'
    result_set = list( filter(lambda x: r1 in x.reactants, gReactions) )
    if not result_set:
        print('No reactions contain %s as a reactant' % r1)
        return
    print('%d reactions contain reactant: %s' % (len(result_set), r1))

    if len(result_set) < 2:
        print('At least two reactions need to be found. Quitting...')
        return
        
    #Print
    counter = 0
    for rxn in result_set:
        counter += 1
        print('%d. %s Tmax=%0.0f' %(counter, rxn.reaction, rxn.Tmax))
    
    picks = input('Pick 2 reactions by entering two numbers: ').strip()
    rxn_numbers = picks.split(' ')
    
    R1 = result_set[ int(rxn_numbers[0]) - 1 ]
    R2 = result_set[ int(rxn_numbers[1]) - 1 ]
    
    print('Reaction1: %s Tmax=%0.0f \nReaction2: %s Tmax=%0.0f' %(R1.reaction, R1.Tmax, R2.reaction, R2.Tmax))
    
    Temp = int( input('Choose a temperature in the range %d - %d: ' % (min(R1.Tmin, R2.Tmin), max(R1.Tmax, R2.Tmax))))
    
    #Gillespie(R1, R2, Temp)
    Runge_Kutta(R1, R2, Temp)
    return
    
def Gillespie(R1, R2, Temp):
    #Seed random by time
    random.seed()

    time_max = 1000 #time of simulation in seconds
    initial_population = 200
    
    k1 = R1.rate_constant(Temp)
    k2 = R2.rate_constant(Temp)
    #print(k1, k2)
    
    #Set initial populations size
    population_sizes = dict()
    for reactant in R1.reactants:
        population_sizes[reactant] = [ initial_population ]
    for reactant in R2.reactants:
        population_sizes[reactant] = [ initial_population ]
    
    for product in R1.products:
        population_sizes[product] = [ 0 ]
    for product in R2.products:
        population_sizes[product] = [ 0 ]
    
    #print(population_sizes)
    
    time = 0.0
    tpoints = [0.0]
    while time < time_max:
        rate1 = k1        
        h1 = 1        
        for reactant in R1.reactants:
            rate1 *= population_sizes[reactant][-1]
            h1 *= population_sizes[reactant][-1]
            
        rate2 = k2
        h2 = 1
        for reactant in R2.reactants:
            rate2 *= population_sizes[reactant][-1]
            h2 *= population_sizes[reactant][-1] #h is the total number of reactant molecule combinations
            
        rate_total = rate1 + rate2
        #print(rate1, rate2)
        if rate_total == 0.0:
            break
        prob_R1 = rate1/rate_total #c1
        prob_R2 = rate2/rate_total #c2
        
        a1 = h1 * prob_R1
        a2 = h2 * prob_R2
        a0 = a1 + a2
        
        rand1 = random.random()
        tau = (1.0 / a0) * math.log(1.0 / rand1)
        
        #TODO
        rand2 = random.random()
        #print(rand2, prob_R1, prob_R2)
        if rand2 < prob_R1:
            for specie, population in population_sizes.items():
                if specie in R1.reactants:
                    population.append(population[-1] - 1)
                elif specie in R1.products:
                    population.append(population[-1] + 1)
                else: # Specie involved in R2 only, repeat last value
                    population.append(population[-1])
        else:
            for specie, population in population_sizes.items():
                if specie in R2.reactants:
                    population.append(population[-1] - 1)
                elif specie in R2.products:
                    population.append(population[-1] + 1)
                else: # Specie involved in R1 only, repeat last value
                    population.append(population[-1])
                
        time += tau
        tpoints.append(time)
    #end while loop
    
    #print(tpoints)
    for specie, population in population_sizes.items():
        #print(specie, len(tpoints), len(population))
        #print(population)
        pylab.plot(tpoints, population, label=specie)
    pylab.legend(loc='upper right')
    
    
def Runge_Kutta(R1, R2, Temp):
    
    a = 0.0
    b = 10.0
    N = 1000
    h = (b-a)/N
    
    def f(R1, R2, population):
        '''
        population: { 'specie1' : last_pop }
        return { 'specie1' : diff }
        '''
        r1 = R1.rate_constant(Temp)
        r2 = R2.rate_constant(Temp)
        
        pop_diff = defaultdict(float)
        for specie in R1.reactants:
            term = -r1
            for r in R1.reactants:
                term *= population[r]
            pop_diff[specie] += term
                
        for specie in R1.products:
            term = +r1
            for r in R1.reactants:
                term *= population[r]
            pop_diff[specie] += term

        for specie in R2.reactants:
            term = -r2
            for r in R2.reactants:
                term *= population[r]
            pop_diff[specie] += term
                
        for specie in R2.products:
            term = +r2
            for r in R2.reactants:
                term *= population[r]
            pop_diff[specie] += term
        
        return pop_diff
    
    
    def vector_add(pop_dict1, pop_dict2):
        result = dict()
        for specie, population in pop_dict1.items():
            result[specie] = pop_dict1[specie] + pop_dict2[specie]
        return result
    
    def vector_mul(pop_dict1, pop_dict2):
        result = dict()
        for specie, population in pop_dict1.items():
            result[specie] = pop_dict1[specie] * pop_dict2[specie]
        return result
        
    def vector_add_scalar(pop_dict, num):
        result = dict()
        for specie, population in pop_dict.items():
            result[specie] = pop_dict[specie] + num
        return result
        
    def vector_mul_scalar(pop_dict, num):
        result = dict()
        for specie, population in pop_dict.items():
            result[specie] = pop_dict[specie] * num
        return result
                    
    initial_population = 200    
    
    #Set initial populations size
    population_sizes = dict()
    last_population_sizes = dict()
    for reactant in R1.reactants:
        population_sizes[reactant] = [ initial_population ]
        last_population_sizes[reactant] = initial_population
    for reactant in R2.reactants:
        population_sizes[reactant] = [ initial_population ]
        last_population_sizes[reactant] = initial_population
    
    for product in R1.products:
        population_sizes[product] = [ 0 ]
        last_population_sizes[product] = 0
    for product in R2.products:
        population_sizes[product] = [ 0 ]
        last_population_sizes[product] = 0
        
    time = a
    tpoints = [time]
    while time < b:
        
        k1 = vector_mul_scalar(f(R1, R2, last_population_sizes), h)
        
        r_plus_05_k1 = vector_add(last_population_sizes, vector_mul_scalar(k1, 0.5))
        k2 = vector_mul_scalar(f(R1, R2, r_plus_05_k1), h)
        
        r_plus_05_k2 = vector_add(last_population_sizes, vector_mul_scalar(k2, 0.5))
        k3 = vector_mul_scalar(f(R1, R2, r_plus_05_k2), h)
        
        r_plus_k3 = vector_add(last_population_sizes, k3)
        k4 = vector_mul_scalar(f(R1, R2, r_plus_k3), h)
        
        change = vector_mul_scalar( vector_add( vector_add(k1, vector_mul_scalar(k2,2)), vector_add(vector_mul_scalar(k3,2), k4) ), 1.0/6.0)
        last_population_sizes = vector_add(last_population_sizes, change)
                
        #Append population size
        for specie, population in last_population_sizes.items():
            population_sizes[specie].append(population)

        
        time += h
        tpoints.append(time)    
        pass        
        '''C2Cl4points.append(r[0])
        Clpoints.append(r[1])
        C2Cl5points.append(r[2])
        
        k1 = h*f(r,t)
        k2 = h*f(r+0.5*k1,t+0.5*h)
        k3 = h*f(r+0.5*k2,t+0.5*h)
        k4 = h*f(r+k3,t+h)
        r += (k1+2*k2+2*k3+k4)/6
        '''
    
    for specie, population in population_sizes.items():
        #print(specie, len(tpoints), len(population))
        #print(population)
        pylab.plot(tpoints, population, label=specie)
        pylab.legend(loc='upper right')
    
    
    
def main():
    read_reaction_data('Scraper/reactions.csv')
    read_user_input()
    pass
    
if __name__ == "__main__":
    main()

    
