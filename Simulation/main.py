import csv

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
        

def read_reaction_data(filename):
    #Reaction	Order	Rall	R1	R2	Pall	P1	P2	P3	P4	Tmin	Tmax	A	n	Ea    
    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        global gReactions        
        for row in reader:
            r = Reaction(row)
            gReactions.append(r)            
            #print(row['Reaction'], row['Ea'])
        #end for loop
    #end with
    print('Loaded %d reactions' % len(gReactions))

def read_user_input():
    global gReactions
    #r1 = input('Enter Reactant 1: ').strip()
    r1 = 'ClO4'
    result_set = list( filter(lambda x: x.order == 2 and r1 in x.reactants, gReactions) )
    if not result_set:
        print('No reactions contain %s as a reactant' % r1)
        return
    print('%d 2nd order reactions contain reactant: %s' % (len(result_set), r1))
    #TODO: Print
    
    #r2 = input('Enter Reactant 2: ').strip()
    r2 = 'Cl'
    
    result_set = list( filter(lambda x: x.order == 2 and r2 in x.reactants, result_set) )
    if not result_set:
        print('Could not find 2nd order reaction with reactants: %s + %s' % (r1, r2))
        return
    elif len(result_set) > 1:
        print('Warning: More than one reaction found' )
    
    r_obj = result_set[0]
    print('Found 2nd order reactions: %s' % (r_obj.reaction))
    if r_obj.Tmin != r_obj.Tmax:
        Temp = int( input('Choose a temperature in the range %d - %d: ' % (r_obj.Tmin, r_obj.Tmax)) )
    else:
        Temp = r_obj.Tmin
    
    #Find 1st order equations involving the products
    product_reactions = dict()    
    for p in r_obj.products:
        result_set = list( filter(lambda x: x.order == 1 and r1 in x.products, gReactions) )
        if not result_set:
            continue
        if len(result_set) > 1:
            print('Warning: More than one reaction found' )
        product_reactions[p] = result_set[0]
        print('Found 1nd order reaction involving product %s : %s' % (p, result_set[0].reaction))
    

def main():
    read_reaction_data('Scraper/reactions.csv')
    read_user_input()
    pass
    
if __name__ == "__main__":
    main()

    
