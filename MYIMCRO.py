# %%
import random
import math
import copy
from typing import List,Tuple

# %% [markdown]
# # MY IMCRO

# %%
# Initialization

# %% [markdown]
# # Molecule Class
# 
# 

# %%
class MoleCule:
    def __init__(self, structure: List[float] ,supersequence: List[int]) -> None:
        '''
            PE : Potential Energy (This is the number we focus on during the algorithm)

            structure : Molecule Structe store as 1d array with 2 elements ge
                generated randomly (For example: [0.3860863116352774 0.4017794415965995])

        '''
        self.pe = 0
        self.ke = 0
        self.opt = 9999999
        self.num_of_hits = 0
        self.structure = copy.deepcopy(structure)
        self.supersequence = copy.deepcopy(supersequence)

    def update(self) -> None:
        """
        This is called whenever a Operator is performed
        If this molecule has a lower energy, reset num_of_hits.
        """
        if self.pe < self.opt:
            self.opt = self.pe
            self.num_of_hits = 0
    def __str__(self) -> str:
        reVal = 'Structure'
        reVal +='[ '
        for i in self.structure:
            reVal+=str(i)+" "
        reVal +=']  Supersequence[ '
        for i in self.supersequence:
            reVal+=str(i)+" "
        reVal +=' ]'
        return reVal

# %% [markdown]
# # Initialization

# %%
def insert_symbol(src_string,inserted_string,pos):
    return ''.join(src_string[:pos] + inserted_string + src_string[pos:])
# Given set of strings and population size for SCS problem
def supersequence_generate(set_of_strings):

    '''
        Make a copy of the set_of_strings parameter for maintaining the original 
        set
    '''
    copied_set_of_strings = copy.deepcopy(set_of_strings)
    supersequence = ''.join(copied_set_of_strings.pop(random.randint(0,len(set_of_strings)-1)))

    for i in range(len(copied_set_of_strings)):
        # print("i = ",i)
        counter = 0
        for j in copied_set_of_strings[i]:
            inserted_pos = random.randint(counter,len(supersequence))
            # print("j and counter and supersequence length and inserted index",j," ",counter," ",len(supersequence)," ",inserted_pos)
            if inserted_pos == len(supersequence) or j != supersequence[inserted_pos]:
                supersequence=insert_symbol(supersequence,j,inserted_pos)
            counter = inserted_pos +1
            # print(supersequence)

    return supersequence


def population_generation(pop_size,set_of_strings):
    l=[]
    for i in range(pop_size):
        l.append(supersequence_generate(set_of_strings))
    return l



# %%
def encoding_population(initial_population):
    '''
        Encoding base on the class of the data
    '''
    dict = {
        'a' : 0,
        'c' : 1,
        'g' : 2,
        't' : 3
    }
    l=[]
    for i in initial_population:
        k=''
        for j in i:
            k+=str(dict[j])
        l.append(k)
    return l

def convert_sequences(sequences):
    nucleotide_dict = {
        'a': 0,
        'c': 1,
        'g': 2,
        't': 3
    }
    
    converted_sequences = []
    for sequence in sequences:
        converted_sequence = [nucleotide_dict[nucleotide] for nucleotide in sequence]
        converted_sequences.append(converted_sequence)

    return converted_sequences


# %%
def createMolecule(string):
    l=[]
    for i in string:
        l.append(int(i))
    return MoleCule(structure=[random.random(), random.random()],supersequence=l)

# %%
def is_subsequence(sub, main):
    it = iter(main)
    return all(any(item == x for x in it) for item in sub)

def is_supersequence(seq):
    subsequences = [[0, 1, 2], [1, 0, 3], [2, 3, 3], [3, 2, 1]]
    for i in subsequences:
        if not is_subsequence(i,seq):
            return False
    return True

# Example usage
sequence = [ 0, 2, 1, 0, 2, 0, 1, 2, 3, 0, 0, 2]
result = is_supersequence(sequence)
print(f"The sequence {sequence} is a supersequence: {result}")


# %%
def initialization(pop_size,set_of_strings) ->List[MoleCule]:
    initial_population =  population_generation(pop_size,set_of_strings)
    print("Population Generation:")
    print(initial_population)
    encoded_population = encoding_population(initial_population)
    print("Encodeding:")
    print(encoded_population)
    l=[]
    for i in encoded_population:
        l.append(createMolecule(i))
    return l

molecules=initialization(20,['acg', 'cat', 'gtt','tgc'])
# 0 1 2 | 1 0 3 | 2 3 3 | 3 2 1
for i in molecules:
    print(i)

# %%
class Operators():
    
    # OnWall Ineffective Colision
    def OnWall (self,molecule : List[int]) -> List[int]:
        print("On Wall: ",molecule)
        '''
        Objective: 
        -  This focuses on instances of molecules colliding with container walls, resulting in structural transformations.
        - A "one-difference operator" is used to make a single change in the molecule's composition to achieve this.

        Input:
        - molecule (list): the input molecule and it represent by list

        Output:
        - The method returns a new list. 
        '''
        #Initial
        m = molecule[:] 
        v_low = 0
        v_upper = 3   
        i = random.randint(0, len(molecule)-1)
        j = random.randint(v_low, v_upper)

        if (molecule[i] + j <= v_upper):    
            m[i] = molecule[i] + j
        else:
            if(molecule[i] - j < 0):
                m[i] = 0
            else:
                m[i] = molecule[i] - j
        
        #Test     
        # print(i, j)
        return m

    # Intermolecular Ineffective Colision
    def Intermolecular(self,molecule1 : List[int], molecule2 : List[int]) -> Tuple[List[int],List[int]]:
        print("Intermolecule: ",molecule1," Length: ",len(molecule1)," ",molecule2," Length: ",len(molecule2))
        '''
        Objective: 
        -  The purpose is to introduce significant changes to enhance local search capabilities and prevent getting stuck in local optimization by promoting diversity.
        - A crossover operator is used in genetic or evolutionary algorithms for optimization. It selects two molecules from the population and uses a two-step mechanism to generate two new solutions.
        - It is a two step process: the first step is to crossover between two molecules, and the second step is to crossover inside the molecule itself

        Input:
        - molecule (list): the input molecule and it represent by list

        Output:
        - The method returns a tuple (m1, m2), where m1 and m2 are the two molecules and m1, m2 are also list. 
        '''
        #Initial 
        # Length of molecule
        # print(molecule1,molecule2)
        length1 = len(molecule1)
        length2 = len(molecule2)
    
        #Two new molecule in first crossover
        # Copy 2 molecules
        m1 = molecule1.copy()
        m2 = molecule2.copy()
        #Two new molecule in second crossover
        m11 = [0] * length1
        m22 = [0] * length2
        
        limit = min(length2, length1) 
        
        #Random numbers x1, x2 generation
        x1 = random.randint(0, limit-2)
        x2 = random.randint(x1+1, limit-1)
    
        # Randormly choose form molecule1 or molecule2
        # Crossover 1
        for i in range(0,limit):
            if (i<x1 or i>x2):  #if odd segments
                m1[i] = molecule1[i]
                m2[i] = molecule2[i]
            elif (i>=x1 and x2>=i): # if even segment
                m1[i] = molecule2[i]
                m2[i] = molecule1[i]
        
        # Crossover 2
        # Crossover 2 for molecule m1
        for j in range(0,length1):
            if (j < x1):  #if odd segments
                m11[length1 - x1 + j] = m1[j]
                
            elif (j>=x1 and x2>=j): # if even segment
                m11[length1 - x1 - x2 + j - 1] = m1[j]
            else:
                m11[j - x2-1] = m1[j]
        
        # Crossover 2 for molecule m2
        for j in range(0,length2):
            if (j < x1):  #if odd segments
                m22[length2 - x1 + j] = m2[j]
            elif (j>=x1 and x2>=j): # if even segment
                m22[length2 - x1 - x2 + j- 1] = m2[j]
            else:
                m22[j - x2-1] = m2[j]
        
        #Test
        
        # print(limit)
        # print(x1, x2)
        # print(m1)
        # print(m2)
        print(" After Intermolecule: ",m11," ",m22)
        return m11,m22
    
    # Decomposition
    def Decomposition(self, molecule : List[int]) -> Tuple[List[int],List[int]]:
        print("Decomposition")
        '''
        Objective: 
        - The decomposition involves randomly selecting two numbers 'a' and 'b', and then splitting the input molecule into two new molecules, 'm1' and 'm2', based on the selected numbers. 
        - The negative number −a is used for shifting to the left a steps. 
        - The positive number j is used for shifting to the right j steps.

        Input:
        - molecule (list): the input molecule and it represent by list

        Output:
        - The method returns a tuple (m1, m2), where m1 and m2 are the two molecules and m1, m2 are also list. 
        '''
        # Step 1: Select two numbers a and b randomly
        a = random.randint(-len(molecule), 0)
        b = random.randint(0, len(molecule))
        
        # Initialize m1 and m2
        m1 = [0] * len(molecule)
        m2 = [0] * len(molecule)

        # Step 2: Decomposition of molecule into m1
        for i in range(len(molecule)):
            if i + 1 <= -a:
                m1[len(molecule) - (-a) + i] = molecule[i]
            else:
                m1[i - (-a)] = molecule[i]

        # Step 3: Decomposition of molecule into m2
        for j in range(len(molecule)):
            if j + 1 <= len(molecule) - b:
                m2[j + b] = molecule[j]
            else:
                m2[j - len(molecule) + b] = molecule[j]
                
        #Test
        # print(molecule)
        # print(a, b)
        
        return m1, m2

    # Synthesis
    def Synthesis(self, molecule1 : str, molecule2)-> List[int]:
        print("Systhesis")
        """
        Objective:
        - Generates a new list by combining two input lists in a way that preserves the frequency of the symbols used in each input list.

        Input:
        - molecule1 (list): The first input list.
        - molecule2 (list): The second input list.

        Output:
        - The method returns a new list. 
        """
        # Generate dictionaries for the frequencies of the symbols used in each input list.
        array1 = {}
        for symbol in molecule1:
            if symbol not in array1:
                array1[symbol] = 0
            array1[symbol] += 1

        array2 = {}
        for symbol in molecule2:
            if symbol not in array2:
                array2[symbol] = 0
            array2[symbol] += 1

        # Initialize the output list.
        length1 = len(molecule1)
        length2 = len(molecule2)
        limit = min(length1, length2)
        
        if(length1 < length2):
            m_prime = molecule2.copy()
        else:
            m_prime = molecule1.copy()

        # Iterate over the symbols in the first input list.
        for i in range(limit):
            symbol1 = molecule1[i]
            symbol2 = molecule2[i]

            frequency_in_array1 = array1.get(symbol1, 0)
            frequency_in_array2 = array2.get(symbol2, 0)

            if frequency_in_array1 >= frequency_in_array2:
                m_prime[i] = symbol1
            # Otherwise, add the symbol from the second input list to the output list.
            else:
                m_prime[i] = symbol2
        #test
        
        # print(molecule1)
        # print(molecule2)
        # print(array1)
        # print(array2)
        
        return m_prime

# %% [markdown]
# ![alt text](image.png)

# %%
# Iteration
class IMCRO():
    optimal : MoleCule = None
    moleColl = 0.2
    alpha = random.randint(10, 100)
    beta = random.randint(10, 100)
    KELossRate = 0.6
    init_ke = 100
    buffer = 0
    pop : List[MoleCule]=[]
    ops : Operators = None
    checkSeqeunce : List[List[int]] = None
    '''
        pop : popsize the List Of Molecule Instance
        ops : for loa
        optimal : This attribute will hold the Molecule Instance with the lowest PE
        Which is also the output of this algorithm
    '''

    def __init__(self,pop : List[MoleCule],checkSeqeunce : List[List[int]]) -> None:
        '''
            Initialize the algorithm with this constructor 
            The parameter (Or the input) is the popsize which mean a List of Instance of Molecule Class
        '''
        self.pop=pop
        self.ops = Operators() # load in the operator for supersequence attribute of the molecule instance
        self.checkSeqeunce = checkSeqeunce

        for mol in self.pop:
            '''
                The fit_func function is self-defined based on the data
                In this example the PE formula is defined as math.sin(m.structure[0]) + math.cos(m.structure[1]) 
                Which i have deveried from this class
            '''
            mol.pe = self.fit_func(mol)
            mol.ke = self.init_ke
            mol.update()
            if self.optimal is None:
                self.optimal = copy.deepcopy(mol)
            elif mol.pe < self.optimal.pe:
                self.optimal = copy.deepcopy(mol)
    
    def run(self) -> None:
        '''
            The Algorithm starts here
            randomly pick which reaction uni(randomly take 1 molecule from popsize) 
            or inter (randomly take 2 molecule from popsize) 

            The number of iteration is abritrary 
        '''
        i=0
        while i!=1000:
            i+=1
            t = random.random()
            if t > self.moleColl or len(self.pop) < 2:
                self.uni_moleReact()
            else :
                self.inter_moleReact()
    
    def is_subsequence(sub, main) -> bool:
        it = iter(main)
        return all(any(item == x for x in it) for item in sub)

    def is_supersequence(self,seq : List[int]) -> bool:
        subsequences = self.checkSeqeunce
        for i in subsequences:
            if not is_subsequence(i,seq):
                return False
        return True
    
    def update(self, m : MoleCule) -> None:
        """
            If m is the current optimal solution, save it to the optimal. 
            This is called after an operator is executed
        """
        if m.pe < self.optimal.pe and is_supersequence(m.supersequence):
            self.optimal = copy.deepcopy(m)

    def uni_moleReact(self) -> None:
        '''
            Randomly pick a molecule from popsize
        '''
        m = random.choice(self.pop)
        if m.num_of_hits > self.alpha:
            self.decomposition(m)
        else :
            self.on_wall(m)

    def inter_moleReact(self) -> None:
        '''
            Randomly pick 2 molecule from popsize
        '''
        m1 , m2 = random.sample(self.pop,2)
        if m1.ke <= self.beta and m2.ke <= self.beta:
            self.synthesis(m1,m2)
        else :
            self.interaction(m1,m2)
    
    def decomposition(self,m : MoleCule) -> None:
        m.num_of_hits += 1

        # You should implement this function in your derived class
        new1, new2 = self.dec(m)
        new1.pe = self.fit_func(new1)
        new2.pe = self.fit_func(new2)
        tmp = m.pe + m.ke - new1.pe - new2.pe
        if tmp >= 0 or tmp + self.buffer >= 0:
            if tmp >= 0:
                q = random.random()
                new1.ke = tmp * q
                new2.ke = tmp * (1 - q)
            else:
                new1.ke = (tmp + self.buffer) * random.random() * random.random()
                new2.ke = (tmp + self.buffer - new1.ke) * random.random() * random.random()
                self.buffer = self.buffer + tmp - new1.ke - new2.ke
            new1.update()
            new2.update()
            self.pop.remove(m)
            self.pop.append(new1)
            self.pop.append(new2)
            self.update(new1)
            self.update(new2)

    def on_wall(self, m : MoleCule) -> None:
        m.num_of_hits += 1
        # You should implement this function in your derived class
        new = self.wall(m)
        new.pe = self.fit_func(new)
        tmp = m.pe + m.ke - new.pe
        if tmp >= 0:
            q = random.uniform(self.KELossRate, 1)
            new.ke = tmp * q
            new.update()
            self.buffer += tmp * (1 - q)
            self.pop.remove(m)
            self.pop.append(new)
            self.update(new)
            
    def interaction(self, m1 : MoleCule, m2 : MoleCule) -> None:
        m1.num_of_hits += 1
        m2.num_of_hits += 1

        # You should implement this function in your derived class
        new1, new2 = self.inter(m1, m2)
        new1.pe = self.fit_func(new1)
        new2.pe = self.fit_func(new2)
        tmp = m1.pe + m2.pe + m1.ke + m2.ke - new1.pe - new2.pe
        if tmp >= 0:
            q = random.random()
            new1.ke = tmp * q
            new2.ke = tmp * (1 - q)
            new1.update()
            new2.update()
            self.pop.remove(m1)
            self.pop.remove(m2)
            self.pop.append(new1)
            self.pop.append(new2)
            self.update(new1)
            self.update(new2)

    def synthesis(self, m1 : MoleCule, m2 : MoleCule) -> None:
        m1.num_of_hits += 1
        m2.num_of_hits += 1

        # You should implement this function in your derived class
        new = self.syn(m1, m2)
        new.pe = self.fit_func(new)
        tmp = m1.pe + m2.pe + m1.ke + m2.ke - new.pe
        if tmp >= 0:
            new.ke = tmp
            new.update()
            self.pop.remove(m1)
            self.pop.remove(m2)
            self.pop.append(new)
            self.update(new)

# %%
class MYIMCRO(IMCRO):
    def __init__(self,pop : List[MoleCule],checkSeqeunce : List[List[int]]) -> None:
        super().__init__(pop,checkSeqeunce)


    def fit_func(self, m : MoleCule):
            '''
            fit_func function for deciding PE base on the struct attribute
                     
            '''
            return math.sin(m.structure[0]) + math.cos(m.structure[1])
    
    def dec(self, m : MoleCule) -> List[MoleCule]:
        new1 = copy.deepcopy(m)
        new2 = copy.deepcopy(m)
        new1.supersequence,new2.supersequence = self.ops.Decomposition(m.supersequence)
        new1.structure[0] += random.random()
        new2.structure[1] += random.random()
        return [new1, new2]
    
    def wall(self, m : MoleCule) -> MoleCule:
        new = copy.deepcopy(m)
        new.structure[0], new.structure[1] = new.structure[1], new.structure[0]
        new.supersequence = self.ops.OnWall(m.supersequence)
        return new
    
    def inter(self, m1 : MoleCule, m2 : MoleCule) -> List[MoleCule]:

        new1 = copy.deepcopy(m1)
        new2 = copy.deepcopy(m2)
        new1.supersequence,new2.supersequence = self.ops.Intermolecular(m1.supersequence,m2.supersequence)
        new1.structure[0] = m2.structure[0]
        new1.structure[1] = m1.structure[1]
        new2.structure[0] = m1.structure[0]
        new2.structure[1] = m2.structure[1]
        return [new1, new2]
    
    def syn(self, m1 : MoleCule, m2 : MoleCule) -> MoleCule:
        new = copy.deepcopy(m1)
        new.supersequence = self.ops.Synthesis(m1.supersequence,m2.supersequence)
        if random.random() < 0.5:
            new.structure[0] = m2.structure[0]
        else:
            new.structure[1] = m2.structure[1]
        return new

# %%
initial = initialization(pop_size=20,set_of_strings=['acg', 'cat', 'gtt','tgc'])

myimcro = MYIMCRO(pop=initial,checkSeqeunce=convert_sequences(['acg', 'cat', 'gtt','tgc']))
myimcro.run()
print(myimcro.optimal)


