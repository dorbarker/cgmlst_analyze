import random
import pandas
from numpy import arange
from matplotlib import pyplot as plt

### Sequence ###

# list of kmers
    # each cell
    # shannon
    # contig count
        # penalty * contigs

class Kmer(object):
    
    def __init__(self, shannon, contigs):

        self.shannon = shannon
        self.contigs = contigs
        
    def __repr__(self):
        return ','.join([str(self.shannon), str(self.contigs)])

### Chromosome values ###

# 0 entry point 
# 1 width left
    # maybe probability instead of absolute value
# 2 width right
    # maybe probability instead of absolute value
# 3 probability to ignore contig

class Individual(object):
    
    alleles = (0,1)
    functions = {'entry': 0, 'left': 1, 'right': 2, 'contig_pass': 3}

    def __init__(self,  sequence, penalty, crossover_rate, 
                 point_mutation_rate, chromosome = None):
        
        self.motif = None 
        self.region_indices = None
        
        self.crossover_rate = crossover_rate
        self.point_mutation_rate = point_mutation_rate
        self.penalty = penalty
        self.sequence = sequence 
        self.chromosome = self.create_chromosome(chromosome)
        self.score = 0
        
        self.calculate_score()
    
    def create_chromosome(self, chromosome):

        seq_len = len(self.sequence)
        bin_seq_len = bin(seq_len).lstrip('0b')
        
        motif = len(bin_seq_len) 
        self.motif = motif
        
        # entry point, left and right widths, and contig pass 
        # encoded in regions of length motif
        if chromosome == None:
            initialize = lambda x: [random.choice(self.alleles) for i in range(x)]   
            l = [initialize(motif) for x in range(4)]
            out = [x for y in l for x in y]
        else:
            out = chromosome
        
        return out 

    def point_mutation(self, gene):
        """Randomly mutates a gene.
        May bit flip or return the starting allele
        """

        self.chromosome[gene] = random.choice(self.alleles)

    def crossover(self, other):

        point = random.randrange(len(self.chromosome))
        self_fragment = self.chromosome[point :]
        other_fragment = other.chromosome[point :]
       
        
        child1_chromosome = self.chromosome[: point] + other_fragment
        
        child1 = Individual(self.sequence, self.penalty, self.crossover_rate,
                        self.point_mutation_rate, chromosome = child1_chromosome)
        
        child2_chromosome = other.chromosome[: point] + self_fragment
        
        child2 = Individual(self.sequence, self.penalty, self.crossover_rate,
                        self.point_mutation_rate, chromosome = child2_chromosome)

        return [child1, child2]
    
    def decimalize_functions(self):

        sub = lambda thing: self.chromosome_function_indices(thing)
        subseq = lambda thing: self.chromosome[sub(thing)[0] : sub(thing)[1]]
        
        return {thing: self.decimalize_binary(subseq(self.functions[thing])) for thing in self.functions}

    def decimalize_binary(self, binary_list):
        """Converts binary subregion of chromosome to decimal"""    
        
        return int(''.join(map(str, binary_list)), 2)

    def chromosome_function_indices(self, thing):
        
        chunk = thing * self.motif
        
        return chunk, chunk + self.motif
        
    def get_region(self):
    
        def check_contigs(self, region, contig_pass):
            
            def get_edge(region, end, step, contig_pass):
                
                halfway = len(region) / 2
                edge = None
                
                for kmer in range(halfway, end, step):
                    edge = kmer
                    broken = False
                    if region[kmer].contigs > 0:
                        for contig in range(region[kmer].contigs):
                            if random.random() < contig_pass:
                                edge = kmer + -step
                                broken = True
                                break
                        if broken:
                            break
               
                if step == 1 and edge != None:
                    edge += 1 
                return edge

            start = get_edge(region, -1, -1, contig_pass)
            stop  = get_edge(region, len(region) - 1, 1, contig_pass)
            
            return start, stop
        
        decimalized = self.decimalize_functions()

        # get left edge of inital region
        if decimalized['entry'] - decimalized['left'] < 0:
            left = 0
        else:
            left = decimalized['entry'] - decimalized['left']

        # out of bounds doesn't matter - it's treated as len(sequence) + 1
        right = decimalized['entry'] + decimalized['right']

        # get probability of extending past a contig
        max_bin = float(self.decimalize_binary(1 for x in range(self.motif)))
        contig_pass = decimalized['contig_pass'] / max_bin 
        
        region = self.sequence[left : right + 1]
        length = len(region)
        
        try: 
            start, stop  = check_contigs(self, region, contig_pass)
            out = region[start : stop], (left + start), (right - (length - stop))
        except IndexError:
            out = [], None, None
        except TypeError:
            out = [], None, None

        # returns selected region and its location within overall sequence
        
        return out

    def calculate_score(self):
        
        region, l_ind, r_ind = self.get_region()
        
        self.score += sum(kmer.shannon for kmer in region)
        self.score -= sum(kmer.contigs * self.penalty for kmer in region)
        self.region_indices = (l_ind, r_ind)

    def clone(self):

        offspring = self.__class__(self.sequence,
                                   self.penalty,
                                   self.crossover_rate,
                                   self.point_mutation_rate,
                                   self.chromosome[:])
        return offspring

    def __repr__(self):

        return ','.join([self.score] + map(str, self.chromosome))

    def __cmp__(self, other):
        
        #return cmp(self.score, other.score)
        return cmp(other.score, self.score)
class Population(object):

    def __init__(self, sequence, size, crossover_rate,
                 point_mutation_rate, penalty, elite):
        
        self.size = size 
        self.sequence = sequence
        self.crossover_rate = crossover_rate
        self.point_mutation_rate = point_mutation_rate
        self.elite = elite
        self.generations = 0  

        self.population = self.generate_population(sequence,
                                                   size, 
                                                   crossover_rate,
                                                   point_mutation_rate,
                                                   penalty)

    def generate_population(self, sequence, size, crossover_rate,
                            point_mutation_rate, penalty) :

        population =  [Individual(sequence, penalty, crossover_rate, point_mutation_rate)
                for x in range(size)]
        
        return population
   
    def selection(self):

        # transfers the best performers to the next population
        next_population = [x.clone() for x in self.population[:self.elite]]
        
        while len(next_population) < len(self.population):
            parent1 = self.tournament_select()
            if random.random() < self.crossover_rate:
                parent2 = self.tournament_select()
                offspring = parent1.crossover(parent2)
            else:
                offspring = [parent1.clone()]
            
            next_population.extend(offspring)
       
        next_population = next_population[:self.size] # hack to force correct pop size
        self.population = next_population

    def step(self):

        self.population.sort(reverse = True)
        self.selection()
        self.evolve()
        self.generations += 1
        
        # diagnostic
        if self.generations % 1000 == 0:  
            
            reporter = StepReporter(self.sequence, self.population)
            reporter.plot_agents()

            
    def evolve(self):
        
        for individual in self.population[self.elite:]:
            for i, gene in enumerate(individual.chromosome):
                if random.random() < self.point_mutation_rate:
                    individual.point_mutation(i)

    def tournament_select(self, size = 10, p = 0.90):
        """Selects survivors for next step of the simulation
        Based on https://en.wikipedia.org/wiki/Tournament_selection
        """

        subset = [random.choice(self.population) for x in range(size)]
        subset.sort()
        
        if random.random() < p:
            return subset[0]

        while True: 
            for i in range(1, len(subset) + 1):
                if random.random() < (p * ((1.0 - p) ** i)):
                    return subset[i]

class StepReporter(object):

    def __init__(self, sequence, population):
        
        self.population = population
        self.sequence = sequence
        self.regions = self.get_regions()
        #self.dist = self.distribution()
        self.data = self.convert_to_data_frame()

    def get_regions(self):
        """Returns start and stop positions of all agents in the population"""
        
        return [individual.region_indices for individual in sorted(self.population)]

    def distribution(self):
        """How many agents include each kmer"""

        dist = []
        
        for start, stop in self.regions:
            dist.extend(range(start, stop + 1))
        
        return dist
    
    def convert_to_data_frame(self):
        
        starts, stops = zip(*self.regions)
        data = zip( range(len(self.sequence)),
                    [x.shannon for x in self.sequence],
                    [x.contigs for x in self.sequence],
                    starts,
                    stops)

        df = pandas.DataFrame(data = data,
                              columns = ["kmer", "shannon", "contigs",
                                         "starts", "stops"])
        
        return df

    def plot_agents(self):
       
        # float() and nested max() to dodge type errors
        maxvalue = float(max(max(self.data["contigs"]), max(self.data["shannon"])))
        slices = maxvalue / len(self.population)
        yrange = arange(1.0, maxvalue + 1.0, slices)

        plt.plot(self.data["kmer"], self.data["shannon"], color = 'blue')
        plt.fill_between(self.data["kmer"], self.data["shannon"], 0, alpha = 0.3) 

        plt.bar(self.data["kmer"], self.data["contigs"], color = 'red', alpha = 0.3)
        
        plt.hlines(y = yrange, xmin = self.data["starts"], xmax = self.data["stops"])


        plt.xlim(0, max(self.data["kmer"]))
        plt.show()
