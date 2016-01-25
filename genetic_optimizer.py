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

    def __init__(self,  sequence, penalty, crossover_rate, 
                 point_mutation_rate, chromosome = None):
        
        self.crossover_rate = crossover_rate
        self.point_mutation_rate = point_mutation_rate
        self.penalty = penalty
        self.sequence = sequence 
        self.chromosome = chromosome or self.create_chromosome()
        self.region_indices = None
        self.score = 0
        
        self.calculate_score()
    
    def create_chromosome(self):

        m = len(self.sequence)
        
        entry = random.randrange(0,m)
        r_width = random.randrange(0, m)
        l_width = random.randrange(0, m)
        contig_pass = random.random()

        return [entry, r_width, l_width, contig_pass]

    def point_mutation(self, gene):
        """Alters the value of a single gene by a small amount
        Larger changes are less likely than smaller ones"""
        
        def delta_probabilities():
            out = []
            for i in range(1,6):
                out.extend([i for x in range(100 / i) * 100])
            out.extend([-x for x in out])
            return out

        deltas = delta_probabilities()
        delta_float =  map(lambda x: x / 100., deltas)
        
        if type(gene) == int:
            self.chromosome[gene] += random.choice(deltas)

        else:
            self.chromosome[gene] += random.choice(delta_float)

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

    def fix_chromosome(self):

        # get left position
        if self.chromosome[1] < 0:
            self.chromosome[1] = 0 

        # fix contig extension probabilities
        if self.chromosome[3] > 1.0:
            self.chromosome[3] = 1.0

        if self.chromosome[3] < 0.0:
            self.chromosome[3] = 0.0

    def get_region(self):
    
        def check_contigs(self, region):

            def get_edge(region, end, step):

                halfway = len(region) / 2
                
                for kmer in range(halfway, end, step):
                    edge = kmer
                    broken = False
                    if region[kmer].contigs > 0:
                        for contig in range(region[kmer].contigs):
                            if random.random() > self.chromosome[3]:
                                edge = kmer + -step
                                broken = True
                                break
                        if broken:
                            break
               
                if step == 1:
                    edge += 1 
                return edge

            start = get_edge(region, -1, -1)
            stop  = get_edge(region, len(region) - 1, 1)
            
            return start, stop
        
        # get left edge of inital region
        left = self.chromosome[0] - self.chromosome[1]
        if left < 0: left = 0

        # out of bounds doesn't matter - it's treated as len(sequence) + 1
        right = self.chromosome[0] + self.chromosome[2]

        region = self.sequence[left : right + 1]
        length = len(region)

        start, stop  = check_contigs(self, region)
       
        # returns selected region and its location within overall sequence
        return region[start : stop], (left + start), (right - (length - stop))

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
        
        return cmp(self.score, other.score)

class Population(object):

    def __init__(self, sequence, size, generations, crossover_rate,
                 point_mutation_rate, penalty, elite):
        
        self.sequence = sequence
        self.crossover_rate = crossover_rate
        self.point_mutation_rate = point_mutation_rate
        self.elite = elite
        self.generations = generations

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
        
        self.population = next_population

    def step(self):

        self.population.sort()
        self.selection()
        self.evolve()
        self.generations += 1
        
        # diagnostic
        if self.generations % 100 == 0:  
            print [x.chromosome for x in self.population]
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
        plt.xlim(0, max(self.data["kmer"]))
        plt.show()
