import copy
import pandas as pd
from optimizers.base_optimizer import Optimizer

class GeneticOptimizer(Optimizer):
    def __init__(self, params):
        self.tourn_size = params.tourn_size
        self.scorer_type = params.scorer_type

        self.selection_type = params.selection_type.lower()
        if self.selection_type not in ['tournament', 'roulette']:
            raise ValueError(f"Unknown optima type {self.selection_type}. Available options are: {['tournament', 'roulette']}")

        self.optima_type = params.optima_type.lower()
        if self.optima_type not in ['minima', 'maxima']:
            raise ValueError(f"Unknown optima type {self.optima_type}. Available options are: {['minima', 'maxima']}")

    def tournament_selection(self, selection_pool):
        selection_pool[["fitness"]] = selection_pool[["fitness"]].apply(pd.to_numeric)

        if self.optima_type == "minima":
            return selection_pool.fitness.idxmin()
        elif self.optima_type == "maxima":
            return selection_pool.fitness.idxmax()

    def select_non_elite(self, population, size):
        self.population = copy.deepcopy(population)
        if self.scorer_type == "LogPTestCase":
            selected_population = pd.DataFrame(columns=['compound_id', 'smiles', 'generation', 'chromosome', 'fitness'])
        elif self.scorer_type == "SeaLikeTanimoto":
            selected_population = pd.DataFrame(columns=['compound_id', 'smiles', 'generation', 'chromosome', 'fitness', 'avg_sea_like_TC'])

        while len(selected_population) < size:
            #Setup the pool of individuals for the tournament selection to be a random sampling of the population
            #Without replacement means that the same individual will not appear in the sampling more than once
            if self.tourn_size > len(self.population):
                selection_pool = self.population.sample(len(self.population), replace=False) #TODO: Can create an alternative WITH replacement
            else:
                selection_pool = self.population.sample(self.tourn_size, replace=False) #TODO: Can create an alternative WITH replacement
           
            bool_fit = pd.notnull(selection_pool['fitness'])
            selection_pool = selection_pool[bool_fit]

            if len(selection_pool) == 0 and self.tourn_size > len(self.population):
                selection_pool = self.population.sample(len(self.population), replace=False)

                while not selection_pool['fitness'].notnull().any():
                    selection_pool = self.population.sample(len(self.population), replace=False)

                bool_fit = pd.notnull(selection_pool['fitness'])
                selection_pool = selection_pool[bool_fit]

            elif len(selection_pool) == 0:
                selection_pool = self.population.sample(self.tourn_size, replace=False)

                while not selection_pool['fitness'].notnull().any():
                    selection_pool = self.population.sample(self.tourn_size, replace=False)

                bool_fit = pd.notnull(selection_pool['fitness'])
                selection_pool = selection_pool[bool_fit]


            if self.selection_type == "tournament":
                selected_individual = self.tournament_selection(selection_pool)
            elif self.selection_type == "roulette":
                selected_individual = self.roulette_selection(selection_pool)

            #Now we can add the selected individual to the selected_population and remove it from the self.population
            # so that we do not have repeated individuals in our selected population. This can lead to false convergence
            #selected_population = pd.concat([selected_population,selection_pool.loc[selected_individual]])
            
            selected_population.loc[len(selected_population.index)] = selection_pool.loc[selected_individual]
            #TODO: Do we need to allow repeated individuals as the population converges? 
            self.population.drop([selected_individual], inplace=True) 

        #Reset the self.population to contain the selected individuals that will be used for creation of next generation
        self.population = selected_population.reset_index(drop=True)
        print("Population_size after selection ", len(self.population))

        return selected_population

    def select_elite_pop(self, population, size):
        population = copy.deepcopy(population)
        
        if self.optima_type == "minima":
            sort_order = True
        elif self.optima_type == "maxima":
            sort_order = False

        selected_population = population.sort_values(by=['fitness'], ascending=sort_order)
        selected_population = selected_population.head(size)

        return selected_population
