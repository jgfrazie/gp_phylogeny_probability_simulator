"""
File: phylogeny_simulator.py
Author: James Frazier
DlM: 11/08/2023
DESC: Will simuate evolution for individuals in some evolutionary algorithm and use tournament selection as the method of selecting
      new individuals to reproduce. Each individual is a list where each element is a tuple of generation the error comes from and
      the error for a specified case.
"""

import random

### Size of tournaments in tournament selection
TOURNAMENT_SIZE = 7

def generate_initial_error_vector(num_cases):
    """Generate one individual for the first generation with num_cases"""
    return [(0, random.randint(0, 9999)) for _ in range(num_cases)]


def generate_initial_population(p, num_cases):
    """Generate p individuals for the first generation with num_cases"""
    return [generate_initial_error_vector(num_cases) for _ in range(p)]


def tournament_selection(population):
    """Performs tournament selection on a population

    Args:
        population ([[(int, int)]]): A population of individuals

    Returns:
        [(int, int)]: The individual with the best error from the tournament selection
    """
    candidates = random.sample(population, TOURNAMENT_SIZE)
    ### lambda function is to gather total error for an individual
    return max(candidates, key=lambda individual: sum([case_states[1] for case_states in individual]))


def downsample(population, rate=0.1):
    """Returns the case ids of what cases are to be evaluated for a generation

    Args:
        population ([[(int, int)]]): A population of individuals
        rate (float): The percent rate of downsampling for case evaluation

    Returns:
        [int]: Case ids to "evaluate" for a generation
    """
    num_downsampled_cases = int(rate * len(population[0]))
    cases_to_replace = random.sample([case_index for case_index in range(len(population[0]))], num_downsampled_cases)
    cases_to_replace.sort() ### Sorted to hold assertions for other functions
    return cases_to_replace


def get_new_generation_age(population):
    """Find the next generation iteration count based on the current population"""
    return max([case_specs[0] for case_specs in population[0]]) + 1


def select_and_vary(population, cases_to_replace):
    """Using tournament selection, will select and vary individuals in a population and replace specified cases
       with new errors. This function simulates mutation genetic operators only

    Args:
        population ([[(int, int)]]): A population of individuals
        cases_to_replace ([int]): Ids of cases to re-evaluate for the generation

    Returns:
        [[(int, int)]]: A new population of individuals with the cases specified replaced with new errors and the current
                        generation's iteration number
    """
    selected_to_vary = [tournament_selection(population) for _ in range(len(population))]
    new_population = []
    new_generation_age = get_new_generation_age(population)
    for parent in selected_to_vary:
        child = []
        for case_id in range(len(parent)):
            ### If a case is one to be replaced, we do so. Otherwise, we simply inherit from the parent
            child.append((new_generation_age, random.randint(0, 9999)) if case_id in cases_to_replace else parent[case_id])
        new_population.append(child)

    return new_population


def simulate_one_generation(population, downsample_rate=0.1):
    """Simualtes a single generation of a EA"""
    return select_and_vary(population, downsample(population, downsample_rate))


def simualte_evolution(num_generations, pop_size, num_cases, downsample_rate=0.1):
    """Simualtes the evolutionary cycle for some problem and returns all populations as evolution continued on

    Args:
        num_generations (int): Number of generations to run
        pop_size (int): The size of the popuations for each generation
        num_cases (int): The number of cases each individual is "evaluated" on
        downsample_rate (float): The rate for downsampling cases for evaluation

    Returns:
        [[[(int, int)]]]: A list of populations where the element at index 'i' is the population after select_and_vary
                          is applied which is for generation 'i'
    """
    populations_of_generations = [generate_initial_population(pop_size, num_cases)]
    for _ in range(num_generations):
        populations_of_generations.append(simulate_one_generation(populations_of_generations[len(populations_of_generations) - 1], downsample_rate))

    return populations_of_generations


def cases_from_ancestor(individual, ancestor):
    """Get cases inherited from an ancestor"""
    return list(filter(lambda case: case[0] == ancestor, individual))


def cases_from_ancestors(individual, interested_ancestors):
    """Get cases from all interested ancestors inherited"""
    return [cases_from_ancestor(individual, ancestor) for ancestor in interested_ancestors]


def percent_cases_from_ancestors_at_specified_generation(results, target_generation, interested_ancestors):
    """Gets information on a per individual basis for how likely it is to select a case from an individual from
       a specified ancestor

    Args:
        results (simulation_results): Results from the simulation function above
        target_generation (int): The generation to inspect
        interested_ancestors ([int]): All ancestor generations of interest

    Returns:
        [[(ancestor_id, float)]]: For each individual, get the ancestor of interest and the likeliness to select a case from
                                  said ancestor
    """
    interested_ancestors.sort()
    targeted_population = results[target_generation]
    num_cases = len(targeted_population[0])
    ### Gathers all cases from specified ancestors
    ancestory_cases_per_individual = [cases_from_ancestors(individual, interested_ancestors) for individual in targeted_population]
    ancestory_rates = []
    for individual_ancestory in ancestory_cases_per_individual:
        ancestory_rates_per_ancestor = []
        for ancestory_cases_id in range(len(individual_ancestory)):
            ### Appends the average rate of being from a specified ancestor
            ancestory_rates_per_ancestor.append((interested_ancestors[ancestory_cases_id], len(individual_ancestory[ancestory_cases_id]) / num_cases))
        ancestory_rates.append(ancestory_rates_per_ancestor)

    return ancestory_rates


def percent_cases_above_ancestor(results, target_generation, targeted_ancestors):
    """Calculates the probability of any case for an individual being above a specified ancestor target

    Args:
        results (simulation_results): Results from the simulation function above
        target_generation (int): The generation to inspect
        targeted_ancestors (int): The ancestor generation to look above (NOTE: target_generation == targeted_ancestors is valid, should always
                                  return downsample rate of simulation)

    Returns:
        float: Percent chance that any case's error is inherited from an ancestor above the target_ancestor genereration
    """
    ### List constructor is to get all ancestors between targeted ancestors and the target generation
    percents_from_ancestor_to_target = percent_cases_from_ancestors_at_specified_generation(results, target_generation, [i for i in range(targeted_ancestors, target_generation + 1)])
    ### To make lists available to insert specified ancestory probabilities into
    percentages_in_nontarget_ancestory_range = [[] for _ in range(len(percents_from_ancestor_to_target[0]))]
    for individual in percents_from_ancestor_to_target:
        for ancestor_id in range(len(individual)):
            percentages_in_nontarget_ancestory_range[ancestor_id].append(individual[ancestor_id][1])
    

    average_ancestor_chance_from_target_ancestor_to_target_generation = [sum(ancestor) / len(ancestor) for ancestor in percentages_in_nontarget_ancestory_range]
    return 1 - sum(average_ancestor_chance_from_target_ancestor_to_target_generation)


def phylogeny_threshold_estimation(d, n):
    """Calculates the threshold estimation from developed formula. NOTE: n should be target_generation - targeted_ancestors + 1
       as the formula is inclusive and the function is exculsive."""
    assert(n >= 0)
    return 1.0 - ((1.0 - d) ** n)


def main():
    """Controls the flow of the program"""
    number_of_generations = 100
    population_size = 100
    number_of_cases = 100
    downsample_rate = 0.1 # 10%

    ### These are simulation results for some arbitrary GP run
    simulation_results = simualte_evolution(number_of_generations, population_size, number_of_cases, downsample_rate)

    target_gen = 100
    target_ancestors = 99

    print(percent_cases_above_ancestor(simulation_results, target_generation=target_gen, targeted_ancestors=target_ancestors))
    print(1 - phylogeny_threshold_estimation(downsample_rate, (target_gen - target_ancestors + 1)))


if __name__ == "__main__":
    main()