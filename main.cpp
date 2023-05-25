#include <iostream>
#include <fstream>
#include "GeneticAlgorithm.h"

ifstream f("/Users/Andrei/CLionProjects/AlgoritmGenetic/input");
ofstream g("/Users/Andrei/CLionProjects/AlgoritmGenetic/Evolutie.txt");

int main()
{
    GeneticAlgorithm alg;
    // we get the inputs from a file and compute the length of the chromosomes
    // which is done in the override of ">>" operator:
    f>>alg;
    f.close();
    //Step 1: we randomly generate the first generation and write it in the file:
    alg.generatePopulation(g);
    alg.printCurrentPopulationToFile(g);
    bool firstGeneration = true;
    for(int i = 0; i < alg.getGenerations(); ++i)
    {
        //Step 2: we compute the probability of selection for each chromosome:
        alg.getProbability(g, firstGeneration);
        //Step 3: compute selection intervals:
        alg.computeProbabilityIntervals(g, firstGeneration);
        //Step 4: roulette strategy:
        //which consists in choosing a random number between [0,1]
        //and then look at what probability interval that number belongs
        //that chromosome goes into the next generation
        alg.rouletteStrategy(g, firstGeneration);
        //show population after selecting it:
        if(firstGeneration)
        {
            g<<"\n\nAfter selection: ";
            alg.printSelectedPopulationToFile(g);
        }
        //Step 5:
        //we select which population shall crossover:
        alg.crossOverCandidates(g, firstGeneration);
        //Step 6:
        //the crossover:
        alg.crossOver(g, firstGeneration);
        //print population after crossover for the first generation:
        if(firstGeneration)
            alg.printSelectedPopulationToFile(g);
        //Step 8:
        //mutation:
        alg.mutation(g, firstGeneration);
        //print population after mutation for the first generation:
        if(firstGeneration)
        {
            g<<"\nAfter mutation:";
            alg.printSelectedPopulationToFile(g);
        }
        //Step 9:
        //elitist criteria:
        alg.elitistCriteria();
        if(firstGeneration)
        {
            g<<"\n\nElitist criteria:";
            alg.printSelectedPopulationToFile(g);
        }
        //Step 10:
        //We print the maximum element evolution:
        if(firstGeneration)
        {
            g<<"\n\nEvolution of Maximum Value:\n";
        }
        int fittestIndice = alg.getFittestIndice();
        g<<"Maximum of generation "<<i + 1<<" is for x = "<<setprecision(10)<<alg.getValue(alg.getPopulation()[fittestIndice]);
        g<<"\nValue of f("<<alg.getValue(alg.getPopulation()[fittestIndice])<<") = "<<alg.fitness(alg.getPopulation()[fittestIndice]);
        g<<"\nAverage for generation "<<i + 1<<" for f(x): "<<setprecision(10)<<alg.computeAverageFitness();
        //Preparing for next generation: nextGeneration => current generation;
        alg.prepareForNextGeneration();
        firstGeneration = false;
    }
    return 0;
}
