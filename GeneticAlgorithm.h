//
// Created by Andrei on 06.05.2023.
//

#ifndef ALGORITMGENETIC_GENETICALGORITHM_H
#define ALGORITMGENETIC_GENETICALGORITHM_H

#include <string>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include "Chromosome.h"
#include <chrono>

using namespace std;

class GeneticAlgorithm
{
private:
    int populationSize; //dimensiunea populatiei
    int x, y; //domeniul de definitie al functiei al carei maxim vrem sa-l calculam
    int a, b, c; //coeficientii polinomului
    int p; //precizia
    double crossoverProbability, mutationProbability; //probabilitatile de recombinare/mutatie
    int generations;//numar etape
    int length; //lungimea unu cromozom;
    double discreteStep; //pas discretizare;
    vector<double> selectionInterval;
    vector<Chromosome> population;
    vector<Chromosome> nextPopulation;
    vector<int> toBeCrossedOver;
public:
    GeneticAlgorithm(int populationSize, int x, int y, int a, int b, int c, int p, double crossoverProbability,
                     double mutationProbability, int generations) : populationSize(populationSize), x(x), y(y), a(a),
                                                                   b(b), c(c), p(p),
                                                                   crossoverProbability(crossoverProbability),
                                                                   mutationProbability(mutationProbability),
                                                                   generations(generations) {}

    GeneticAlgorithm(){}

    friend ifstream& operator>>(ifstream& f, GeneticAlgorithm& aux)
    {
        f>>aux.populationSize>>aux.x>>aux.y>>aux.a>>aux.b>>aux.c>>aux.p
            >>aux.crossoverProbability>>aux.mutationProbability>>aux.generations;

        aux.length = aux.chromosomeLength();
        aux.discreteStep = aux.discrete_Step();
        //cout<<aux.length;
        return f;
    }

    //functia care imi calculeaza lungimea cromozomului; din curs:
    int chromosomeLength()
    {
        double result = log2((b - a) * pow(10, p));
        return ceil(result);
    }

    double discrete_Step()
    {
        double result = (b - a) /(pow(2,length));
        return result;
    }

    //functia care imi transforma numerele din binar in baza 10
    double getValue(Chromosome aux)
    {
        bitset<64> binary(aux.getChromosomeValue());
        double lowerBoundCounter = binary.to_ulong();

        return(a + lowerBoundCounter * discreteStep);
    }

    //functia de fitness care in cazul nostru pur si simplu calculeaza
    //valoarea functiei pentru un anumit cromozom
    double fitness(Chromosome aux)
    {
        double x = getValue(aux);
        double result = a * pow(x, 2) + (b * x) + (c);
        return result;
    }

    void generatePopulation(ofstream& g)
    {
        population.resize(populationSize);
        for(int i = 0; i < populationSize; ++i)
        {
            population[i] = *new Chromosome(this->length);
        }
        g<<"Initial Population: ";
    }

    void printCurrentPopulationToFile(ofstream& g)
    {
        Chromosome aux;
        for(int i = 0; i < population.size(); ++i)
        {
            g<<"\n"<<"Chromosome "<<i + 1<< ": ";

            aux = population[i];

            g<<"x = "<<getValue(aux)<<" f = "<< fitness(aux);
        }
    }

    void printSelectedPopulationToFile(ofstream& g)
    {
        Chromosome aux;
        for(int i = 0; i < nextPopulation.size(); ++i)
        {
            g<<"\n"<<"Chromosome "<<i + 1<< ": ";

            aux = nextPopulation[i];

            g<<"x = "<<getValue(aux)<<" f = "<< fitness(aux);
        }
    }

    vector<double> getFitnessOfAllChromosomes()
    {
        vector<double> fitnessOfAllChromosomes;

        for(int i = 0; i < populationSize; ++i)
            fitnessOfAllChromosomes.push_back(fitness(population[i]));

        return fitnessOfAllChromosomes;
    }

    void getProbability(ofstream& g, bool firstGeneration)
    {
        vector<double> fitnessOfAllChromosomes = getFitnessOfAllChromosomes();
        double total = 0;

        for(int i = 0; i < populationSize; ++i)
            total += fitnessOfAllChromosomes[i];

        for(int i = 0; i < populationSize; ++i)
            population[i].setSelectionProbability(fitness(population[i]) / total);

        if(firstGeneration)
        {
            g<<"\n\nSelection Probabilities: ";

            for(int i = 0; i < populationSize; ++i)
                g<<"\nChromosome "<<i + 1<<": "<<population[i].getSelectionProbability();
        }
    }

    void computeProbabilityIntervals(ofstream& g, bool firstGeneration)
    {
        vector<double> fitnessOfAllChromosomes = getFitnessOfAllChromosomes();
        double curent = 0.0;

        for(int i = 0; i < populationSize; ++i)
        {
            selectionInterval.push_back(curent);
            curent += population[i].getSelectionProbability();
        }

        selectionInterval.push_back(1);

        if(firstGeneration)
        {
            g<<"\n\nSelection Intervals:";
            for(int i = 0; i <= populationSize; ++i)
                g<<"\n"<<selectionInterval[i];
        }

    }

    double get_random()
    {
        std::mt19937_64 rng;
        uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
        rng.seed(ss);
        std::uniform_real_distribution<double> unif(0, 1);
        double randomNumber = unif(rng);
        return randomNumber;
    }

    double get_random(int a, int b)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(a, b);
        return dis(gen);
    }

    int binarySearch(double candidate, const std::vector<double>& selectionInterval)
    {
        int low = 0;
        int high = selectionInterval.size() - 1;

        while (low <= high)
        {
            int mid = low + (high - low) / 2;

            if (selectionInterval[mid] < candidate)
                low = mid + 1;
            else if (selectionInterval[mid - 1] >= candidate)
                high = mid - 1;
            else
                return mid;
        }

        return -1; // Candidate not found (handle this case in your code accordingly)
    }

    void rouletteStrategy(ofstream& g, bool firstGeneration)
    {
        double candidate;
        int indice;

        if(firstGeneration)
        {
            g<<"\n\nChosen population:";
        }

        for(int i = 1; i <= populationSize; ++i)
        {
            candidate = get_random();

            indice = binarySearch(candidate, selectionInterval);

            if(firstGeneration)
                   g<<"\nu = "<<candidate<<" we select chromosome "<<indice;

            nextPopulation.push_back(population[indice - 1]);
        }

//        double candidate1, candidate2;
//        int indice1, indice2;
//
//        if(firstGeneration)
//        {
//            g<<"\n\nChosen population:";
//        }
//
//        for(int i = 1; i <= populationSize; ++i)
//        {
//            indice1 = 1;
//            indice2 = 1;
//
//            candidate1 = get_random();
//            candidate2 = get_random();
//
//            while(!(selectionInterval[indice1 - 1] <= candidate1 and
//                    candidate1 <= selectionInterval[indice1])) indice1++;
//
//            while(!(selectionInterval[indice2 - 1] <= candidate2 and
//                    candidate2 <= selectionInterval[indice2])) indice2++;
//
//            if(fitness(population[indice1 - 1]) > fitness(population[indice2 - 1]))
//            {
//                if(firstGeneration)
//                    g<<"\nu = "<<candidate1<<" we select chromosome "<<indice1;
//
//                nextPopulation.push_back(population[indice1 - 1]);
//            }
//
//            else
//            {
//                if(firstGeneration)
//                    g<<"\nu = "<<candidate2<<" we select chromosome "<<indice2;
//
//                nextPopulation.push_back(population[indice2 - 1]);
//            }
//
//        }
    }

    //choose random number from [0,1] for each chromosome
    //if random number < crossoverProbability => added to the list of chromosomes
    //that will be crossed over;
    void crossOverCandidates(ofstream& g, bool firstGeneration)
    {
        double odds;

        if(firstGeneration)
            g<<"\n\nCrossover odds: "<<crossoverProbability;

        for(int i = 0; i < populationSize; ++i)
        {
            odds = get_random();

            if(firstGeneration)
            {
                g << "\nChromosome " << i + 1 << ": " <<
                  nextPopulation[i].getChromosomeValue() <<
                  " u = " << odds;

                if (odds < crossoverProbability)
                    g << " < " << crossoverProbability << " => participates";
            }

            if (odds < crossoverProbability)
                toBeCrossedOver.push_back(i + 1);
        }
        g<<"\n";
    }

    void swapPartOfTheGenome(int start, Chromosome& c1, Chromosome& c2)
    {
        string suffix1 = c1.getChromosomeValue().substr(c1.getChromosomeValue().length() - start);
        string suffix2 = c2.getChromosomeValue().substr(c2.getChromosomeValue().length() - start);

        string aux1 = c1.getChromosomeValue(), aux2 = c2.getChromosomeValue();

        aux1.replace(c1.getChromosomeValue().length() - start, start, suffix2);
        aux2.replace(c2.getChromosomeValue().length() - start, start, suffix1);

        c1.setChromosomeValue(aux1);
        c2.setChromosomeValue(aux2);

        Chromosome c;
        c = c1;
        c1 = c2;
        c2 = c;

    }

    void crossOver(ofstream& g, bool firstGeneration)
    {
        int first, second, randomPoint;
        while(toBeCrossedOver.size() >= 2)
        {
            if(!firstGeneration)
            {
                randomPoint = get_random(0, length - 1);

                first = get_random(0, toBeCrossedOver.size() - 1);
                int element1 = toBeCrossedOver[first];
                toBeCrossedOver.erase(this->toBeCrossedOver.begin() + first);

                second = get_random(0, toBeCrossedOver.size() - 1);
                int element2 = toBeCrossedOver[second];
                toBeCrossedOver.erase(this->toBeCrossedOver.begin() + second);

                swapPartOfTheGenome(length - randomPoint, nextPopulation[element1],
                                    nextPopulation[element2]);
            }

            if(firstGeneration)
            {
                randomPoint = get_random(0, length - 1);

                first = get_random(0, toBeCrossedOver.size() - 1);
                int element1 = toBeCrossedOver[first];
                toBeCrossedOver.erase(this->toBeCrossedOver.begin() + first);

                g<<"\nCrossover between Chromosome "<<element1;

                second = get_random(0, toBeCrossedOver.size() - 1);
                int element2 = toBeCrossedOver[second];
                toBeCrossedOver.erase(this->toBeCrossedOver.begin() + second);

//                while(first == second)
//                    second = get_random(0, toBeCrossedOver.size() - 2);

                g<<" and Chromosome "<<element2<<":\n"<<
                nextPopulation[element1].getChromosomeValue()<< " || "<<
                nextPopulation[element2].getChromosomeValue()<<
                " swapped genome from the point: " << randomPoint << "\n";


                swapPartOfTheGenome(length - randomPoint, nextPopulation[element1],
                                    nextPopulation[element2]);

                g<<"Result:\n"<<nextPopulation[element1].getChromosomeValue() <<
                 " || "<<nextPopulation[element2].getChromosomeValue();

            }
        }

        g<<"\n";
    }

    void mutation(ofstream& g, bool firstGeneration)
    {
        double odds;

        if(firstGeneration)
        {
            g<<"\n\nMutation probability: "<<mutationProbability<<"\nThe following chromosomes have been changed: \n";

            for(int i = 0; i < nextPopulation.size(); ++i)
            {
                odds = get_random();

                if(odds < mutationProbability)
                {
                    g<<"Chromosome "<< i + 1<<"\n";
                    //we choose a random bit and negate it:
                    int randomBit = get_random(0, length);

                    string chromosomeValue = nextPopulation[i].getChromosomeValue();

                    if(chromosomeValue[randomBit] == '0')
                        chromosomeValue[randomBit] = '1';

                    else
                        chromosomeValue[randomBit] = '0';

                    nextPopulation[i].setChromosomeValue(chromosomeValue);
                }
            }
        }

        else
        {
            for(int i = 0; i < populationSize; ++i)
            {
                odds = get_random();

                if(odds < mutationProbability)
                {
                    //we choose a random bit and negate it:
                    int randomBit = get_random(0, length);

                    string chromosomeValue = nextPopulation[i].getChromosomeValue();

                    if(chromosomeValue[randomBit] == '0')
                        chromosomeValue[randomBit] = '1';

                    else
                        chromosomeValue[randomBit] = '0';

                    nextPopulation[i].setChromosomeValue(chromosomeValue);
                }
            }
        }
    }

    int getFittestIndice()
    {
        int fittestIndice = 0;

        for(int i = 0; i < population.size(); ++i)
            if(fitness(population[fittestIndice]) < fitness(population[i]))
                fittestIndice = i;

        return fittestIndice;
    }

    void elitistCriteria()
    {
        int fittestIndice = 0, leastFitIndice = 0;

        for(int i = 0; i < population.size(); ++i)
            if(fitness(population[fittestIndice]) < fitness(population[i]))
                fittestIndice = i;

        for(int i = 0; i < nextPopulation.size(); ++i)
            if(fitness(nextPopulation[leastFitIndice]) > fitness(nextPopulation[i]))
                leastFitIndice = i;

        nextPopulation[leastFitIndice] = population[fittestIndice];
    }

    double computeAverageFitness()
    {
        double total = 0;

        for(int i = 0; i < nextPopulation.size(); ++i)
            total += fitness(nextPopulation[i]);

        return (total / nextPopulation.size());
    }

    void prepareForNextGeneration()
    {
        for(int i = 0; i < population.size(); ++i)
            population[i] = nextPopulation[i];
    }

    int getGenerations()
    {
        return this->generations;
    }

    const vector<Chromosome> &getNextPopulation() const {
        return nextPopulation;
    }

    const vector<Chromosome> &getPopulation() const{
        return population;
    }
};

#endif //ALGORITMGENETIC_GENETICALGORITHM_H
