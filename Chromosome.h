//
// Created by Andrei on 06.05.2023.
//

#ifndef ALGORITMGENETIC_CHROMOSOME_H
#define ALGORITMGENETIC_CHROMOSOME_H
#include<vector>
#include<iostream>
using namespace std;

class Chromosome
{
private:
    string chromosomeValue;
    int length;
    float selectionProbability;

public:
    Chromosome(int length)
    {
        this->length = length;

        int aux;

        for(int i = 0; i < this->length; ++i)
        {
            aux = std::rand() % 2;
            chromosomeValue += std::to_string(aux);
        }
    }

    Chromosome()
    {

    }

    Chromosome& operator=(const Chromosome& aux)
    {
        if(this != &aux)
        {
            chromosomeValue = aux.chromosomeValue;
            length = aux.length;
        }
        return *this;
    }

    Chromosome(const Chromosome& other)
            : chromosomeValue(other.chromosomeValue), length(other.length),
              selectionProbability(other.selectionProbability) {}

    void setChromosomeValue(const string &chromosomeValue) {
        Chromosome::chromosomeValue = chromosomeValue;
    }

    const string &getChromosomeValue() const {
        return chromosomeValue;
    }

    float getSelectionProbability() const {
        return selectionProbability;
    }

    void setSelectionProbability(float selectionProbability) {
        Chromosome::selectionProbability = selectionProbability;
    }

};


#endif //ALGORITMGENETIC_CHROMOSOME_H
