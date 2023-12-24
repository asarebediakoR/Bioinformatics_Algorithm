//Multiple Seq. Alignment

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

//Defining the structure of sequence (FASTA format)
struct Sequence {
    std::string header;
    std::string data;
};

//Reading sequences from a FASTA file
std::vector<Sequence> readFASTA(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open the FASTA file,Kindly check the file." << std::endl;
        exit(1);
    }

    std::vector<Sequence> sequences;
    Sequence currentSequence;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue; 

        if (line[0] == '>') { 
            if (!currentSequence.header.empty()) {
                sequences.push_back(currentSequence);
                currentSequence.header.clear();
                currentSequence.data.clear();
            }
            currentSequence.header = line.substr(1); //Removing '>' in the FASTA format
        } else {
            currentSequence.data += line;
        }
    }

    if (!currentSequence.header.empty()) {
        sequences.push_back(currentSequence);
    }

    return sequences;
}

//Function to calculate the consensus sequence
std::string calculateConsensus(const std::vector<Sequence>& sequences) {
    int sequenceLength = sequences[0].data.length();
    std::string consensus(sequenceLength, ' '); //Initialize consensus with spaces

    for (int i = 0; i < sequenceLength; ++i) {
        std::vector<int> counts(26, 0);       //Counting occurrences of each letter 

        for (const auto& sequence : sequences) {
            char residue = sequence.data[i];
            if (residue >= 'A' && residue <= 'Z') {
                counts[residue - 'A']++;
            }
        }

        int maxCount = 0;
        char consensusResidue = ' ';

        for (int j = 0; j < 26; ++j) {
            if (counts[j] > maxCount) {
                maxCount = counts[j];
                consensusResidue = static_cast<char>('A' + j);
            }
        }
            //Agreement and disagreement statements
        if (maxCount == sequences.size()) {
            consensus[i] = consensusResidue;
        } else {
            consensus[i] = '*'; 
        }
    }

    return consensus;
}

int main() {

    //Ask user to input the filename
    std::cout << "Enter the name of the FASTA file: ";
    std::string filename;
    std::cin >> filename;

    //Reading sequences from a FASTA file 
    std::vector<Sequence> sequences = readFASTA(filename);

    //Calculate the consensus sequence
    std::string consensus = calculateConsensus(sequences);

    //Printout the consensus sequence
    std::cout << " The Consensus Sequence is:" << std::endl;
    std::cout << consensus << std::endl;
    std::cout << "The * in seqeunce reps position all seqs agree for residue" << std:: endl;

    return 0;
}
