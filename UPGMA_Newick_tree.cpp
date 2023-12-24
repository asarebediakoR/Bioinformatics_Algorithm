//UPGMA tree/Newick tree

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>


using namespace std;

//Function to calculate the distance between two sequences
double calculateDistance(const string& seq1, const string& seq2) {
    int mismatchCount = 0;
    int length = min(seq1.length(), seq2.length());
    for (int i = 0; i < length; i++) {
        if (seq1[i] != seq2[i]) {
            mismatchCount++;
        }
    }
    return static_cast<double>(mismatchCount)/length;
}


//Function to generate a consensus sequence
string generateConsensus(const vector<string>& sequences) {
    int length = sequences[0].length();
    string consensus;
    for (int i = 0; i < length; i++) 
    {
        char consensusChar = sequences[0][i];
        bool mismatch = false;
        for (const string& seq : sequences) {
            if (seq[i] != consensusChar) {
                mismatch = true;
                break;
            }
        }
        if (mismatch) {
            consensus += '*';  //disagreement
        } else {
            consensus += consensusChar;  //agreement
        }
    }
    return consensus;
}

//Creating the Newick tree with  calculated distances
string createNewickTree(const vector<string>& sequences, const vector<string>& sequenceNames) {
    if (sequences.size() == 1) {
        return sequenceNames[0];
    }

    int minDistIndex1 = 0, minDistIndex2 = 1;
    double minDistance = calculateDistance(sequences[0], sequences[1]);

    for (int i = 0; i < sequences.size(); i++) {
        for (int j = i + 1; j < sequences.size(); j++) {
            double dist = calculateDistance(sequences[i], sequences[j]);
            if (dist < minDistance) {
                minDistance = dist;
                minDistIndex1 = i;
                minDistIndex2 = j;
            }
        }
    }

    string newick1 = createNewickTree({sequences[minDistIndex1]}, {sequenceNames[minDistIndex1]});
    string newick2 = createNewickTree({sequences[minDistIndex2]}, {sequenceNames[minDistIndex2]});

    double roundedDistance = round(minDistance * 100.0) / 100.0;

    return "(" + newick1 + ":" + to_string(roundedDistance) + "," + newick2 + ":" + to_string(roundedDistance) + ")";
}

int main() {
    string fileName;
    cout << "Please enter the name of your FASTA sequence file,Kindly check spellings: ";
    cin >> fileName;

    ifstream inputFile(fileName);
    if (!inputFile.is_open()) {
        cout << "Failed to open FASTA file." << endl;
        return 1;
    }

    vector<string> sequences;
    vector<string> sequenceNames;
    string line, sequence;
    while (getline(inputFile, line)) {
        if (line[0] == '>') {
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                sequence.clear();
            }
            sequenceNames.push_back(line.substr(1));
        } else {
            sequence += line;
        }
    }
    if (!sequence.empty()) {
        sequences.push_back(sequence);
    }


    //Calculate distances, build tree, and generate consensus
    string newickTree = createNewickTree(sequences, sequenceNames);
    string consensus = generateConsensus(sequences);

    //Print the results
    cout << "The Consensus Sequence is: " << consensus << endl;
    cout << "Cordinates for Newick Tree: " << newickTree << endl;
    cout << "The * in seq represents position all seq agree in residue" << endl;
    return 0;
}

