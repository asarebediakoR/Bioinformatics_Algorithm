
//Pairwise Alignment

#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <string>
#include <algorithm>
#include <climits>
#include <cstdlib>

//Reading sequences in a FASTA file
 std::vector< std::string> readFastaFile(const  std::string& filename) {
     std::vector< std::string> sequences;
     std::ifstream file(filename);
     std::string line, sequence;
    while(!file.is_open())
    {
         std::cerr << "Invaild file, please check the file format: " <<  std::endl;
         std::string fileName;
         std::cin >> fileName;
        file.open(fileName);
    }
    while (getline(file, line)) {
         std::cout << line <<  std::endl;
        if (line.empty())
            continue;
        if (line[0] == '>') {
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                sequence.clear();
            }
        } else {
            sequence += line;
        }
    }
    if (!sequence.empty()) {
        sequences.push_back(sequence);
    }
    sequences.push_back(sequence);
    file.close();
    return sequences;
}

//Global sequence alignment with affine gap penalties
void globalSequenceAlignment(
    const  std::string& sequence1, const std::string& sequence2,
    int matchScore, int mismatchScore, int gapOpenPenalty, int gapExtensionPenalty
) {
    int n = sequence1.length();
    int m = sequence2.length();

    //Initialize matrices for Dynamic Programming (DP)
     std::vector< std::vector<int>> dp(n + 1,  std::vector<int>(m + 1));
     std::vector< std::vector<int>> gapOpen(n + 1,  std::vector<int>(m + 1));
     std::vector< std::vector<int>> gapExtend(n + 1,  std::vector<int>(m + 1));

    //Initialize matrices with appropriate values based on affine gap penalties
     std::vector< std::vector<int>> gap_open_matrix(n,  std::vector<int>(m, 0));
     std::vector< std::vector<int>> gap_extension_matrix(n,  std::vector<int>(m, 0));

   //Initialize the 1st row of the gap_open_matrix
   for (int i = 1; i < n; ++i) {
    gap_open_matrix[i][0] = gapOpenPenalty + (i - 1) * gapExtensionPenalty;
   }

    //Initialize the 1st column of the gap_open_matrix
   for (int j = 1; j < m; ++j) {
    gap_open_matrix[0][j] = gapOpenPenalty + (j - 1) * gapExtensionPenalty;
   }

    //Initialize the row (sequence1) of the gap_extension_matrix with a large negative value
   for (int i = 1; i < n; ++i) {
    gap_extension_matrix[i][0] = INT_MIN;
   }

    //Initialize the column (sequence2) of the gap_extension_matrix with a large negative value
   for (int j = 1; j < m; ++j) {
    gap_extension_matrix[0][j] = INT_MIN;
   }

    //Filling DP matrices
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            //Calculate match/mismatch score
            int match = dp[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchScore);

            //Calculate gap penalties
            int gap1 = dp[i - 1][j] + gapOpen[i - 1][j] + gapExtensionPenalty;
            int gap2 = dp[i][j - 1] + gapOpen[i][j - 1] + gapExtensionPenalty;

            //maximum scores
            dp[i][j] = std::max({match, gap1, gap2});

            //Update gap extension matrices
            gapOpen[i][j] = std::max(gapOpen[i - 1][j] + gapExtensionPenalty, dp[i - 1][j] + gapOpenPenalty);
            gapExtend[i][j] =  std::max(gapExtend[i][j - 1] + gapExtensionPenalty, dp[i][j - 1] + gapOpenPenalty);
        }
    }

    // Backtracking to find the alignment
     std::string alignedSeq1, alignedSeq2;
    int i = n, j = m;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchScore)) {
            alignedSeq1 += sequence1[i - 1];
            alignedSeq2 += sequence2[j - 1];
            i--;
            j--;
        } else if (i > 0 && dp[i][j] == gapExtend[i][j]) {
            alignedSeq1 += sequence1[i - 1];
            alignedSeq2 += '-';
            i--;
        } else {
            alignedSeq1 += '-';
            alignedSeq2 += sequence2[j - 1];
            j--;
        }
    }

    reverse(alignedSeq1.begin(), alignedSeq1.end());
    reverse(alignedSeq2.begin(), alignedSeq2.end());

    //Print the alignment
     std::cout << "Sequence 1: " << alignedSeq1 <<  std::endl;
     std::cout << "Sequence 2: " << alignedSeq2 <<  std::endl;
}


int gapPenalty; 

std::vector<std::vector<int>> readPAM100Matrix(const std::string &pam100file) {
    std::string fileName = pam100file;
    std::ifstream file(fileName);
    while (!file.is_open()) {
        std::cout << "Error: Unable to open the PAM100 matrix file, kindly check filename." << std::endl;
        // std::cin >> fileName;
    }
    std::string val;
    std::vector<std::vector<int>> pam100(20, std::vector<int>(20, 0));
    int x = 0;
    int y = 0;
   while(file >> val)
   {
        if(std::isalpha(val[0])) continue;
        pam100[y][x] = std::stoi(val);
    // std::cout <<  pam100[y][x] << std::endl;
        x++;
        if(x >= 20)
        {
            x = 0;
            y++;
        }
   }
    file.close();
    return pam100;
}


// Global sequence alignment with affine gap penalties
void globalAlignment(const std::string &sequence1, const std::string &sequence2, const std::vector<std::vector<int>> &pam100) {
    int n = sequence1.length();
    int m = sequence2.length();

    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1));

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int match = dp[i - 1][j - 1] + pam100[sequence1[i - 1] - 'A'][sequence2[j - 1] - 'A'];
            int gap1 = dp[i - 1][j] - gapPenalty;  // Gaps in sequence1
            int gap2 = dp[i][j - 1] - gapPenalty;  // Gaps in sequence2
            dp[i][j] = std::max({0, match, gap1, gap2});
        }
    }

    // Backtracking to find the alignment
    std::string alignedSeq1, alignedSeq2;
    int i = n, j = m;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + pam100[sequence1[i - 1] - 'A'][sequence2[j - 1] - 'A']) {
            alignedSeq1 += sequence1[i - 1];
            alignedSeq2 += sequence2[j - 1];
            i--;
            j--;
        } else if (i > 0 && dp[i][j] == dp[i - 1][j] - gapPenalty) {
            alignedSeq1 += sequence1[i - 1];
            alignedSeq2 += '-';
            i--;
        } else {
            alignedSeq1 += '-';
            alignedSeq2 += sequence2[j - 1];
            j--;
        }
    }

    std::reverse(alignedSeq1.begin(), alignedSeq1.end());
    std::reverse(alignedSeq2.begin(), alignedSeq2.end());

    // Print the alignment
    std::cout << "Sequence 1: " << alignedSeq1 << std::endl;
    std::cout << "Sequence 2: " << alignedSeq2 << std::endl;
}



int main()
 {
     std::string filename;
     char sequenceType;
     std::string pam100File = "pam2.txt";  // Replace pam2 matrix file
    
     std::cout << "Enter the name of the input FASTA file: ";
     std::cin >> filename;

     std::cout << "Are the sequences in the file nucleotide (N) or protein (P)? ";
     std::cin >> sequenceType;

     int matchScore, mismatchScore, gapOpenPenalty, gapExtensionPenalty;
     std::vector< std::string> sequences = readFastaFile(filename);
     if (sequences.size() <= 2)
     {
         std::cerr << "Error: The FASTA file must contain at least two sequences." << std::endl;
         return 1;
     }

     if (sequenceType == 'P' || sequenceType == 'p')
     {
       // For Protein sequences, ask the user for scoring parameters
         std::cout << "Enter gapPenalty score: ";
         std::cin >> gapPenalty;
         
        std::string pam100File;

        std::cout << "Enter the name of the PAM100 matrix file: ";
        std::cin >> pam100File;

        std::vector<std::vector<int>> pam100 = readPAM100Matrix(pam100File);
    
        globalAlignment(sequences[0],sequences[1],pam100);
     } else {
         // For nucleotide sequences, ask the user for scoring parameters
         std::cout << "Enter match score: ";
         std::cin >> matchScore;
         std::cout << "Enter mismatch score: ";
         std::cin >> mismatchScore;
         std::cout << "Enter gap opening penalty: ";
         std::cin >> gapOpenPenalty;
         std::cout << "Enter gap extension penalty: ";
         std::cin >> gapExtensionPenalty;
         globalSequenceAlignment(
         sequences[0], sequences[1],
         matchScore, mismatchScore, gapOpenPenalty, gapExtensionPenalty
     );
     }
     
     return 0;
 }

