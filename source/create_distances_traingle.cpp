#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// ekvivalentno five_two_adic_distance
double compute_distance(std::string& s1, std::string& s2){
    double distance = 0.0;
    int min_len = std::min(s1.size(), s2.size());

    for (int i = 0; i < min_len; i += 3){

        if(s1[i] != s2[i]){
            distance += 1.0;
        }
        else if(s1[i+1] != s2[i+1]){
            distance += 0.2;
        }
        else if(s1[i+2] != s2[i+2]){
            if( (s1[i+2] == 'A' && s2[i+2] == 'G') 
            || (s1[i+2] == 'G' && s2[i+2] == 'A')
            || (s1[i+2] == 'C' && s2[i+2] == 'T')
            || (s1[i+2] == 'T' && s2[i+2] == 'C') ){
                distance += 0.02;
            }
            else{
                distance += 0.04;
            }
        }
    }

    // uzmi u obzir i duzinu sekvenci
    // povecaj distancu za broj tripleta
    double size_diff = (double)s1.size() - (double)s2.size();
    distance += std::abs(size_diff) / 3.0;
    return distance;
}

int main(int argc, char **argv){

    if(argc < 2){
        std::cerr << "no destination file" << std::endl;
        return 1;
    }

    std::ifstream myfile("../data/sequences");

    std::string myline;
    std::vector<std::string> sequences;
    if ( myfile.is_open() ) {

        while ( myfile ) {
            std::getline (myfile, myline);
            sequences.push_back(myline);
        }
    }
    myfile.close();

    sequences.pop_back();
    int n = sequences.size();

    std::ofstream outFile(argv[1]);
    for (int i = 0; i < n-1; i++){
        
        for(int j = i+1; j < n; j++){
            double distance = compute_distance(sequences[i], sequences[j]);
            outFile << distance << ", ";
        }
        outFile << std::endl;
    }
    outFile.close();

    return 0;
}