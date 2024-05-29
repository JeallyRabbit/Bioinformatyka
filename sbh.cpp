// sbh_algorithm.cpp
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <ctime>
#include <deque>

using namespace std;

// Data structure to represent an oligonucleotide and its connections
struct OligoNode {
    string sequence;
    unordered_set<OligoNode*> edges;
    OligoNode(string seq) : sequence(seq) {}
};

// Function to generate oligonucleotides from a DNA sequence
vector<string> generateOligonucleotides(const string& dna_sequence, int k) {
    vector<string> oligonucleotides;
    for (size_t i = 0; i <= dna_sequence.length() - k; ++i) {
        oligonucleotides.push_back(dna_sequence.substr(i, k));
    }
    return oligonucleotides;
}

// Function to generate the graph from the oligonucleotide sequences
unordered_map<string, OligoNode*> generateGraph(const vector<string>& oligonucleotides) {
    unordered_map<string, OligoNode*> graph;
    for (const auto& oligo : oligonucleotides) {
        if (graph.find(oligo) == graph.end()) {
            graph[oligo] = new OligoNode(oligo);
        }
    }
    return graph;
}

// Function to add edges between oligonucleotides
void addEdges(unordered_map<string, OligoNode*>& graph, int k) {
    for (auto& item : graph) {
        const string& seq = item.first;
        OligoNode* node = item.second;
        string prefix = seq.substr(0, k - 1);
        for (auto& other_item : graph) {
            const string& other_seq = other_item.first;
            OligoNode* other_node = other_item.second;
            string suffix = other_seq.substr(1, k - 1);
            if (prefix == suffix) {
                node->edges.insert(other_node);
            }
        }
    }
}

// Function to simulate errors (positive and negative)
vector<string> simulateErrors(const vector<string>& oligonucleotides, int n, bool positive, bool negative) {
    vector<string> erroneous;
    unordered_set<string> seen;
    for (const auto& oligo : oligonucleotides) {
        if (negative && rand() % 2 == 0) continue; // Simulate negative error by omitting the oligo
        erroneous.push_back(oligo);
        seen.insert(oligo);
    }
    if (positive) {
        for (int i = 0; i < n / 10; ++i) {
            string random_oligo = oligonucleotides[rand() % oligonucleotides.size()];
            if (seen.find(random_oligo) == seen.end()) {
                erroneous.push_back(random_oligo);
                seen.insert(random_oligo);
            }
        }
    }
    return erroneous;
}

// Function to reconstruct the DNA sequence using Eulerian path approach
string reconstructDNA(unordered_map<string, OligoNode*>& graph, const string& start_oligo, int k) {
    string dna_sequence = start_oligo;
    OligoNode* current = graph[start_oligo];
    unordered_set<OligoNode*> visited;
    deque<OligoNode*> dfs_stack;
    dfs_stack.push_back(current);

    while (!dfs_stack.empty()) {
        OligoNode* node = dfs_stack.back();
        if (node->edges.empty()) {
            dfs_stack.pop_back();
            if (!dfs_stack.empty()) {
                dna_sequence += dfs_stack.back()->sequence.back();
            }
        }
        else {
            auto neighbor = *node->edges.begin();
            node->edges.erase(neighbor);
            dfs_stack.push_back(neighbor);
            dna_sequence += neighbor->sequence.back();
        }
    }

    return dna_sequence;
}

string generateDNASequence(int n) {
    string dna_sequence;
    const char nucleotides[] = { 'A', 'T', 'G', 'C' };
    for (int i = 0; i < n; ++i) {
        dna_sequence += nucleotides[rand() % 4];
    }
    return dna_sequence;
}


int main() {
    srand(time(0)); // Seed for random number generation

    int n = 100; // Length of DNA
    int k = 3; // Length of oligonucleotides

    // Provided DNA sequence for testing
    string dna_sequence = generateDNASequence(n);
    cout << "Generated DNA sequence: " << dna_sequence << endl;

    // Use the first k nucleotides as the starting oligonucleotide
    string start_oligo = dna_sequence.substr(0, k);

    // Generate oligonucleotides from the DNA sequence
    vector<string> oligonucleotides = generateOligonucleotides(dna_sequence, k);
    /*
    cout << "Generated Oligonucleotides:" << endl;
    for (const auto& oligo : oligonucleotides) {
        cout << oligo << " ";
    }
    cout << endl;*/
    

    // Generate the graph
    auto graph = generateGraph(oligonucleotides);
    addEdges(graph, k);

    // Debug graph
    /*cout << "Graph:" << endl;
    for (const auto& item : graph) {
        cout << "Node: " << item.first << " Edges: ";
        for (const auto neighbor : item.second->edges) {
            cout << neighbor->sequence << " ";
        }
        cout << endl;
    }*/
    

    // Simulate errors
    bool positive = true;
    bool negative = true;
    auto erroneous_oligos = simulateErrors(oligonucleotides, n, positive, negative);

    // Generate the graph from erroneous oligos
    auto erroneous_graph = generateGraph(erroneous_oligos);
    addEdges(erroneous_graph, k);

    // Debug erroneous graph
    /*cout << "Erroneous Graph:" << endl;
    for (const auto& item : erroneous_graph) {
        cout << "Node: " << item.first << " Edges: ";
        for (const auto neighbor : item.second->edges) {
            cout << neighbor->sequence << " ";
        }
        cout << endl;
    }*/
    

    // Reconstruct the DNA sequence
    string dna_sequence_reconstructed = reconstructDNA(erroneous_graph, start_oligo, k);
    cout << "Reconstructed DNA sequence: " << dna_sequence_reconstructed << endl;

    return 0;
}
