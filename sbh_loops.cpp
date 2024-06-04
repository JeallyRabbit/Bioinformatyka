#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <time.h>
#include <fstream>

using namespace std;

// Funkcja obliczająca odległość Levenshteina między dwoma ciągami
int levenshteinDistance(const string& s1, const string& s2) {
    int m = s1.size();
    int n = s2.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (i == 0)
                dp[i][j] = j;
            else if (j == 0)
                dp[i][j] = i;
            else if (s1[i - 1] == s2[j - 1])
                dp[i][j] = dp[i - 1][j - 1];
            else
                dp[i][j] = 1 + min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]});
        }
    }
    return dp[m][n];
}

// funkcja generująca ciąg DNA o długości n
string generateDNA(int n) {
    string dna = "";
    string bases = "ATCG";
    for (int i = 0; i < n; ++i) {
        int randIndex = rand() % 4; // losowy indeks do wyboru z bases
        dna += bases[randIndex];
    }
    return dna;
}

// funkcja dzieląca ciąg DNA na fragmenty o długości k
vector<string> cutDNA(string dna, int k) {
    vector<string> fragments;
    for (int i = 0; i <= dna.length() - k; ++i) {
        fragments.push_back(dna.substr(i, k));
    }
    return fragments;
}

// funkcja dodająca błędy pozytywne i negatywne
string addErrors(string dna, double p_positive, double p_negative) {
    // Dodawanie błędów pozytywnych
    for (int i = 0; i < dna.length(); ++i) {
        if (rand() / double(RAND_MAX) < p_positive) {
            int randIndex = rand() % 4; // losowy indeks do wyboru z bases
            dna.insert(i, 1, "ATCG"[randIndex]);
        }
    }
    // Dodawanie błędów negatywnych
    for (int i = dna.length() - 1; i >= 0; --i) {
        if (rand() / double(RAND_MAX) < p_negative) {
            dna.erase(i, 1);
        }
    }
    return dna;
}

// funkcja znajdująca scieżkę przechodzącą przez wszystkie krawędzie w grafie
void findPath(vector<vector<int>>& graph, vector<int>& path, vector<int>& visited, int current) {
    for (int i = 0; i < graph.size(); ++i) {
        if (graph[current][i] > 0 && visited[i] == 0) {
            visited[i] = 1;
            findPath(graph, path, visited, i);
            path.push_back(i);
        }
    }
}

int main() {
    srand(time(0));

    fstream plik;
    plik.open("bioinf.csv");
    plik<<"n;k;bledy_pozotywne;bledy_negatywne;miara_levenshteina\n";

    for(int n=300;n<=700;n+=20)
    {
        for(int k=6;k<=10;k++)
        {
            for(int p_positive_num=0;p_positive_num<40;p_positive_num++)
            {
                for(int p_negative_num=0;p_negative_num<40;p_negative_num++)
                {
                    double p_positive=p_negative_num*0.01;
                    double p_negative=p_negative_num*0.01;
                    string originalDNA = generateDNA(n);

                    // Krok 2: Dzielenie DNA na fragmenty
                    vector<string> fragments = cutDNA(originalDNA, k);

                    // Krok 3: Dodawanie błędów
                    for (int i = 0; i < fragments.size(); ++i) {
                        fragments[i] = addErrors(fragments[i], p_positive, p_negative);
                    }

                    // Krok 4: Tworzenie grafu
                    vector<vector<int>> graph(fragments.size(), vector<int>(fragments.size(), 0));
                    for (int i = 0; i < fragments.size(); ++i) {
                        for (int j = 0; j < fragments.size(); ++j) {
                            if (i != j) {
                                int commonPrefix = 0;
                                while (commonPrefix < k && fragments[i][commonPrefix] == fragments[j][commonPrefix]) {
                                    commonPrefix++;
                                }
                                graph[i][j] = k - commonPrefix;
                            }
                        }
                    }

                    // Krok 5: Znajdowanie ścieżki przez graf
                    vector<int> path;
                    vector<int> visited(fragments.size(), 0);
                    visited[0] = 1;
                    findPath(graph, path, visited, 0);
                    path.push_back(0); // Dodaj pierwszy wierzchołek, aby zamknąć pętlę

                    // Krok 6: Odtwarzanie DNA z pierwszych liter wierzchołków
                    string reconstructedDNA = "";
                    for (int i = path.size() - 1; i >= 0; --i) {
                        reconstructedDNA += fragments[path[i]].substr(0, 1);
                    }

                    // Krok 7: Odtwarzanie komplementarnej nici DNA
                    string complementaryDNA = "";
                    for (char base : reconstructedDNA) {
                        switch (base) {
                            case 'A':
                                complementaryDNA += 'T';
                                break;
                            case 'T':
                                complementaryDNA += 'A';
                                break;
                            case 'C':
                                complementaryDNA += 'G';
                                break;
                            case 'G':
                                complementaryDNA += 'C';
                                break;
                        }
                    }
                    plik<<n<<";"<<k<<";"<<p_positive<<";"<<p_negative<<";"<<levenshteinDistance(reconstructedDNA, originalDNA)<<"\n";
                    //cout<<n<<" "<<k<<" "<<p_positive<<" "<<p_negative<<" "<<levenshteinDistance(reconstructedDNA, originalDNA)<<endl;

                }
            }
        }
    }
    // Parametry
    int n = 30; // długość DNA
    int k = 3;  // długość fragmentów
    double p_positive = 0.1; // prawdopodobieństwo błędu pozytywnego
    double p_negative = 0.1; // prawdopodobieństwo błędu negatywnego

    // Krok 1: Generowanie DNA
    
    

    // Krok 8: Wypisywanie wyników
    /*cout << "Pierwotne DNA: " << originalDNA << endl;
    cout << "Zrekonstruowane DNA: " << reconstructedDNA << endl;
    cout << "Komplementarne do zrekonstruowanego: " << complementaryDNA << endl;

     cout << "Odległość Levenshteina między zrekonstruowanym a pierwotnym DNA: "
         << levenshteinDistance(reconstructedDNA, originalDNA) << endl;*/
    

    return 0;
}
