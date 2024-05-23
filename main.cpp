#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <unordered_set>

using namespace std;

// O(E)
vector<int> podziel(const string &linia, const string &sep = " ")
{
    vector<int> lista;
    int start = 0;
    int koniec = linia.find(sep);
    while (koniec != -1)
    {
        lista.push_back(stoi(linia.substr(start, koniec - start)));
        start = koniec + sep.size();
        koniec = linia.find(sep, start);
    }
    lista.push_back(stoi(linia.substr(start, koniec - start)));
    return lista;
}

vector<string> readFasta(string fileName)
{
    std::fstream plik;
    vector<string> sekwencje = {};
    plik.open(fileName, ios_base::in);
    string line = "";
    string sekwencja = "";

    if (plik.good())
    {
        getline(plik, line);
        while (getline(plik, line))
        {
            if (line == "" || line[0] == '\n')
            {
                continue;
            }
            else if (line[0] == '>')
            {
                sekwencje.push_back(sekwencja);
                sekwencja = "";
            }
            else
            {
                sekwencja += line;
            }
        }
        sekwencje.push_back(sekwencja);
        plik.close();
        return sekwencje;
    }
    else
    {
        cerr << "Blad otwierania pliku podczas wczytywania: " << errno;
        return {};
    }
}

vector<vector<int>> readQual(string fileName)
{
    std::fstream plik;
    vector<vector<int>> jakosci = {};
    plik.open(fileName, ios_base::in);
    string line = "";
    string sekwencja = "";

    if (plik.good())
    {
        getline(plik, line);

        while (getline(plik, line))
        {
            if (line == "" || line[0] == '\n')
            {
                continue;
            }
            else if (line[0] == '>')
            {
                sekwencja.pop_back();
                jakosci.push_back(podziel(sekwencja));
                sekwencja = "";
            }
            else
            {
                sekwencja += (line + " ");
            }
        }
        sekwencja.pop_back();
        jakosci.push_back(podziel(sekwencja));
        plik.close();
        return jakosci;
    }
    else
    {
        cerr << "Blad otwierania pliku podczas wczytywania: " << errno;
        return {};
    }
}

void removeLowQuals(vector<string> *sekwencje, vector<vector<int>> *quals, vector<vector<int>> *pozycje, int minQual)
{
    for (int i = 0; i < (*quals).size(); i++)
    {
        for (int j = (*quals)[i].size() - 1; j >= 0; j--)
        {
            if ((*quals)[i][j] < minQual)
            {
                (*quals)[i].erase((*quals)[i].begin() + j);
                (*sekwencje)[i].erase((*sekwencje)[i].begin() + j);
                (*pozycje)[i].erase((*pozycje)[i].begin() + j);
            }
        }
    }
}

typedef struct spectrumElem
{
    int nrSekwencji = 0;
    int nrPozyzji = 0;
    string oligo = "";
    int stopien = 0;
} spectrumElem;

bool isEmpty(vector<vector<bool>> matrix)
{
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = i; j < matrix[i].size(); j++)
        {
            if (matrix[i][j] || matrix[j][i])
                return false;
        }
    }
    return true;
}

vector<int> verticeDegreeVector(vector<vector<bool>> matrix)
{
    vector<int> degrees(matrix.size(), 0);
    for (int i = 0; i < matrix.size(); i++)
    {
        degrees[i] += count(matrix[i].begin(), matrix[i].end(), true);
        for (int j = i; j < matrix[i].size(); j++)
        {
            if (matrix[j][i])
                degrees[i]++;
        }
    }
    return degrees;
}

double graphDenisty(vector<vector<bool>> matrix)
{
    double suma = 0.0;
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = i; j < matrix[i].size(); j++)
        {
            if (matrix[j][i])
                suma++;
        }
    }

    return suma / (double)matrix.size();
}

// https://par.nsf.gov/servlets/purl/10175443
vector<spectrumElem> CharikarGreedyPeel(vector<vector<bool>> G, vector<spectrumElem> spectrum)
{
    auto densestG = G;
    auto densestSpectrum = spectrum;
    auto H = G;
    auto tmpSpectrum = spectrum;
    unordered_set<int> etykietySekwencji;

    // O(V)
    for (auto &&item : spectrum)
    {
        etykietySekwencji.insert(item.nrSekwencji);
    }

    int iloscSekwencji = etykietySekwencji.size();
    vector<int> iloscOligoZSekwencji(iloscSekwencji, 0);

    // O(V)
    for (int i = 0; i < spectrum.size(); i++)
    {
        iloscOligoZSekwencji[spectrum[i].nrSekwencji]++;
    }

    // O(V)
    while (H.size() > 0)
    {
        // O(V^2)
        // Find the vertex u ∈ H with minimum degH(u)
        vector<int> tmpVerticeDegrees = verticeDegreeVector(H);
        auto minVertex = min_element(tmpVerticeDegrees.begin(), tmpVerticeDegrees.end());
        int minVertexIndex = distance(tmpVerticeDegrees.begin(), minVertex);

        // sprawdz czy to nie ostatni wierzcholek z ktores sekwencji
        if (iloscOligoZSekwencji[tmpSpectrum[minVertexIndex].nrSekwencji] == 2)
        {
            break;
        }

        // O(V)
        //  Remove u and all its adjacent edges uv from H
        for (int i = 0; i < H.size(); i++)
        {
            H[i].erase(H[i].begin() + minVertexIndex);
        }
        H.erase(H.begin() + minVertexIndex);
        iloscOligoZSekwencji[tmpSpectrum[minVertexIndex].nrSekwencji]--;
        tmpSpectrum.erase(tmpSpectrum.begin() + minVertexIndex);
        if (graphDenisty(H) > graphDenisty(densestG))
        {
            densestG = H;
            densestSpectrum = tmpSpectrum;
        }
    }

    for (auto &&row : densestG)
    {
        for (auto &&item : row)
        {
            if (item)
            {
                cout << "1 ";
            }
            else
            {
                cout << "0 ";
            }
        }
        cout << endl;
    }

    return densestSpectrum;
}

int main()
{
    string fastaFileName = "";
    cout << "Nazwa pliku Fasta:\t";
    cin >> fastaFileName;

    string qualFileName = "";
    cout << "Nazwa pliku Qual:\t";
    cin >> qualFileName;

    int minQual = 0;
    cout << "Min qual [0-40]:\t";
    cin >> minQual;

    int frameWidth = 0;
    cout << "Szerokosc okna [4-9]:\t";
    cin >> frameWidth;

    // Parametry do testow
    // fastaFileName = "instancja5.fasta";
    // qualFileName = "instancja5.qual";
    // minQual = 0;
    // frameWidth = 4;

    // Wczytanie sekwencji oraz wartosci qual
    auto rawSekwencje = readFasta(fastaFileName);
    auto rawQuals = readQual(qualFileName);

    // tworzenie kopii sekwencji wartosci qual
    auto sekwencje = rawSekwencje;
    auto quals = rawQuals;

    // dodawanie sztucznego elementu na pierwszej pozycji aby numeracja rozpoczynała sie od 1
    for (auto &&sekwencja : sekwencje)
    {
        sekwencja.insert(sekwencja.begin(), 'X');
    }
    for (auto &&qual : quals)
    {
        qual.insert(qual.begin(), 999);
    }

    // wypelnianie vectora zawierajacego oryginalne pozycje nukleotydow
    vector<vector<int>> pozycje((int)sekwencje.size());
    for (int i = 0; i < sekwencje.size(); i++)
    {
        for (int j = 0; j < sekwencje[i].size(); j++)
        {
            pozycje[i].push_back(j);
        }
    }

    // uzuwanie z sekwencji i listy qual elementow o warosci qual nizszej niz prog minimalny
    removeLowQuals(&sekwencje, &quals, &pozycje, minQual);

    // tworzenie listy ktora zawiera wszystkie podciagi
    vector<spectrumElem> spectrum;
    for (int i = 0; i < sekwencje.size(); i++)
    {
        if (sekwencje[i].size() < frameWidth)
            continue;
        for (int j = 1; j < sekwencje[i].size() - frameWidth + 1; j++)
        {
            spectrumElem tmpElem = {i, pozycje[i][j], sekwencje[i].substr(j, frameWidth), 0};
            spectrum.push_back(tmpElem);
        }
    }

    if (spectrum.size() < sekwencje.size())
    {
        cout << "Nie udalo stworzyc sie grafu z uzyciem podanych sekwencji i parametrow." << endl;
        return -1;
    }

    vector<vector<bool>> matrix(spectrum.size(), vector<bool>(spectrum.size(), false));
    for (int i = 0; i < spectrum.size(); i++)
    {
        for (int j = i; j < spectrum.size(); j++)
        {
            if (spectrum[i].nrSekwencji != spectrum[j].nrSekwencji && spectrum[i].oligo == spectrum[j].oligo && abs(spectrum[i].nrPozyzji - spectrum[j].nrPozyzji) <= frameWidth * 10)
            {
                matrix[i][j] = true;
                matrix[j][i] = true;
            }
        }
    }

    // O(V^3)
    auto densestSubgraph = CharikarGreedyPeel(matrix, spectrum);

    // O(N log N)
    std::sort(densestSubgraph.begin(), densestSubgraph.end(), [](spectrumElem const &a, spectrumElem const &b)
              { return a.nrSekwencji < b.nrSekwencji; });

    // for (int i = 0; i < densestSubgraph.size(); i++)
    // {
    //     cout << densestSubgraph[i].nrSekwencji + 1 << " " << densestSubgraph[i].nrPozyzji << " " << densestSubgraph[i].oligo << endl;
    // }

    vector<vector<spectrumElem>> podgrupy(sekwencje.size());
    vector<spectrumElem> wynik(sekwencje.size());
    for (auto &&item : densestSubgraph)
    {
        podgrupy[item.nrSekwencji].push_back(item);
    }

    // sprobuj znalezc klike
    for (int i = 0; i < podgrupy[0].size(); i++)
    {
        wynik.clear();
        wynik.push_back(podgrupy[0][i]);
        for (int j = 1; j < podgrupy.size(); j++)
        {
            for (int k = 0; k < podgrupy[j].size(); k++)
            {
                if (podgrupy[0][i].oligo == podgrupy[j][k].oligo && abs(podgrupy[0][i].nrPozyzji - podgrupy[j][k].nrPozyzji) <= frameWidth * 10)
                // if (podgrupy[0][i].oligo == podgrupy[j][k].oligo)
                {
                    wynik.push_back(podgrupy[j][k]);
                    break;
                }
            }
        }
        if (wynik.size() == sekwencje.size())
        {
            break;
        }
    }

    // jesli nie ma kliki sprobuj znalezc strukture podobna do kliki
    if (wynik.size() != sekwencje.size())
    {
        for (int i = 0; i < podgrupy[0].size(); i++)
        {
            wynik.clear();
            wynik.push_back(podgrupy[0][i]);
            for (int j = 1; j < podgrupy.size(); j++)
            {
                for (int k = 0; k < podgrupy[j].size(); k++)
                {
                    // if (podgrupy[0][i].oligo == podgrupy[j][k].oligo && abs(podgrupy[0][i].nrPozyzji - podgrupy[j][k].nrPozyzji) <= frameWidth * 10)
                    if (podgrupy[0][i].oligo == podgrupy[j][k].oligo)
                    {
                        wynik.push_back(podgrupy[j][k]);
                        break;
                    }
                }
            }
            if (wynik.size() == sekwencje.size())
            {
                break;
            }
        }
    }

    // wypisz wynik
    if (wynik.size() == sekwencje.size())
    {
        cout << "Wynik: " << endl;
        for (int i = 0; i < wynik.size(); i++)
        {
            cout << wynik[i].nrSekwencji + 1 << " " << wynik[i].nrPozyzji << " " << wynik[i].oligo << endl;
        }
    }
    else
    {
        cout << "Nie udalo sie znalezc motywu." << endl;
    }

    return 0;
}
