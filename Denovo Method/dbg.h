#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>

using namespace std;

vector<string> makeShortRead(string myDna);
map<string, vector<string>> makeEdgeforDBG(string Dna, vector<string> SR);
map<string, int> makeDegreeforDBG(string Dna, vector<string> SR, map<string, vector<string>> node);
vector<string> EulerianGraph(map<string, vector<string>> node, map<string, int> degree, string Dna);
vector<string> findEulerianPath(map<string, vector<string>> node, string start);
double CheckAccuracy(string myDna, string reseqDna);
