//
// Created by Benedikt Völker on 11/04/2017.
//

#include "Würfelergebnis.h"


Würfelergebnis::Würfelergebnis() {
    ergebnisse = new vector<int>;
}

void Würfelergebnis::sortieren(){
    sort(ergebnisse->begin(), ergebnisse->end(), std::greater<int>());
}

void Würfelergebnis::würfeln(){
    ergebnisse->push_back(m.natuerlicheZahlinN(6));
}

void Würfelergebnis::würfeln(int n) {
    for(int i = 0; i < n; i++){
        würfeln();
    }
}

void Würfelergebnis::loeschen() {
    ergebnisse->clear();
}

void Würfelergebnis::ergebnisVorgeben(int a, int b, int c) {
    ergebnisse->clear();
    ergebnisse->push_back(a);
    ergebnisse->push_back(b);
    ergebnisse->push_back(c);
}

