//
// Created by Benedikt Völker on 09.05.17.
//

#include "Wuerfel.h"


Wuerfel::Wuerfel() {
    historie = new vector<int>;
}

void Wuerfel::sortieren(){
    sort(historie->begin(), historie->end(), std::greater<int>());
}

void Wuerfel::würfeln(){
    historie->push_back(m.natuerlicheZahlinN(6));
}

void Wuerfel::würfeln(int n) {
    for(int i = 0; i < n; i++){
        würfeln();
    }
}

void Wuerfel::loeschen() {
    historie->clear();
}

