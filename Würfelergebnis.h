//
// Created by Benedikt Völker on 11/04/2017.
//

#ifndef AWR_PRAKTIKUM_WÜRFELERGEBNIS_H
#define AWR_PRAKTIKUM_WÜRFELERGEBNIS_H

#include <iostream>
#include <algorithm>
#include "CMyMath.h"

class Würfelergebnis {
public:
    vector<int>* ergebnisse;
    CMyMath m = CMyMath();

    Würfelergebnis ();

    void sortieren();
    void würfeln();
    void würfeln(int n);
    void loeschen();
    void ergebnisVorgeben(int a, int b, int c = 0);
};

#endif //AWR_PRAKTIKUM_WÜRFELERGEBNIS_H
