//
// Created by Benedikt Völker on 09.05.17.
//

#ifndef AWR_PRAKTIKUM_WUERFEL_H
#define AWR_PRAKTIKUM_WUERFEL_H

#include <iostream>
#include <algorithm>
#include "CMyMath.h"

class Wuerfel {
public:
    std::vector<int>* historie;
    CMyMath m = CMyMath();

    Wuerfel();

    void sortieren();
    void würfeln();
    void würfeln(int n);
    void loeschen();
};


#endif //AWR_PRAKTIKUM_WUERFEL_H
