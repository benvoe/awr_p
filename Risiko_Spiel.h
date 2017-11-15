//
// Created by Benedikt Völker on 13/04/2017.
//

#ifndef AWR_PRAKTIKUM_RISIKO_SPIEL_H
#define AWR_PRAKTIKUM_RISIKO_SPIEL_H

#include "Wuerfel.h"
#include "MathTools.h"

class Risiko_Spiel {
private:
    Wuerfel* angreifer;
    Wuerfel* verteidiger;

    int _angreifer_siege, _verteigiger_siege, _k, _l, _n;

    void auswerten();
    void initialisieren();
    void ausgeben();

public:

    Risiko_Spiel();

    double spielen ( int k, int l, int n = 1, bool ausgabe = true);
};


#endif //AWR_PRAKTIKUM_RISIKO_SPIEL_H
