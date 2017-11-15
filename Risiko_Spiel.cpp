//
// Created by Benedikt Völker on 13/04/2017.
//

#include "Risiko_Spiel.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Risiko_Spiel::Risiko_Spiel() {
    angreifer = new Wuerfel();
    verteidiger = new Wuerfel();

    initialisieren();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Risiko_Spiel::spielen(int k, int l, int n, bool ausgabe) {
    initialisieren();

    _n = n;
    _k = k;
    _l = l;

    for(int i = 0; i < n; i++){
        angreifer->loeschen();
        verteidiger->loeschen();

        angreifer->würfeln(k);
        verteidiger->würfeln(l);
        auswerten();
    }
    if(ausgabe){
        ausgeben();
    }

    double diff = _angreifer_siege / double(_n) - 0.371656;
    if(diff < 0) diff *= -1;

    return diff / n;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Risiko_Spiel::auswerten () {
    int anz = 0;
    if(_k > _l)
        anz = _l;
    else
        anz = _k;

    angreifer->sortieren();
    verteidiger->sortieren();

    int angreifer_teilsiege = 0;
    int verteidiger_teilsiege = 0;

    for(int i = 0; i < anz; i++){
        if(angreifer->historie->at(i) > verteidiger->historie->at(i)){
            angreifer_teilsiege++;
        }
        else {
            verteidiger_teilsiege++;
        }
    }

    if( angreifer_teilsiege > verteidiger_teilsiege ){
        _angreifer_siege++;
    }
    else if ( angreifer_teilsiege < verteidiger_teilsiege ) {
        _verteigiger_siege++;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Risiko_Spiel::initialisieren() {
    angreifer->loeschen();
    verteidiger->loeschen();

    _angreifer_siege = 0;
    _verteigiger_siege = 0;
    _k = 0;
    _l = 0;
    _n = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Risiko_Spiel::ausgeben() {
    cout << "Auswertung des Spiels mit k = " << _k << " und l = " << _l << endl;
    cout << "Der Angreifer gewinnt zu " << _angreifer_siege / double(_n) * 100 << " %" << endl;
    cout << "Der Verteidiger gewinnt zu " << _verteigiger_siege / double(_n) * 100 << " %" << endl;
    cout << "Zum Unentschieden kam es zu " << (_n - _angreifer_siege - _verteigiger_siege ) / double(_n) * 100 << " %" << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////